/*
--------------------------------------------------------------------------------
filter_sens_pthread_sharded.c

Two‑pass sensitivity filtering (no volume weighting) using pthreads + atomics,
with:
  • Fast stdio (setvbuf(..., 8 MB))
  • Per‑block weight precompute (w = dval^q)
  • Sharded output accumulation to reduce atomic contention on hot columns

Math
  Forward (density) filter:
      x̃_i = ( Σ_j H_{i j} x_j ) / ( Σ_j H_{i j} )

  Adjoint sensitivity filter:
      df/dx_i = Σ_j [ H_{j i} / ( Σ_k H_{j k} ) ] * (df/dx̃_j)

Files (current directory, text, 1‑based indices)
  drow.dat : row j per line
  dcol.dat : col i per line
  dval.dat : weight H_{j i} per line
  All have equal line counts (nnz).

API
  void filterSensitivity_buffered_mt(const double *SensIn,  // df/dx̃ (size ne)
                                     double *SensOut,       // df/dx  (size ne)
                                     int ne,
                                     long long nnz_total,   // <=0 → auto-count
                                     double q);             // typically 1.0

Build
  gcc -O3 -march=native -ffast-math -pthread -fopenmp \
      filter_sens_pthread_sharded.c -o filter_sens_sharded

Notes
  • Keep SHARDS a power‑of‑two (e.g., 4, 8, 16). Memory = SHARDS * ne * 8 bytes.
  • If row sums also become a hotspot, you can shard row_sum similarly.
--------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

#ifdef _OPENMP
  #include <omp.h>   // allows #pragma omp atomic if available
#endif

// ---- Tunables ---------------------------------------------------------------

// ≈ 32 MB of triplets per block: (r,c,v,w) ~ 8 bytes + 8 + 8 + 8 + 4 + 4 ~ 32B
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 100000000
#endif

// Output shards to reduce contention (must be power-of-two: 4..16 typical)
#ifndef SHARDS
#define SHARDS 8
#endif

// ---- Atomic add macro (OpenMP atomic or GCC atomics fallback) ---------------
#ifdef _OPENMP
  #define ATOMIC_ADD(var, val) do { _Pragma("omp atomic") (var) += (val); } while(0)
#else
  #define ATOMIC_ADD(var, val) __atomic_fetch_add(&(var), (val), __ATOMIC_RELAXED)
#endif

// ---- Thread arguments -------------------------------------------------------
typedef struct {
    int thread_id, num_threads;
    int *drow, *dcol;          // buffered indices (1-based in files)
    double *w;                 // buffered weights (already ^q)
    int block_read, ne;

    // Pass 1:
    double *row_sum;           // size ne: row_sum[j] = Σ_k H_{j k}^q

    // Pass 2 (sharded output):
    const double *SensIn;      // df/dx̃ (size ne)
    double **out_shard;        // SHARDS × ne arrays
} ThreadArgs;

// ---- Small helpers ----------------------------------------------------------
static inline int inb(int x,int n){ return (unsigned)x < (unsigned)n; }

// Count nnz by scanning files once (lockstep), with big stdio buffers
static long long count_nnz_from_files(void){
    FILE *fr=fopen("drow.dat","r");
    FILE *fc=fopen("dcol.dat","r");
    FILE *fv=fopen("dval.dat","r");
    if(!fr||!fc||!fv){
        fprintf(stderr,"ERROR: open drow/dcol/dval: %s\n", strerror(errno));
        if(fr) fclose(fr); if(fc) fclose(fc); if(fv) fclose(fv);
        return -1;
    }
    setvbuf(fr,NULL,_IOFBF,8<<20);
    setvbuf(fc,NULL,_IOFBF,8<<20);
    setvbuf(fv,NULL,_IOFBF,8<<20);

    long long cnt=0; int r,c; double v;
    for(;;){
        int okr=fscanf(fr,"%d",&r);
        int okc=fscanf(fc,"%d",&c);
        int okv=fscanf(fv,"%lf",&v);
        if(okr==1 && okc==1 && okv==1) ++cnt; else break;
    }
    fclose(fr); fclose(fc); fclose(fv);
    return cnt;
}

// ---- Workers ----------------------------------------------------------------

// Pass 1: build donor row sums
static void *worker_build_row_sums(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;

    for (int t = s; t < e; ++t) {
        int j = a->drow[t] - 1;           // donor row
        if (!inb(j,a->ne)) continue;      // 'i' not needed here
        ATOMIC_ADD(a->row_sum[j], a->w[t]);
    }
    return NULL;
}

// Pass 2: accumulate sensitivities with sharded output
static void *worker_accumulate_sens(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;

    for (int t = s; t < e; ++t) {
        int j = a->drow[t] - 1;                   // donor row
        int i = a->dcol[t] - 1;                   // design column
        if (!inb(j,a->ne) || !inb(i,a->ne)) continue;

        double denom = a->row_sum[j];             // Σ_k H_{j k}^q
        if (denom <= 0.0) continue;

        double contrib = (a->w[t] / denom) * a->SensIn[j];

        // Choose shard by low bits (SHARDS must be power-of-two)
        int shard = i & (SHARDS - 1);
        ATOMIC_ADD(a->out_shard[shard][i], contrib);
    }
    return NULL;
}

// ---- Public API -------------------------------------------------------------
void filterSensitivity_buffered_mts(const double *SensIn,  // df/dx̃[j]
                                   double *SensOut,       // df/dx[i] (output)
                                   int ne,
                                   long long nnz_total)   // <=0 → auto-count)              // usually 1.0
{

    double q = 1.0;
    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }

    // #threads from OMP_NUM_THREADS or default 4
    int num_threads = 4;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env && *env) {
        int tmp = atoi(env);
        if (tmp > 0) num_threads = tmp;
    }
    printf("filterSensitivity_buffered_mt: using %d thread(s)\n", num_threads);

    // Determine nnz_total if needed
    if (nnz_total <= 0) {
        nnz_total = count_nnz_from_files();
        if (nnz_total <= 0) {
            fprintf(stderr,"ERROR: could not determine nnz_total from files.\n");
            exit(EXIT_FAILURE);
        }
    }

    printf("Opening filter matrix files for donor row sums (PASS 1)...");
    // ---------------- PASS 1: donor row sums ----------------
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        fprintf(stderr,"ERROR: opening triplet files: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    setvbuf(frow,NULL,_IOFBF,8<<20);
    setvbuf(fcol,NULL,_IOFBF,8<<20);
    setvbuf(fval,NULL,_IOFBF,8<<20);
    printf("done!\n");

    printf("Allocating block buffers and threads...");
    // Block buffers
    int    *drow_block = (int*)   malloc((size_t)BLOCK_SIZE * sizeof(int));
    int    *dcol_block = (int*)   malloc((size_t)BLOCK_SIZE * sizeof(int));
    double *dval_block = (double*)malloc((size_t)BLOCK_SIZE * sizeof(double));
    double *w_block    = (double*)malloc((size_t)BLOCK_SIZE * sizeof(double));
    if (!drow_block || !dcol_block || !dval_block || !w_block) {
        fprintf(stderr,"ERROR: block buffer allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Donor normalization vector
    double *row_sum = (double*)calloc((size_t)ne, sizeof(double));
    if (!row_sum) { fprintf(stderr,"ERROR: row_sum allocation failed\n"); exit(EXIT_FAILURE); }

    // Thread state
    pthread_t  *threads = (pthread_t*) malloc((size_t)num_threads * sizeof(pthread_t));
    ThreadArgs *targs   = (ThreadArgs*)malloc((size_t)num_threads * sizeof(ThreadArgs));
    if (!threads || !targs) { fprintf(stderr,"ERROR: thread allocation failed\n"); exit(EXIT_FAILURE); }
    printf("done \n");

    printf("starting stream pass 1...");
    // Stream pass‑1
    {
        long long total_read = 0;
        while (total_read < nnz_total) {
            long long rem = nnz_total - total_read;
            int n = (int)(rem < BLOCK_SIZE ? rem : (long long)BLOCK_SIZE);

            for (int i=0; i<n; ++i) {
                if (fscanf(frow,"%d",&drow_block[i])!=1 ||
                    fscanf(fcol,"%d",&dcol_block[i])!=1 ||
                    fscanf(fval,"%lf",&dval_block[i])!=1) {
                    fprintf(stderr,"ERROR: triplet read (pass 1) at %lld\n", total_read+i);
                    exit(EXIT_FAILURE);
                }
            }

            // Precompute w = (dval)^q (fast‑paths for q=1,2)
            if (q == 1.0) {
                memcpy(w_block, dval_block, (size_t)n*sizeof(double));
            } else if (q == 2.0) {
                for (int i=0;i<n;++i) w_block[i] = dval_block[i] * dval_block[i];
            } else {
                for (int i=0;i<n;++i) w_block[i] = pow(dval_block[i], q);
            }

            for (int t = 0; t < num_threads; ++t) {
                targs[t] = (ThreadArgs){
                    .thread_id=t, .num_threads=num_threads,
                    .drow=drow_block, .dcol=dcol_block, .w=w_block,
                    .block_read=n, .ne=ne,
                    .row_sum=row_sum, .SensIn=NULL, .out_shard=NULL
                };
                if (pthread_create(&threads[t], NULL, worker_build_row_sums, &targs[t]) != 0) {
                    perror("pthread_create pass1"); exit(EXIT_FAILURE);
                }
            }
            for (int t=0; t<num_threads; ++t) pthread_join(threads[t], NULL);

            total_read += n;
        }
    }
    fclose(frow); fclose(fcol); fclose(fval);

    printf("done \n");
    // ---------------- PASS 2: accumulate sensitivities (sharded) --------------
    // Allocate sharded outputs and zero

    printf("Starting pass 2 with shards...");
    double *out_shard[SHARDS];
    for (int s=0; s<SHARDS; ++s) {
        out_shard[s] = (double*)calloc((size_t)ne, sizeof(double));
        if (!out_shard[s]) { fprintf(stderr,"ERROR: out_shard alloc failed\n"); exit(EXIT_FAILURE); }
    }

    // Reopen files for pass‑2
    frow = fopen("drow.dat", "r");
    fcol = fopen("dcol.dat", "r");
    fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        fprintf(stderr,"ERROR: reopening triplet files (pass 2): %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    setvbuf(frow,NULL,_IOFBF,8<<20);
    setvbuf(fcol,NULL,_IOFBF,8<<20);
    setvbuf(fval,NULL,_IOFBF,8<<20);

    // Stream pass‑2
    {
        long long total_read = 0;
        while (total_read < nnz_total) {
            long long rem = nnz_total - total_read;
            int n = (int)(rem < BLOCK_SIZE ? rem : (long long)BLOCK_SIZE);

            for (int i=0; i<n; ++i) {
                if (fscanf(frow,"%d",&drow_block[i])!=1 ||
                    fscanf(fcol,"%d",&dcol_block[i])!=1 ||
                    fscanf(fval,"%lf",&dval_block[i])!=1) {
                    fprintf(stderr,"ERROR: triplet read (pass 2) at %lld\n", total_read+i);
                    exit(EXIT_FAILURE);
                }
            }

            if (q == 1.0) {
                memcpy(w_block, dval_block, (size_t)n*sizeof(double));
            } else if (q == 2.0) {
                for (int i=0;i<n;++i) w_block[i] = dval_block[i] * dval_block[i];
            } else {
                for (int i=0;i<n;++i) w_block[i] = pow(dval_block[i], q);
            }

            for (int t = 0; t < num_threads; ++t) {
                targs[t] = (ThreadArgs){
                    .thread_id=t, .num_threads=num_threads,
                    .drow=drow_block, .dcol=dcol_block, .w=w_block,
                    .block_read=n, .ne=ne,
                    .row_sum=row_sum, .SensIn=SensIn,
                    .out_shard=out_shard
                };
                if (pthread_create(&threads[t], NULL, worker_accumulate_sens, &targs[t]) != 0) {
                    perror("pthread_create pass2"); exit(EXIT_FAILURE);
                }
            }
            for (int t=0; t<num_threads; ++t) pthread_join(threads[t], NULL);

            total_read += n;
        }
    }
    fclose(frow); fclose(fcol); fclose(fval);

    // Reduce shards → SensOut
    for (int i=0; i<ne; ++i) {
        double sum = 0.0;
        for (int s=0; s<SHARDS; ++s) sum += out_shard[s][i];
        SensOut[i] = sum;
    }

    // Cleanup
    for (int s=0; s<SHARDS; ++s) free(out_shard[s]);
    free(drow_block); free(dcol_block); free(dval_block); free(w_block);
    free(row_sum); free(threads); free(targs);

    printf("done!");
}
