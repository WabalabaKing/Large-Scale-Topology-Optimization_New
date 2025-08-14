// filter_sens_omp_atomic.c
// Sensitivity filtering (no volume weighting) with OpenMP atomics.
// Files read from CWD: drow.dat, dcol.dat, dval.dat
// Math: df/dx_i = sum_j [ H_{j i} / sum_k H_{j k} ] * (df/dx_tilde_j)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE (1<<20)  // ~1M triplets per block; tune as needed
#endif

static inline int inb(int x, int n){ return (unsigned)x < (unsigned)n; }

// Count nnz by scanning triplet text files in lockstep.
static long long count_nnz_from_files(void){
    FILE *fr=fopen("drow.dat","r");
    FILE *fc=fopen("dcol.dat","r");
    FILE *fv=fopen("dval.dat","r");
    if(!fr||!fc||!fv){
        fprintf(stderr,"ERROR: cannot open drow.dat/dcol.dat/dval.dat: %s\n", strerror(errno));
        if(fr) fclose(fr); if(fc) fclose(fc); if(fv) fclose(fv);
        return -1;
    }
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

// Public API:
// df_dxtilde: in, size ne  (∂f/∂x̃_j)
// df_dx:      out, size ne  (∂f/∂x_i), zeroed inside
// ne:         number of elements
// nnz_total:  total nnz in H; if <=0, auto-count
// q:          exponent on H weights (usually 1.0; q==1 is fast-pathed)
void filter_sensitivities_buffered_mt(const double *df_dxtilde,
                                     double *df_dx,
                                     int ne,
                                     long long nnz_total,
                                     double *q_ptr)
{

    

    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }

    if (nnz_total <= 0) {
        nnz_total = count_nnz_from_files();
        if (nnz_total <= 0) {
            fprintf(stderr,"ERROR: could not determine nnz_total from files.\n");
            exit(EXIT_FAILURE);
        }
    }

    // Allocate accumulators and I/O blocks
    double *row_sum = (double*)calloc((size_t)ne, sizeof(double));
    if(!row_sum){ fprintf(stderr,"alloc row_sum failed\n"); exit(EXIT_FAILURE); }
    for (int i=0;i<ne;++i) df_dx[i] = 0.0;

    int *drow = (int*)malloc(BLOCK_SIZE * sizeof(int));
    int *dcol = (int*)malloc(BLOCK_SIZE * sizeof(int));
    double *dval = (double*)malloc(BLOCK_SIZE * sizeof(double));
    if(!drow||!dcol||!dval){ fprintf(stderr,"alloc triplet blocks failed\n"); exit(EXIT_FAILURE); }

    const int use_pow = (q != 1.0);

    // ---------------- PASS 1: build donor row sums ----------------
    {
        FILE *fr=fopen("drow.dat","r");
        FILE *fc=fopen("dcol.dat","r");
        FILE *fv=fopen("dval.dat","r");
        if(!fr||!fc||!fv){
            fprintf(stderr,"ERROR: cannot open drow.dat/dcol.dat/dval.dat: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }

        long long tot=0;
        while (tot < nnz_total) {
            int n = (int)((nnz_total - tot) < BLOCK_SIZE ? (nnz_total - tot) : BLOCK_SIZE);
            for (int i=0;i<n;++i) {
                if (fscanf(fr,"%d",&drow[i])!=1 ||
                    fscanf(fc,"%d",&dcol[i])!=1 ||
                    fscanf(fv,"%lf",&dval[i])!=1) {
                    fprintf(stderr,"triplet read error pass1 at %lld\n", tot+i);
                    exit(EXIT_FAILURE);
                }
            }

            // Parallel over the buffered block; atomic add to row_sum[j]
            #pragma omp parallel for schedule(static)
            for (int t=0; t<n; ++t) {
                int j = drow[t]-1;   // donor row (in H)
                int i = dcol[t]-1;   // receiver col
                if (!inb(j,ne) || !inb(i,ne)) continue;
                double w = use_pow ? pow(dval[t], q) : dval[t];
                #pragma omp atomic
                row_sum[j] += w;
            }

            tot += n;
        }
        fclose(fr); fclose(fc); fclose(fv);
    }

    // ---------------- PASS 2: accumulate df_dx ----------------
    {
        FILE *fr=fopen("drow.dat","r");
        FILE *fc=fopen("dcol.dat","r");
        FILE *fv=fopen("dval.dat","r");
        if(!fr||!fc||!fv){
            fprintf(stderr,"ERROR: cannot open drow.dat/dcol.dat/dval.dat: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }

        long long tot=0;
        while (tot < nnz_total) {
            int n = (int)((nnz_total - tot) < BLOCK_SIZE ? (nnz_total - tot) : BLOCK_SIZE);
            for (int i=0;i<n;++i) {
                if (fscanf(fr,"%d",&drow[i])!=1 ||
                    fscanf(fc,"%d",&dcol[i])!=1 ||
                    fscanf(fv,"%lf",&dval[i])!=1) {
                    fprintf(stderr,"triplet read error pass2 at %lld\n", tot+i);
                    exit(EXIT_FAILURE);
                }
            }

            // Parallel over the buffered block; atomic add to df_dx[i]
            #pragma omp parallel for schedule(static)
            for (int t=0; t<n; ++t) {
                int j = drow[t]-1; // donor j
                int i = dcol[t]-1; // design index i
                if (!inb(j,ne) || !inb(i,ne)) continue;
                double denom = row_sum[j];
                if (denom <= 0.0) continue;
                double w = use_pow ? pow(dval[t], q) : dval[t];
                double contrib = (w / denom) * df_dxtilde[j];
                #pragma omp atomic
                df_dx[i] += contrib;
            }

            tot += n;
        }
        fclose(fr); fclose(fc); fclose(fv);
    }

    free(row_sum); free(drow); free(dcol); free(dval);
}
