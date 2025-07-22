#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "CalculiX.h"

typedef struct {
    int thread_id;
    int ne0;
    int ne_start, ne_end;
    double *elCentroid;
    double rmin_local;
    int *fnnzassumed;
    int *filternnzElems;

    int *drow;
    int *dcol;
    double *dval;
    int *nnz_count;  // local count array (length = ne_end - ne_start)
    int max_buffer;
} ThreadArgs;

void *filter_thread_worker(void *args_ptr) {
    ThreadArgs *args = (ThreadArgs *)args_ptr;
    int count = 0;

    for (int i = args->ne_start; i < args->ne_end; ++i) {
        int local_nnz = 0;

        double xi = args->elCentroid[3 * i + 0];
        double yi = args->elCentroid[3 * i + 1];
        double zi = args->elCentroid[3 * i + 2];

        for (int j = 0; j < args->ne0; ++j) {
            if (i == j) continue;

            double xj = args->elCentroid[3 * j + 0];
            double yj = args->elCentroid[3 * j + 1];
            double zj = args->elCentroid[3 * j + 2];

            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);

            if (dist <= args->rmin_local) {
                double w = args->rmin_local - dist;

                int idx = count + local_nnz * 2;
                if (idx + 1 >= args->max_buffer) {
                    fprintf(stderr, "ERROR: Thread %d buffer overflow at element %d\n", args->thread_id, i);
                    exit(EXIT_FAILURE);
                }

                args->drow[idx]     = i + 1;  // 1-based indexing
                args->dcol[idx]     = j + 1;
                args->dval[idx]     = w;

                args->drow[idx + 1] = j + 1;
                args->dcol[idx + 1] = i + 1;
                args->dval[idx + 1] = w;

                local_nnz++;
            }
        }

        args->filternnzElems[i] = local_nnz;  // still uses global index
        args->nnz_count[i - args->ne_start] = local_nnz * 2;  // must subtract ne_start!
        count += local_nnz * 2;

        if (local_nnz > *args->fnnzassumed) {
            printf("WARNING: Element %d has %d neighbors. Increase fnnzassumed.\n", i, local_nnz);
            exit(EXIT_FAILURE);
        }
    }

    return NULL;
}

void densityfilterFast_mt(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                          ITG *ne, double *ttime, double *timepar,
                          ITG *mortar, double *rmin, ITG *filternnz,
                          ITG *filternnzElems, ITG itertop, ITG *fnnzassumed)
{
    int num_threads = 1;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) num_threads = atoi(env);
    if (num_threads <= 0) num_threads = 1;

    ITG ne0 = *ne;
    double time = timepar[1];

    // Compute centroids
    double *elCentroid = NULL;
    NNEW(elCentroid, double, 3 * ne0);
    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, elCentroid);

    int elems_per_thread = (ne0 + num_threads - 1) / num_threads;
    int max_neighbors = *fnnzassumed;
    int buffer_per_thread = 2 * max_neighbors * elems_per_thread;
    int total_buffer = buffer_per_thread * num_threads;

    int *global_drow = malloc(total_buffer * sizeof(int));
    int *global_dcol = malloc(total_buffer * sizeof(int));
    double *global_dval = malloc(total_buffer * sizeof(double));
    int *global_nnz_count = calloc(ne0, sizeof(int));

    if (!global_drow || !global_dcol || !global_dval || !global_nnz_count) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = malloc(num_threads * sizeof(ThreadArgs));

    if (!threads || !args) {
        fprintf(stderr, "Thread control allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    for (int t = 0; t < num_threads; ++t) {
        int start = t * elems_per_thread;
        int end = (t + 1) * elems_per_thread;
        if (end > ne0) end = ne0;

        int offset = t * buffer_per_thread;

        args[t] = (ThreadArgs){
            .thread_id = t,
            .ne0 = ne0,
            .ne_start = start,
            .ne_end = end,
            .elCentroid = elCentroid,
            .rmin_local = *rmin,
            .fnnzassumed = fnnzassumed,
            .filternnzElems = filternnzElems,
            .drow = &global_drow[offset],
            .dcol = &global_dcol[offset],
            .dval = &global_dval[offset],
            .nnz_count = &global_nnz_count[start],  // per-thread section
            .max_buffer = buffer_per_thread
        };

        if (pthread_create(&threads[t], NULL, filter_thread_worker, &args[t]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", t);
            exit(EXIT_FAILURE);
        }
    }

    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], NULL);
    }

    // Write to files
    FILE *frow = fopen("drow.dat", "w");
    FILE *fcol = fopen("dcol.dat", "w");
    FILE *fval = fopen("dval.dat", "w");
    FILE *fdnnz = fopen("dnnz.dat", "w");

    int total_nnz = 0;
    for (int i = 0; i < ne0; ++i) {
        fprintf(fdnnz, "%d\n", filternnzElems[i]);
        total_nnz += global_nnz_count[i];
    }

    for (int t = 0; t < num_threads; ++t) {
        int offset = t * buffer_per_thread;
        int start = args[t].ne_start;
        int end = args[t].ne_end;
        int local_idx = 0;

        for (int i = start; i < end; ++i) {
            int row_nnz = global_nnz_count[i];
            for (int k = 0; k < row_nnz; ++k, ++local_idx) {
                fprintf(frow, "%d\n", global_drow[offset + local_idx]);
                fprintf(fcol, "%d\n", global_dcol[offset + local_idx]);
                fprintf(fval, "%.6f\n", global_dval[offset + local_idx]);
            }
        }
    }

    *filternnz = total_nnz;

    fclose(frow); fclose(fcol); fclose(fval); fclose(fdnnz);
    SFREE(elCentroid);
    free(global_drow); free(global_dcol); free(global_dval); free(global_nnz_count);
    free(threads); free(args);
}
