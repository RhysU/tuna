/**
 * Copyright (C) 2013, 2026 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * A Tuna example for dynamically investigating block size trade offs for
 * simple blocked matrix-matrix multiplication.  Rather than tuning over
 * algorithms, here tuning over algorithmic parameters is demonstrated.
 *
 * \warning This problem is chosen because it is a simple, self-contained use
 * case.  However, one should prefer an architecture-specific BLAS
 * implementation over naive code like this example.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
# include <omp.h>
#endif

#include <tuna.h>

// Normally site and chunks might be inside blockedmm() but
// they are global to permit querying them in main().
static tuna_site site;
static const char* labels[] = { "block___1",
                                "block___2",
                                "block___4",
                                "block___8",
                                "block__16",
                                "block__32",
                                "block__64",
                                "block_128"
                              };
static tuna_chunk chunks[tuna_countof(labels)];

// Preform a blocked matrix-matrix multiply using autotuned blocking.
static
void
blockedmm(double*       c, // Output C += A*B
          const double* a, // Input A of size N-by-N
          const double* b, // Input B of size N-by-N
          const int log2N)  // Base-2 log of A, B, and C's dimension
{
    // Routine takes log2N so any smaller power-of-2 block size
    // will trivially partition the matrix into uniform submatrices.
    assert(log2N >= 0);
    const int N = 1 << log2N;
    const int log2B = log2N < tuna_countof(chunks) ? log2N + 1 : tuna_countof(chunks);
    tuna_stack stack;
    int B;
    #pragma omp critical (tuna_blockedmm)
    B = 1 << tuna_pre(&site, &stack, chunks, log2B);

    // Multiply logic based upon the documentation (but not the source)
    // from https://code.google.com/p/mm-matrixmultiplicationtool/
    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                const double* sub_a = &a[i + N * k];
                const double* sub_b = &b[k + N * j];
                double*       sub_c = &c[i + N * j];
                for (int l = 0; l < B; ++l) {
                    for (int m = 0; m < B; ++m) {
                        for (int n = 0; n < B; ++n) {
                            sub_c[l + N * m] += sub_a[l + N * n] * sub_b[n + N * m];
                        }
                    }
                }
            }
        }
    }

    // Update autotuning knowledge
    #pragma omp critical (tuna_blockedmm)
    tuna_post(&stack, chunks);
}

int main(int argc, char* argv[])
{
    // Parse any incoming command line arguments
    const int niter = argc > 1 ? atof(argv[1]) : 256; // Iteration count?
    const int log2N = argc > 2 ? atof(argv[2]) :   7; // log2 of matrix size?
    const int N     = 1 << log2N;

    #pragma omp parallel default(none) firstprivate(niter, log2N, N)
    {
        #pragma omp single
        printf("niter=%d, log2N=%d, N=%d\n", niter, log2N, N);

        // Fill matrices with U[0,1] data and repeatedly compute C += A*B
        double* const a = malloc(N * N * sizeof(double));
        double* const b = malloc(N * N * sizeof(double));
        double* const c = malloc(N * N * sizeof(double));
        #pragma omp for schedule(runtime)
        for (int i = 0; i < niter; ++i) {
            for (int i = 0; i < N * N; ++i) {
                a[i] = rand() / (double) RAND_MAX;
            }
            for (int i = 0; i < N * N; ++i) {
                b[i] = rand() / (double) RAND_MAX;
            }
            for (int i = 0; i < N * N; ++i) {
                c[i] = rand() / (double) RAND_MAX;
            }
            blockedmm(c, a, b, log2N);
        }
        free(c);
        free(b);
        free(a);
    }

    // Display observations
    tuna_fprint(stdout, &site, chunks, tuna_countof(chunks), "blockedmm", labels);

    return EXIT_SUCCESS;
}
