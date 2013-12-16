/**
 * Copyright (C) 2013 Rhys Ulerich
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

#include <tuna.h>

static tuna_site si;      // Normally si and ks would be inside blockedmm()
static tuna_chunk ks[12]; // but they are global to permit querying them

// Preform a blocked matrix-matrix multiply using autotuned blocking.
static
void
blockedmm(double       c[], // Output C += A*B
          const double a[], // Input A of size N-by-N
          const double b[], // Input B of size N-by-N
          const int log2N)  // Base-2 log of A, B, and C's dimension
{
    // Sanity check incoming arguments
    assert(log2N >=  0);
    assert(log2N <= tuna_countof(ks));

    // Autotune on the power-of-2 block size over range possible given N.
    // Routine (weirdly) took log2N so (trivially) any smaller
    // power-of-2 block size will evenly divide the matrix extents.
    tuna_stack st;
    const int N = 1 << log2N;
    const int B = 1 << tuna_pre(&si, &st, ks, log2N);

    // Multiply logic based upon the documentation (but not the source)
    // from https://code.google.com/p/mm-matrixmultiplicationtool/
    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                const double *sub_a = &a[i + N*k];
                const double *sub_b = &b[k + N*j];
                double       *sub_c = &c[i + N*j];
                for (int l = 0; l < B; ++l) {
                    for (int m = 0; m < B; ++m) {
                        for (int n = 0; n < B; ++n) {
                            sub_c[l + N*m] += sub_a[l + N*n] * sub_b[n + N*m];
                        }
                    }
                }
            }
        }
    }

    // Record autotuning results
    tuna_post(&st, ks);
}

int main(int argc, char *argv[])
{
    // Parse any incoming command line arguments and echo settings
    const int niter = argc > 1 ? atof(argv[1]) : 64; // Iteration count?
    const int log2N = argc > 2 ? atof(argv[2]) :  8; // log2 of matrix size?
    const int N     = 1 << log2N;
    printf("niter=%d, log2N=%d, N=%d\n", niter, log2N, N);

    // Allocate storage
    double* const a = malloc(N*N*sizeof(double));
    double* const b = malloc(N*N*sizeof(double));
    double* const c = malloc(N*N*sizeof(double));

    // On each iteration, fill matrices with U[0,1] data and computed C += A*B
    for (int i = 0; i < niter; ++i) {
        for (int i = 0; i < N*N; ++i) a[i] = rand() / (double) RAND_MAX;
        for (int i = 0; i < N*N; ++i) b[i] = rand() / (double) RAND_MAX;
        for (int i = 0; i < N*N; ++i) c[i] = rand() / (double) RAND_MAX;
        blockedmm(c, a, b, log2N);
    }

    // Deallocate storage
    free(c);
    free(b);
    free(a);

    // Display observations.  Chunks are labeled by log2 of block size.
    tuna_fprint(stdout, &si, ks, tuna_countof(ks), "blockedmm", NULL);

    return EXIT_SUCCESS;
}
