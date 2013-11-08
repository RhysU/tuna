/**
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <tuna/tuna.h>

static const double U01(void) { return rand() / (double) RAND_MAX; }
static const double N01(void) { return tuna_ltqnorm(U01());        }

// TODO This is currently a work-in-progress to discover the correct APIs.
int main(int argc, char *argv[])
{
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    srand((unsigned int) (tp.tv_sec + tp.tv_nsec));

    // Parse and display any incoming command line arguments in a header
    const int    niter = argc > 1 ? atoi(argv[1]) : 100  ; // Iteration count?
    const double mA    = argc > 2 ? atof(argv[2]) :  10.0; // Case A mean?
    const double sA    = argc > 3 ? atof(argv[3]) :   1.0; // Case A stddev?
    const double mB    = argc > 4 ? atof(argv[4]) :  10.1; // Case B mean?
    const double sB    = argc > 5 ? atof(argv[5]) :   1.0; // Case B stddev?
    const int    debug = argc > 6 ? atoi(argv[6]) :   0  ; // Output debugging?
    printf("# niter=%d, mA=%g, sA=%g, mB=%g, sB=%g\n", niter, mA, sA, mB, sB);

    tuna_kernel k[2] = {/*zero fill*/};
    for (int i = 0; i < niter; ++i) {
        // Which branch should be taken this iteration?
        const int ndx = !tuna_stats_cnt(&k[0].stats) ? 0
                      : !tuna_stats_cnt(&k[1].stats) ? 1
                      : tuna_stats_welch1(&k[0].stats,&k[1].stats) < U01();

        // Take the branch, tracking the "elapsed" time.
        double elapsed;
        switch (ndx) {
            case 0: elapsed = mA + N01()*sA; break;
            case 1: elapsed = mB + N01()*sB; break;
        }
        tuna_kernel_obs(k+ndx, elapsed);

        // Output the chosen branch and behavior for debugging purposes
        if (debug) printf("%6d\t%d\t%g\n", i, ndx, elapsed);
    }

    // Summarize results, including any discarded outliers
    for (int i = 0; i < tuna_countof(k); ++i) {
        tuna_stats s = {};
        tuna_kernel_merge(&s, k + i);
        printf("# m%c=%g, s%c=%g, c%c=%zd\n",
               'A' + i, tuna_stats_avg(&s),
               'A' + i, tuna_stats_std(&s),
               'A' + i, tuna_stats_cnt(&s));
    }

    return EXIT_SUCCESS;
}
