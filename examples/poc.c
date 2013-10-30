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

#include "tuna/ltqnorm.h"
#include "tuna/stats.h"

static const double U01(void) { return rand() / (double) RAND_MAX; }
static const double N01(void) { return tuna_ltqnorm(U01());        }

// A minimum viable proof of concept for dynamic A/B-based autotuning
// using the minimum possible machinery from the libtuna library.
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
    const int    b     = argc > 6 ? atof(argv[6]) :   7  ; // How much burn in?
    const int    debug = argc > 7 ? atoi(argv[7]) :   0  ; // Output debugging?
    printf("# niter=%d, mA=%g, sA=%g, mB=%g, sB=%g, b=%d\n",
           niter, mA, sA, mB, sB, b);

    tuna_stats s[2] = {/*zero fill*/};
    for (int i = 0; i < niter; ++i) {
        // Which branch should be taken this iteration?
        const int ndx = tuna_stats_cnt(s+0) < b ? 0
                      : tuna_stats_cnt(s+1) < b ? 1
                      : tuna_stats_welch1(s+0,s+1) > U01();

        // Take the branch, tracking the "elapsed" time.
        double elapsed;
        switch (ndx) {
            case 0: elapsed = mA + N01()*sA; break;
            case 1: elapsed = mB + N01()*sB; break;
        }
        tuna_stats_obs(s+ndx, elapsed);

        // Output the chosen branch and behavior for debugging purposes
        if (debug) printf("%6d\t%d\t%g\n", i, ndx, elapsed);
    }

    // Summarize results
    printf("# mA=%g, sA=%g, cA=%zd\n",
           tuna_stats_avg(s+0), tuna_stats_std(s+0), tuna_stats_cnt(s+0));
    printf("# mB=%g, sB=%g, cB=%zd\n",
           tuna_stats_avg(s+1), tuna_stats_std(s+1), tuna_stats_cnt(s+1));

    return EXIT_SUCCESS;
}
