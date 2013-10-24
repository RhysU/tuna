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

#include "tuna/stats.h"

const double u01(void) { return rand() / (double) RAND_MAX; }

// This is currently a work-in-progress to discover the correct APIs.
int main(int argc, char *argv[])
{
    srand(time(NULL));

    // Parse and display any incoming command line arguments in a header
    const int    niter = argc > 1 ? atoi(argv[1]) : 25  ; // How many iterations?
    const double m0    = argc > 2 ? atof(argv[2]) : 10.0; // Case 1 mean elapsed?
    const double r0    = argc > 3 ? atof(argv[3]) :  1.0; // Case 1 variability?
    const double m1    = argc > 4 ? atof(argv[4]) : 10.1; // Case 2 mean elapsed?
    const double r1    = argc > 5 ? atof(argv[5]) :  1.0; // Case 2 variability?
    const int    b     = argc > 6 ? atof(argv[6]) :  3  ; // How much burn in?
    printf("# niter=%d, m0=%g, r0=%g, m1=%g, r1=%g, b=%d\n",
           niter, m0, r0, m1, r1, b);

    tuna_stats s[2] = {/*zero fill*/};
    for (int i = 0; i < niter; ++i) {
        // Which branch should be taken this iteration?
        const int ndx = tuna_stats_cnt(s+0) < b ? 0
                      : tuna_stats_cnt(s+1) < b ? 1
                      : tuna_stats_welch1(s+0,s+1) > u01();

        // Take the branch, tracking the "elapsed" time.
        double elapsed;
        switch (ndx) {
            case 0: elapsed = m0 + r0*(u01()-0.5); break;
            case 1: elapsed = m1 + r1*(u01()-0.5); break;
        }
        tuna_stats_obs(s+ndx, elapsed);

        // Output the chosen branch and behavior for debugging purposes
        printf("%6d\t%d\t%g\n", i, ndx, elapsed);
    }

    // Summarize results
    printf("# m0=%g, s0=%g, c0=%zd\n",
           tuna_stats_avg(s+0), tuna_stats_std(s+0), tuna_stats_cnt(s+0));
    printf("# m1=%g, s1=%g, c1=%zd\n",
           tuna_stats_avg(s+1), tuna_stats_std(s+1), tuna_stats_cnt(s+1));

    return EXIT_SUCCESS;
}
