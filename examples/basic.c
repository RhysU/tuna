/**
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <stdio.h>
#include <stdlib.h>

#include <tuna/tuna.h>

static const double N01(void)
{ return tuna_ltqnorm(rand() / (double) RAND_MAX); }

int main(int argc, char *argv[])
{
    // Parse and display any incoming command line arguments in a header
    const int    niter = argc > 1 ? atof(argv[1]) : 100  ; // Iteration count?
    const double mA    = argc > 2 ? atof(argv[2]) :  10.0; // Case A mean?
    const double sA    = argc > 3 ? atof(argv[3]) :   1.0; // Case A stddev?
    const double mB    = argc > 4 ? atof(argv[4]) :  10.1; // Case B mean?
    const double sB    = argc > 5 ? atof(argv[5]) :   1.0; // Case B stddev?
    const int    debug = argc > 6 ? atoi(argv[6]) :   0  ; // Output debugging?
    printf("# niter=%d, mA=%g, sA=%g, mB=%g, sB=%g\n", niter, mA, sA, mB, sB);

    static tuna_state  s;
    static tuna_kernel k[2];
    for (int i = 0; i < niter; ++i) {
        // Which branch should be taken this iteration?
        const int ndx = tuna_pre(&s, k, tuna_countof(k));

        // Take the branch, tracking "elapsed" time in a hypothetical kernel
        double elapsed;
        switch (ndx) {
            default:
            case 0:  elapsed = mA + N01()*sA; break;
            case 1:  elapsed = mB + N01()*sB; break;
        }
        tuna_kernel_obs(k+ndx, elapsed);

        // Output the chosen branch and behavior for debugging purposes
        if (debug) printf("%6d\t%d\t%g\n", i, ndx, elapsed);
    }

    // Display choice summary
    for (int i = 0; i < tuna_countof(k); ++i) {
        tuna_stats o = {};
        tuna_kernel_merge(&o, k + i);
        printf("# m%c=%g, o%c=%g, c%c=%zd\n",
               'A' + i, tuna_stats_avg(&o),
               'A' + i, tuna_stats_std(&o),
               'A' + i, tuna_stats_cnt(&o));
    }

    return EXIT_SUCCESS;
}
