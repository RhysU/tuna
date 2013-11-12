/**
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * Dumps iteration-by-iteration algorithm behavior to aid comparisons.  See
 * \ref tuna_algo_default() for how to change the default algorithm at runtime
 * using the process environment.
 */

#include <stdio.h>
#include <stdlib.h>

#include <tuna.h>

static const double N01(void)
{ return tuna_ltqnorm(rand() / (double) RAND_MAX); }

int main(int argc, char *argv[])
{
    // Parse and display any incoming command line arguments in a header
    const int    niter = argc > 1 ? atof(argv[1]) : 1000  ; // Iteration count?
    const double mA    = argc > 2 ? atof(argv[2]) :   10.0; // Case A mean?
    const double sA    = argc > 3 ? atof(argv[3]) :    1.0; // Case A stddev?
    const double mB    = argc > 4 ? atof(argv[4]) :   10.1; // Case B mean?
    const double sB    = argc > 5 ? atof(argv[5]) :    1.0; // Case B stddev?
    const int    best  = mA < mB  ? 0 : 1;                  // Who should win?
    printf("# niter=%d, mA=%g, sA=%g, mB=%g, sB=%g\n", niter, mA, sA, mB, sB);

    static tuna_site   s;
    static tuna_kernel k[2];
    for (int i = 0; i < niter; ++i) {

        // Simulate one iteration of autotuning over alternatives
        double cost;
        switch (tuna_pre(&s, k, tuna_countof(k))) {
            default: cost = mA + N01()*sA;
                     break;
            case 1:  cost = mB + N01()*sB;
                     break;
        }
        tuna_post_cost(&s, k, cost);

        // Output running information for best alternative (includes outliers).
        //
        // An perfect algorithm has the final column be identically one as this
        // indicates the fastest-in-the-limit-of-infinite-data algorithm has
        // been universally preferred.
        tuna_stats o = {};
        tuna_kernel_merge(&o, k + best);
        printf("m%c\t%d\t%12.8g\t%12.8g\t%12zd\t%12.8g\n", 'A' + best, i,
                tuna_stats_avg(&o), tuna_stats_std(&o), tuna_stats_cnt(&o),
                (double) tuna_stats_cnt(&o) / niter);
    }

    return EXIT_SUCCESS;
}
