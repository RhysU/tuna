/**
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * An example demonstrating bare-bones API usage.
 */

#include <stdio.h>
#include <stdlib.h>

#include <tuna.h>

static const double N01(void)
{ return tuna_ltqnorm(rand() / (double) RAND_MAX); }

int main(int argc, char *argv[])
{
    // Parse any incoming command line arguments
    const int    niter = argc > 1 ? atof(argv[1]) : 1000  ; // Iteration count?
    const double mA    = argc > 2 ? atof(argv[2]) :   10.0; // Case A mean?
    const double sA    = argc > 3 ? atof(argv[3]) :    1.0; // Case A stddev?
    const double mB    = argc > 4 ? atof(argv[4]) :   10.1; // Case B mean?
    const double sB    = argc > 5 ? atof(argv[5]) :    1.0; // Case B stddev?

    static tuna_site   s;
    static tuna_kernel k[2];
    for (int i = 0; i < niter; ++i) {

        // Autotune over the alternatives (simulating kernel-specific costs)
        // To track runtime via TUNA_CLOCK, call tuna_post(&s, k) instead
        double cost;
        switch (tuna_pre(&s, k, tuna_countof(k))) {
            default: cost = mA + N01()*sA;
                     break;
            case 1:  cost = mB + N01()*sB;
                     break;
        }
        tuna_post_cost(&s, k, cost);

    }

    // Display settings and static memory overhead required for autotuning
    printf("niter=%d, mA=%g, sA=%g, mB=%g, sB=%g, memory=%zd bytes\n",
           niter, mA, sA, mB, sB, sizeof(s) + sizeof(k));

    // Display observations from each alternative
    for (int i = 0; i < tuna_countof(k); ++i) {
        tuna_stats o = {};
        tuna_kernel_merge(&o, k + i);
        printf("m%c=%g, s%c=%g, c%c=%zd\n",
               'A' + i, tuna_stats_avg(&o),
               'A' + i, tuna_stats_std(&o),
               'A' + i, tuna_stats_cnt(&o));
    }

    return EXIT_SUCCESS;
}
