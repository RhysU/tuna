/**
 * Copyright (C) 2013, 2026 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * An Tuna example demonstrating bare-bones API usage.
 */

#include <stdio.h>
#include <stdlib.h>

#include <tuna.h>

int main(int argc, char* argv[])
{
    // Parse any incoming command line arguments
    const int    niter = argc > 1 ? atof(argv[1]) : 1000  ; // Iteration count?
    const double mA    = argc > 2 ? atof(argv[2]) :   10.0; // Case A mean?
    const double sA    = argc > 3 ? atof(argv[3]) :    1.0; // Case A stddev?
    const double mB    = argc > 4 ? atof(argv[4]) :   10.1; // Case B mean?
    const double sB    = argc > 5 ? atof(argv[5]) :    1.0; // Case B stddev?

    static tuna_site  site;                 // Notice zero initialization
    tuna_stack        stack;                // Stack-based state
    static tuna_chunk chunks[2];            // Notice zero initialization
    tuna_state state = tuna_state_default(); // Used only to simulate chunk timings
    for (int i = 0; i < niter; ++i) {

        // Autotune over the alternatives (simulating chunk-specific costs)
        // To track runtime via TUNA_CLOCK, call tuna_post(&site, chunks) instead
        double cost;
        switch (tuna_pre(&site, &stack, chunks, tuna_countof(chunks))) {
            default: cost = mA + tuna_rand_n01(&state) * sA;
                break;
            case 1:  cost = mB + tuna_rand_n01(&state) * sB;
                break;
        }
        tuna_post_cost(&stack, chunks, cost);

    }

    // Display settings and static memory overhead required for autotuning
    printf("niter=%d, mA=%g, sA=%g, mB=%g, sB=%g, memory=%zd bytes\n",
           niter, mA, sA, mB, sB, sizeof(site) + sizeof(chunks));

    // Display observations from each alternative
    tuna_fprint(stdout, &site, chunks, tuna_countof(chunks), "basic", NULL);

    return EXIT_SUCCESS;
}
