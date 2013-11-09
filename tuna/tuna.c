/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * \copydoc tuna.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tuna.h"

#include <time.h>

// TODO Do something intelligent with clock_getres(2) information

int tuna_pre(tuna_site* st,
             const tuna_kernel* ks,
             const int nk)
{
    // Ensure a zero-initialize st argument produces good behavior by...
    if (!st->al) {
        // ...providing a default algorithm when not set, and
        st->al = tuna_algo_default();

        if (!st->sd) {
            // ...providing a wall-clock-based seed when not set.
            clock_gettime(CLOCK_REALTIME, &st->ts);
            st->sd = st->ts.tv_sec + st->ts.tv_nsec;
        }
    }

    // Invoke the algorithm recording the result in for tuna_post_cost()
    st->ik = st->al(nk, ks, &st->sd);

    // Glimpse at the clock so we may compute elapsed time in tuna_post()
    clock_gettime(TUNA_CLOCK, &st->ts);

    return st->ik;
}

inline
void tuna_post_cost(tuna_site*  st,
                    tuna_kernel* ks,
                    const double cost)
{
    tuna_kernel_obs(ks + st->ik, cost);
}

double tuna_post(tuna_site*  st,
                 tuna_kernel* ks)
{
    // Glimpse at the clock and compute double-valued elapsed time
    struct timespec te;
    clock_gettime(TUNA_CLOCK, &te);
    double elapsed = te.tv_nsec - st->ts.tv_nsec; // Nanoseconds...
    elapsed *= 1e-9;                              // ...to seconds
    elapsed += te.tv_sec - st->ts.tv_sec;         // ...plus seconds

    // Delegate recording the observation to tuna_post_cost()
    tuna_post_cost(st, ks, elapsed);

    return elapsed;
}
