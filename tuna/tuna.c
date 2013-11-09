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

int tuna_pre(tuna_state* st,
             const tuna_kernel* ks,
             const int nk)
{
    // Ensure a zero-initialize st argument produces good behavior by...
    if (!st->al) {
        // ...providing a default algorithm when not set, and
        st->al = &tuna_algo_welch1_nuinf;

        if (!st->sd) {
            // ...providing a wall-clock-based seed when not set.
            clock_gettime(CLOCK_REALTIME, &st->ts);
            st->sd = st->ts.tv_sec + st->ts.tv_nsec;
        }
    }

    // Invoke the algorithm recording the result in st for tuna_post()
    st->ik = st->al(nk, ks, &st->sd);

    // Glimpse at the clock so we may compute elapsed time in tuna_post()
#ifdef CLOCK_PROCESS_CPUTIME_ID
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &st->ts);
#else
    clock_gettime(CLOCK_REALTIME, &st->ts);
#endif

    return st->ik;
}
