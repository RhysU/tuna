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

#include <assert.h>
#include <time.h>

int tuna_pre(tuna_state* st,
             const tuna_kernel* ks,
             const int nk)
{
    // Provide a time-based seed if trivial
    if (!st->sd) {
        struct timespec tp;
        clock_gettime(CLOCK_REALTIME, &tp);
        st->sd = tp.tv_sec + tp.tv_nsec;
    }

    // Provide a default algorithm if trivial
    if (!st->al) {
        st->al = &tuna_algo_welch1_nuinf;
    }

    // Invoke the algorithm recording the result
    return st->ik = st->al(nk, ks, &st->sd);
}
