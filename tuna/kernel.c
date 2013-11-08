/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * \copydoc kernel.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "countof.h"
#include "kernel.h"

/**
 * Enforces <code>*a < *b<code> as a postcondition.
 * If doing so required swapping *a and *b, return 1.  Otherwise 0.
 */
static inline
int enforce_lt(double* const a, double* const b)
{
    if (*a < *b) {
        return 0;
    } else {
        double t = *a;
        *a = *b;
        *b =  t;
        return 1;
    }
}

tuna_kernel* tuna_kernel_obs(tuna_kernel* const k, double t)
{
    // First, find smallest observation among set {t, k->outliers[0], ... }
    // placing it into storage t while maintaining sorted-ness of k->outliers.
    // The loop is one bubble sort pass with possibility of short-circuiting.
    if (enforce_lt(&t, k->outliers)) {
        for (size_t i = 1;
             i < sizeof(k->outliers) / sizeof(k->outliers[0])
             && enforce_lt(k->outliers - 1 + i, k->outliers + i);
             ++i)
            ;
    }

    // Second, when non-zero, record statistics about the best observation.
    if (t) {
        tuna_stats_obs(&k->stats, t);
    }

    // Together, these two steps cause a zero-initialized tuna_kernel to
    // gather tuna_noutliers pieces of information before beginning to
    // track any statistics.  This effectively provides some "start up"
    // or "burn in" period in addition to preventing highly improbable
    // observations from unduly inflating the discovered variability.

    return k;
}

tuna_stats tuna_kernel_stats(tuna_kernel* const k)
{
    tuna_stats s = k->stats;
    return s;
}

tuna_stats* tuna_kernel_merge(tuna_stats*   const s,
                              const tuna_kernel* const k)
{
    tuna_stats_merge(s, &k->stats);
    tuna_stats_nobs(s, k->outliers, tuna_countof(k->outliers));
    return s;
}
