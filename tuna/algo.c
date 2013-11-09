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

#include "algo.h"

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <strings.h>

#include "welch.h"

// http://agentzlerich.blogspot.com/2011/01/c-header-only-unit-testing-with-fctx.html
static inline
void trim(char * const a)
{
    char *p = a, *q = a;
    while (isspace(*q))            ++q;
    while (*q)                     *p++ = *q++;
    *p = '\0';
    while (p > a && isspace(*--p)) *p = '\0';
}

// TODO Is welch1_nuinf so much better for runtime that it should be default?
//      The behavioral difference can be seen in, e.g., ./examples/smallsort.
tuna_algo tuna_algo_default(void)
{
    char * d = getenv("TUNA_ALGO");
    if (d) {
        trim(d);
        if (!strncasecmp(d, "welch1", sizeof("welch1"))) {
            return &tuna_algo_welch1;
        } else if (!strncasecmp(d, "welch1_nuinf", sizeof("welch1_nuinf"))) {
            return &tuna_algo_welch1_nuinf;
        }
    }
    return &tuna_algo_welch1_nuinf; // Default
}

int tuna_algo_welch1_nuinf(const int nk,
                           const tuna_kernel* ks,
                           tuna_seed* seed)
{
    assert(nk == 2);                            // FIXME Generalize
    const tuna_stats* const a = &ks[0].stats;   // Brevity
    const tuna_stats* const b = &ks[1].stats;
    if (tuna_stats_cnt(a) < 1) {
        return 0;
    } else if (tuna_stats_cnt(b) < 1) {
        return 1;
    } else {
        double p = tuna_welch1_nuinf(tuna_stats_avg(a),
                                     tuna_stats_var(a),
                                     tuna_stats_cnt(a),
                                     tuna_stats_avg(b),
                                     tuna_stats_var(b),
                                     tuna_stats_cnt(b));
        return p < tuna_rand_u01(seed);
    }
}

int tuna_algo_welch1(const int nk,
                     const tuna_kernel* ks,
                     tuna_seed* seed)
{
    assert(nk == 2);                            // FIXME Generalize
    const tuna_stats* const a = &ks[0].stats;   // Brevity
    const tuna_stats* const b = &ks[1].stats;
    if (tuna_stats_cnt(a) < 2) {
        return 0;
    } else if (tuna_stats_cnt(b) < 2) {
        return 1;
    } else {
        double p = tuna_welch1(tuna_stats_avg(a),
                               tuna_stats_var(a),
                               tuna_stats_cnt(a),
                               tuna_stats_avg(b),
                               tuna_stats_var(b),
                               tuna_stats_cnt(b));
        return p < tuna_rand_u01(seed);
    }
}
