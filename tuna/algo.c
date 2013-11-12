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
        } else if (!strncasecmp(d, "zero", sizeof("zero"))) {
            return &tuna_algo_zero;
        }
    }
    return &tuna_algo_welch1_nuinf; // Default
}

int tuna_algo_welch1_nuinf(const int nk,
                           const tuna_kernel* ks,
                           tuna_seed* seed)
{
    assert(nk > 0);
    int i = 0;
    double iavg = tuna_stats_fastavg(&ks[0].stats);
    double ivar = tuna_stats_fastvar(&ks[0].stats);
    double icnt = tuna_stats_cnt    (&ks[0].stats);
    for (int j = 1; j < nk; ++j) {
        const double javg = tuna_stats_fastavg(&ks[j].stats);
        const double jvar = tuna_stats_fastvar(&ks[j].stats);
        const double jcnt = tuna_stats_cnt    (&ks[j].stats);
        const double p    = tuna_welch1_nuinf (iavg, ivar, icnt,
                                               javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int tuna_algo_welch1(const int nk,
                     const tuna_kernel* ks,
                     tuna_seed* seed)
{
    assert(nk > 0);
    int i = 0;
    double iavg = tuna_stats_fastavg(&ks[0].stats);
    double ivar = tuna_stats_fastvar(&ks[0].stats);
    double icnt = tuna_stats_cnt    (&ks[0].stats);
    for (int j = 1; j < nk; ++j) {
        const double javg = tuna_stats_fastavg(&ks[j].stats);
        const double jvar = tuna_stats_fastvar(&ks[j].stats);
        const double jcnt = tuna_stats_cnt    (&ks[j].stats);
        const double p    = tuna_welch1       (iavg, ivar, icnt,
                                               javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int tuna_algo_zero(const int nk,
                   const tuna_kernel* ks,
                   tuna_seed* seed)
{
    (void) nk;
    (void) ks;
    (void) seed;
    return 0;
}
