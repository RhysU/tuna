//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc stats.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "stats.h"

// C99 extern declarations for inlined functions from stats.h
extern size_t tuna_stats_cnt(const tuna_stats * const t);
extern double tuna_stats_avg(const tuna_stats * const t);
extern double tuna_stats_var(const tuna_stats * const t);
extern double tuna_stats_std(const tuna_stats * const t);

tuna_stats* tuna_stats_obs(tuna_stats * const t, const double x)
{
    // Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.
    // Knuth shows better behavior than Welford 1962 on test data.
    const size_t n = ++(t->n);
    if (n > 1) {  // Second and subsequent invocation
        double d  = x - t->m;
        t->m     += d / n;
        t->s     += d * (x - t->m);
    } else {      // First invocation requires special treatment
        t->m = x;
        t->s = 0;
    }
    return t;
}

tuna_stats* tuna_stats_merge(      tuna_stats * const dst,
                             const tuna_stats * const src)
{
    if        (src->n == 0) {  // src contains no data
        // NOP
    } else if (dst->n == 0) {  // dst contains no data
        *dst = *src;
    } else {                   // merge src into dst
        size_t total = dst->n + src->n;
        double dM    = dst->m - src->m;  // Cancellation issues?
        dst->m       = (dst->n * dst->m + src->n * src->m) / total;
        dst->s       = (dst->n == 1 ? 0 : dst->s)
                     + (src->n == 1 ? 0 : src->s)
                     + ((dM * dM) * (dst->n * src->n)) / total;
        dst->n       = total;
    }
    return dst;
}

// Beware this implementation is atrociously bad for small sample sizes!  We
// incorrectly always assume NDOF nu \to \infty and use N(0,1) distribution.
double tuna_stats_welch1(const tuna_stats * const a,
                         const tuna_stats * const b)
{
    const double t_num  = tuna_stats_avg(a) - tuna_stats_avg(b);
    const double t_den2 = tuna_stats_var(a) / tuna_stats_cnt(a)
                        + tuna_stats_var(b) / tuna_stats_cnt(b);
    const double t      = t_num / sqrt(t_den2);
    const double p      = erfc(-t * M_SQRT1_2) / 2;
    return p;
}
