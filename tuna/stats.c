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
extern size_t tuna_stats_cnt(const tuna_stats * const s);
extern double tuna_stats_avg(const tuna_stats * const s);
extern double tuna_stats_var(const tuna_stats * const s);
extern double tuna_stats_std(const tuna_stats * const s);

tuna_stats* tuna_stats_obs(tuna_stats * const s, const double x)
{
    // Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.
    // Knuth shows better behavior than Welford 1962 on test data.
    const size_t n = ++(s->n);
    if (n > 1) {  // Second and subsequent invocation
        double d  = x - s->M;
        s->M     += d / n;
        s->S     += d * (x - s->M);
    } else {      // First invocation requires special treatment
        s->M = x;
        s->S = 0;
    }
    return s;
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
        double dM    = dst->M - src->M;  // Cancellation issues?
        dst->M       = (dst->n * dst->M + src->n * src->M) / total;
        dst->S       = (dst->n == 1 ? 0 :   dst->S)
                     + (src->n == 1 ? 0 : src->S)
                     + ((dM * dM) * (dst->n * src->n)) / total;
        dst->n       = total;
    }
    return dst;
}
