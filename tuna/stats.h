/*
 * Copyright (C) 2012, 2013 Rhys Ulerich
 * Copyright (C) 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_STATS_H
#define TUNA_STATS_H

/** @file
 * Provides statistical accumulators for special cases of interest.
 */

#include <math.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* TODO Defend against overflowing the counter tuna_stats.n */

/**
 * Accumulates running mean and variance details from a data stream.  Fill
 * storage with zeros, e.g. from POD zero initialization, to construct or reset
 * an instance.
 *
 * Adapted from <a
 * href="https://red.ices.utexas.edu/projects/suzerain/wiki">Suzerain</a>'s
 * <code>suzerain::running_statistics</code> class which adapted it from
 * http://www.johndcook.com/standard_deviation.html.  Rewritten in C and
 * extended to permit merging statistics from multiple instances.  Storage
 * overhead reduced relative to Cook's presentation.
 */
typedef struct tuna_stats
{
    double m;
    double s;
    size_t n;
} tuna_stats;

/** Obtain the running number of samples provided thus far. */
static inline
size_t tuna_stats_cnt(const tuna_stats * const t)
{ return t->n; }

/** Obtain the running mean. */
static inline
double tuna_stats_avg(const tuna_stats * const t)
{ return t->n ? t->m : NAN; }

/** Obtain the running sample variance. */
static inline
double tuna_stats_var(const tuna_stats * const t)
{ return t->n ? (t->n > 1 ? t->s / (t->n - 1): 0) : NAN; }

/** Obtain the running sample standard deviation. */
static inline
double tuna_stats_std(const tuna_stats * const t)
{ return sqrt(tuna_stats_var(t)); }

/** Accumulate a new observation \c x into statistics \c t. */
static inline
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

/** Incorporate running information from another instance. */
tuna_stats* tuna_stats_merge(      tuna_stats * const dst,
                             const tuna_stats * const src);

/**
 * Compute a one-sided Welch t-test that \c a is greater than \c b.
 * See http://en.wikipedia.org/wiki/Welch's_t_test for background.
 */
double tuna_stats_welch1(const tuna_stats * const a,
                         const tuna_stats * const b);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_STATS_H */
