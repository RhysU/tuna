/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_WELCH_H
#define TUNA_WELCH_H

/** @file
 * Provides variants of <a
 * href="http://en.wikipedia.org/wiki/Welch's_t_test">Welch's t test</a>.
 */

#include <math.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the Welch t-statistic comparing samples \c A and \c B.
 *
 * @param xA  Mean of \c A
 * @param sA2 Variance of \c A
 * @param nA  Number of observations of A
 * @param xB  Mean of \c B
 * @param sB2 Variance of \c B
 * @param nB  Number of observations of B
 *
 * @return Welch's t-statistic.
 */
static inline
double tuna_welch_t(double xA, double sA2, size_t nA,
                    double xB, double sB2, size_t nB)
{
    double t_num  = xA - xB;
    double t_den2 = sA2 / nA + sB2 / nB;
    double t      = t_num / sqrt(t_den2);
    return t;
}

/**
 * Compute a one-sided Welch t-test that \c A is greater than \c B using a
 * \f$\nu\to\infty\f$ approximation.  Accordingly, the t-statistic is compared
 * against the normal distribution.  This is faster than using the
 * t-distribution, but it produces poor results for small sample sizes.
 *
 * @param xA  Mean of \c A
 * @param sA2 Variance of \c A
 * @param nA  Number of observations of A
 * @param xB  Mean of \c B
 * @param sB2 Variance of \c B
 * @param nB  Number of observations of B
 *
 * @return Approximate p-value.
 */
static inline
double tuna_welch1_nuinf(double xA, double sA2, size_t nA,
                         double xB, double sB2, size_t nB)
{
    return erfc(-tuna_welch_t(xA, sA2, nA, xB, sB2, nB) * M_SQRT1_2) / 2;
}

//// TODO Implement a one-sided test using t-distn facts to broaden variances
// static inline
// double tuna_welch1_approx(double xA, double sA2, size_t nA,
//                           double xB, double sB2, size_t nB)

//// TODO Implement a one-sided test using proper degrees of freedom
// static inline
// double tuna_welch1(double xA, double sA2, size_t nA,
//                    double xB, double sB2, size_t nB)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_WELCH_H */
