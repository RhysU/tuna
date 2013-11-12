/*
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * \mainpage
 *
 * Tuna provides lightweight autotuning over semantically indistinguishable
 * kernels.
 *
 * See the current <a
 * href="https://github.com/RhysU/tuna/blob/master/README.rst">README</a> for a
 * more detailed overview and http://github.com/RhysU/ar for project
 * information.
 */

/** \file
 * Tuna declarations.
 */

#ifndef TUNA_H
#define TUNA_H

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>

/**
 * Macros to count the number of array elements at compile time.  Provides both
 * a type-unsafe C version and a type-safe C++ version.  The C++ version is
 * from Ivan J. Johnson's article "Counting Array Elements at Compile Time"
 * published 6 March 2007 in Dr. Dobb's (http://drdobbs.com/cpp/197800525).
 * @{
 */

#ifndef __cplusplus

/** Count the number of elements in an array at compile time */
#define tuna_countof(x) (sizeof(x)/sizeof((x)[0]))

#else /* __cplusplus */

/** Type safe count the number of elements in an array at compile time */
#define tuna_countof(x)  (                                                  \
    0*sizeof(reinterpret_cast<const ::tuna::BAD_ARGUMENT_TO_COUNTOF*>(x)) + \
    0*sizeof(::tuna::BAD_ARGUMENT_TO_COUNTOF::check_type((x), &(x))     ) + \
    sizeof(x)/sizeof((x)[0])   )

#ifndef TUNA_PARSED_BY_DOXYGEN
namespace tuna {

class BAD_ARGUMENT_TO_COUNTOF
{
public:
    class Is_pointer;
    class Is_array {};
    template<typename T>
    static Is_pointer check_type(const T*, const T* const*);
    static Is_array check_type(const void*, const void*);
};

}
#endif /* TUNA_PARSED_BY_DOXYGEN */

#endif /* __cplusplus */

/** @} */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Provides statistical accumulators for special cases of interest.
 * @{
 */

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
typedef struct tuna_stats {
    double m;
    double s;
    size_t n;
} tuna_stats;

/** Obtain the running number of samples provided thus far. */
size_t
tuna_stats_cnt(const tuna_stats* const t);

/** Obtain the running mean. */
double
tuna_stats_avg(const tuna_stats* const t);

/**
 * Quickly obtain the running mean
 * on precondition <code>tuna_stats_cnt(t) > 0</code>.
 */
double
tuna_stats_fastavg(const tuna_stats* const t);

/** Obtain the running sample variance. */
double
tuna_stats_var(const tuna_stats* const t);

/**
 * Quickly obtain the running sample variance
 * on precondition <code>tuna_stats_cnt(t) > 1</code>.
 */
double
tuna_stats_fastvar(const tuna_stats* const t);

/** Obtain the running sample standard deviation. */
double
tuna_stats_std(const tuna_stats* const t);

/**
 * Quickly obtain the running sample standard deviation
 * on precondition <code>tuna_stats_cnt(t) > 1</code>.
 */
double
tuna_stats_faststd(const tuna_stats* const t);

/** Obtain the running sum. */
double
tuna_stats_sum(const tuna_stats* const t);

/**
 * Accumulate a new observation \c x into statistics \c t
 * on precondition <code>tuna_stats_cnt(t) > 0</code>.
 */
tuna_stats*
tuna_stats_fastobs(tuna_stats* const t,
                   const double x);

/** Accumulate a new observation \c x into statistics \c t. */
tuna_stats*
tuna_stats_obs(tuna_stats* const t,
               const double x);

/**
 * Accumulate \c N distinct observations <code>x[0]</code>, ...,
 * <code>x[N-1]</code> into statistics \c t.
 */
tuna_stats*
tuna_stats_nobs(tuna_stats* const t,
                const double* x,
                size_t N);

/** Incorporate running information from another instance. */
tuna_stats*
tuna_stats_merge(tuna_stats* const dst,
                 const tuna_stats* const src);

/** @} */

/**
 * Gathers elapsed time information for compute kernels.
 * @{
 */

/**
 * Accumulates runtime information about the performance of a compute kernel.
 * Fill storage with zeros, e.g. from POD zero initialization, to construct or
 * reset an instance.
 */
typedef struct tuna_kernel {
    double     outliers[3];  /**< Invariantly-sorted greatest outliers.  */
    tuna_stats stats;        /**< Accumulated statistics sans outliers. */
} tuna_kernel;

/**
 * Record a new cost observation \c t about kernel \c k.  If cost \c t is
 * identically zero, no observation is recorded.  Cost might be elapsed time,
 * but it might also be some other performance metric.  Regardless of what is
 * chosen, smaller should mean better.
 */
tuna_kernel*
tuna_kernel_obs(tuna_kernel* const k,
                double t);

/**
 * Incorporate all cost information recorded about kernel \c k into \c s,
 * including any outliers otherwise discarded from consideration.
 */
tuna_stats*
tuna_kernel_merge(tuna_stats* const s,
                  const tuna_kernel* const k);

/** @} */

/**
 * Statistical building blocks.
 * @{
 */

/**
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * \internal
 * Approximate the inverse of the standard normal per Peter John Acklam
 * (http://home.online.no/~pjacklam/notes/invnorm/) with error handling
 * following Chad Sprouse's implementation given at
 * (http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c).
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 * \endinternal
 */
double
tuna_ltqnorm(const double p);

/**
 * Computes the lower tail probability for Student's t-distribution.
 *
 * \internal
 * The algorithm is AS 3 from Applied Statistics (1968) Volume 17 Page 189 as
 * reported by <a href="http://lib.stat.cmu.edu/apstat/">StatLib</a>.  The <a
 * href="http://lib.stat.cmu.edu/apstat/3">Fortran source</a> has been
 * processed by <a href="http://www.netlib.org/f2c/">f2c</a> and then cleaned.
 * Double precision has been used throughout, reentrant operation made
 * possible, and parameters are passed by value.
 * \endinternal
 *
 * \param t  Threshold of interest
 * \param nu Strictly positive number of degrees of freedom.
 *
 * \return The lower tail probability for \c t given \c nu.
 */
double
tuna_as3(double t, int nu);

/** @} */

/**
 * Pseudo-random number generation.
 * @{
 */

/**
 * What state is required for tuna_u01() and tuna_n01() to be re-entrant safe
 * and orthogonal to all other pseudo-random generators that might be in use?
 */
typedef unsigned int tuna_seed;

/**
 * Retrieves a default seed.  If the whitespace-trimmed environment variable
 * <code>TUNA_SEED</code> can be parsed as a seed, it is used.  Otherwise, a
 * time-based seed is returned.
 */
tuna_seed
tuna_seed_default();

/** Generate a uniform draw from <tt>[0, 1]</tt>. */
double
tuna_rand_u01(tuna_seed* sd);

/** Generate a draw from <tt>N(0, 1)</tt>. */
double
tuna_rand_n01(tuna_seed* sd);

/** @} */

/**
 * Provides variants of <a
 * href="http://en.wikipedia.org/wiki/Welch's_t_test">Welch's t test</a>.
 * @{
 */

/**
 * Compute the Welch t-statistic \c t and degrees of freedom \c nu given
 * samples \c A and \c B.  When only the t-statistic is required, prefer method
 * \ref tuna_welch_t to this one.
 *
 * \param[in ] xA  Mean of \c A
 * \param[in ] sA2 Variance of \c A
 * \param[in ] nA  Number of observations of A
 * \param[in ] xB  Mean of \c B
 * \param[in ] sB2 Variance of \c B
 * \param[in ] nB  Number of observations of B
 * \param[out] t   Welch's t-statistic
 * \param[out] nu  Number of degrees of freedom per
 *                 Welch--Satterthwaite equation
 */
void
tuna_welch(double xA, double sA2, size_t nA,
           double xB, double sB2, size_t nB,
           double* const t, double* const nu);

/**
 * Compute the Welch t-statistic given samples \c A and \c B.
 *
 * \param xA  Mean of \c A
 * \param sA2 Variance of \c A
 * \param nA  Number of observations of A
 * \param xB  Mean of \c B
 * \param sB2 Variance of \c B
 * \param nB  Number of observations of B
 *
 * \return Welch's t-statistic.
 */
double
tuna_welch_t(double xA, double sA2, size_t nA,
             double xB, double sB2, size_t nB);

/**
 * Compute a one-sided Welch t-test that \c A is greater than \c B using a
 * \f$\nu\to\infty\f$ approximation.  Accordingly, the t-statistic is compared
 * against the normal distribution.  This is faster than using the
 * t-distribution, but it produces poor results for small sample sizes.
 *
 * \param xA  Mean of \c A
 * \param sA2 Variance of \c A
 * \param nA  Number of observations of A
 * \param xB  Mean of \c B
 * \param sB2 Variance of \c B
 * \param nB  Number of observations of B
 *
 * \return Approximate p-value.
 */
double
tuna_welch1_nuinf(double xA, double sA2, size_t nA,
                  double xB, double sB2, size_t nB);

/**
 * Compute a one-sided Welch t-test that \c A is greater than \c B using a
 * coarse approximation to the t-distribution.  The quality of the
 * approximation improves as \f$\nu\to\infty\f$.
 *
 * \param xA  Mean of \c A
 * \param sA2 Variance of \c A
 * \param nA  Number of observations of A
 * \param xB  Mean of \c B
 * \param sB2 Variance of \c B
 * \param nB  Number of observations of B
 *
 * \return Approximate p-value.
 */
double
tuna_welch1_approx(double xA, double sA2, size_t nA,
                   double xB, double sB2, size_t nB);

/**
 * Compute a one-sided Welch t-test that \c A is greater than \c B.
 *
 * \param xA  Mean of \c A
 * \param sA2 Variance of \c A
 * \param nA  Number of observations of A
 * \param xB  Mean of \c B
 * \param sB2 Variance of \c B
 * \param nB  Number of observations of B
 *
 * \return Exact p-value.
 */
double
tuna_welch1(double xA, double sA2, size_t nA,
            double xB, double sB2, size_t nB);

/** @} */

/**
 * Available autotuning algorithms.
 * @{
 */

/**
 * Type signature for all autotuning algorithms.  As a precondition, at least
 * two observations are available on each kernel prior to invocation.  This
 * choice permits using branchless query functions like \ref tuna_stats_var().
 *
 * \param[in   ] nk How many alternatives are under consideration?
 * \param[inout] ks Tracks information about \c nk alternatives.
 *                  Must be stored contiguously in memory.
 * \param[inout] sd Localized pseudo-random number generator state.
 *
 * \return The zero-based index of the kernel that has been selected.
 */
typedef int (*tuna_algo)(const int nk,
                         const tuna_kernel* ks,
                         tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1_nuinf. */
int
tuna_algo_welch1_nuinf(const int nk,
                       const tuna_kernel* ks,
                       tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1. */
int
tuna_algo_welch1(const int nk,
                 const tuna_kernel* ks,
                 tuna_seed* sd);

/**
 * An "autotuning" algorithm always selecting index zero.
 * Useful for testing/debugging.  See also \ref tuna_seed_default().
 */
int
tuna_algo_zero(const int nk,
               const tuna_kernel* ks,
               tuna_seed* sd);

/**
 * Retrieve a default algorithm when left unspecified.  If the
 * whitespace-trimmed environment variable <code>TUNA_ALGO</code>
 * case-insensitively names an algorithm without the <code>tuna_algo_</code>
 * prefix, that algorithm will be used.
 */
tuna_algo
tuna_algo_default(void);

/** @} */

/** The clock_gettime(2) used for internally-managed timing. */
#define TUNA_CLOCK CLOCK_PROCESS_CPUTIME_ID
#ifndef CLOCK_PROCESS_CPUTIME_ID
# error "CLOCK_PROCESS_CPUTIME_ID unavailable"
#endif

/**
 * High-level APIs for autotuning.
 * @{
 */

/**
 * Kernel-independent state required for each autotuning site.
 * Contents are internally managed but a non-opaque type
 * is used so the compiler may compute this POD type's size.
 */
typedef struct tuna_site {
    tuna_algo       al;  /**< The chosen tuning algorithm.                */
    tuna_seed       sd;  /**< Random number generator state.              */
    int             ik;  /**< Index of the most recently selected kernel. */
    struct timespec ts;  /**< Records clock_gettime(2) in tuna_pre().     */
} tuna_site;

/**
 * Invoke the currently selected autotuning algorithm.
 *
 * \param[inout] st Information local to one autotuning site.
 * \param[in   ] ks Tracks information about \c nk alternatives.
 *                  Must be stored contiguously in memory.
 * \param[in   ] nk How many alternatives are under consideration?
 *
 * \return The zero-based index of the kernel which should be selected.
 */
int
tuna_pre(tuna_site* st,
         const tuna_kernel* ks,
         const int nk);

/**
 * Record the results from the last autotuned kernel invocation
 * using internally-managed elapsed time via \ref TUNA_CLOCK.
 *
 * \param[inout] st   Information local to one autotuning site.
 * \param[in   ] ks   Tracks information about the alternatives.
 *
 * \return The inclusive time in seconds as measured by \ref TUNA_CLOCK.
 */
double
tuna_post(tuna_site*  st,
          tuna_kernel* ks);

/**
 * Record the last autotuned kernel invocation
 * using a user-provided \c cost metric.
 *
 * \param[inout] st   Information local to one autotuning site.
 * \param[in   ] ks   Tracks information about the alternatives.
 * \param[in   ] cost User-provided measure of the employed kernel's cost.
 *                    This may be elapsed time or some sophisticated measure.
 *                    It should be strictly positive.  Lower means better.
 */
void
tuna_post_cost(tuna_site*  st,
               tuna_kernel* ks,
               const double cost);

/** @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
