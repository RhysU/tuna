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
 * (A) Chunk(-based) Light(weight) (Auto)Tuna.
 *
 * See the current <a
 * href="https://github.com/RhysU/tuna/blob/master/README.rst">README</a> for a
 * more detailed overview and http://github.com/RhysU/tuna for project
 * information.
 */

/** \file
 * Tuna public API.
 * Designed to have the minimial possible set of dependencies.
 */

#ifndef TUNA_H
#define TUNA_H

#include <stddef.h>

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
namespace tuna
{

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

/**
 * Compile-time static assertion macro from Ralf Holly's <a
 * href="http://drdobbs.com/184401873">Compile Time Assertions</a>.  Most
 * compilers report "expected constant expression" or similar on failed
 * assertion.
 */
#define tuna_assert_static(e)             \
    do {                                  \
        enum { assert_static__ = 1/(e) }; \
    } while (0)

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Spinlocking utilities for protecting critical data structures.
 *
 * Spinlocks chosen as critical regions are short and lack OS-intensive calls.
 * Also, our contention should be light.  These definitions may need to be
 * modified depending on the compiler and/or platform used.
 * @{
 */

/** Provides storage necessary to support one spinlock. */
typedef volatile int tuna_spinlock;

/** Lock a \ref tuna_spinlock using a GCC-defined atomic operation. */
#define tuna_lock(spinlock)   while (__sync_lock_test_and_set(&spinlock, 1));

/** Unlock a \ref tuna_spinlock using a GCC-defined atomic operation. */
#define tuna_unlock(spinlock) while __sync_lock_release(&spinlock);

/** @} */

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

/**
 * Quickly obtain the running number of samples provided thus far.
 */
size_t
tuna_stats_fastcnt(const tuna_stats* const t);

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
 * on precondition <code>tuna_stats_fastcnt(t) > 1</code>.
 */
double
tuna_stats_faststd(const tuna_stats* const t);

/** Obtain the running sum. */
double
tuna_stats_sum(const tuna_stats* const t);

/**
 * Quickly obtain the running sum
 * on precondition <code>tuna_stats_fastcnt(t) > 0</code>.
 */
double
tuna_stats_fastsum(const tuna_stats* const t);

/**
 * Quickly accumulate a new observation \c x into statistics \c t on
 * precondition <code>tuna_stats_fastcnt(t) > 0</code>.
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
 * Gather and manipulate cost information for compute chunks.
 * @{
 */

/**
 * Accumulates runtime information about the performance of a compute chunk.
 * Fill storage with zeros, e.g. from POD zero initialization, to construct or
 * reset an instance.
 */
typedef struct tuna_chunk {
    double     outliers[3];  /**< Invariantly-sorted greatest outliers.  */
    tuna_stats stats;        /**< Accumulated statistics sans outliers. */
} tuna_chunk;

/**
 * Record a new cost observation \c t about chunk \c k.  If cost \c t is
 * identically zero, no observation is recorded.  Cost might be elapsed time,
 * but it might also be some other performance metric.  Regardless of what is
 * chosen, smaller should mean better.
 */
tuna_chunk*
tuna_chunk_obs(tuna_chunk* const k,
               double t);

/**
 * Incorporate all cost information recorded about chunk \c k into \c s,
 * including any outliers otherwise discarded from consideration.
 */
tuna_stats*
tuna_chunk_merge(tuna_stats* const s,
                 const tuna_chunk* const k);

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
 * two observations are available on each chunk prior to invocation.  This
 * choice permits using branchless query functions like \ref tuna_stats_var().
 *
 * \param[in   ] nk How many alternatives are under consideration?
 * \param[inout] ks Tracks information about \c nk alternatives.
 *                  Must be stored contiguously in memory.
 * \param[inout] sd Localized pseudo-random number generator state.
 *
 * \return The zero-based index of the chunk that has been selected.
 */
typedef int (*tuna_algo)(const int nk,
                         const tuna_chunk ks[],
                         tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1_nuinf. */
int
tuna_algo_welch1_nuinf(const int nk,
                       const tuna_chunk ks[],
                       tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1. */
int
tuna_algo_welch1(const int nk,
                 const tuna_chunk ks[],
                 tuna_seed* sd);

/**
 * An "autotuning" algorithm always selecting index zero.
 * Useful for testing/debugging.  See also \ref tuna_seed_default().
 */
int
tuna_algo_zero(const int nk,
               const tuna_chunk ks[],
               tuna_seed* sd);

/**
 * Retrieve the name of the algorithm, if known.  Otherwise, return "unknown".
 */
const char *
tuna_algo_name(tuna_algo al);

/**
 * Retrieve a default algorithm when one is left unspecified.  If <code>nk <
 * 2</code>, \ref tuna_algo_zero() is returned.  If the whitespace-trimmed
 * environment variable <code>TUNA_ALGO</code> case-insensitively names an
 * algorithm without the <code>tuna_algo_</code> prefix, that algorithm will be
 * used.  Otherwise, a sensible default which may or may not take into account
 * \c nk is chosen.
 */
tuna_algo
tuna_algo_default(const int nk);

/** @} */

/**
 * High-level APIs for autotuning.
 * @{
 */

/**
 * Chunk-independent state maintained \e once for each autotuning site.
 * Members are managed by Tuna but a non-opaque type is used so the compiler
 * may compute this POD type's size to permit \c static instances.
 */
typedef struct tuna_site {
    tuna_algo al; /**< The chosen tuning algorithm.   */
    tuna_seed sd; /**< Random number generator state. */
} tuna_site;

/**
 * Chunk-independent state maintained for \e every autotuned invocation.
 * Members are managed by Tuna but a non-opaque type is used so the compiler
 * may compute this POD type's size to permit stack-allocated instances.  Stack
 * allocation is recommended to allow reentrant usage of the library.
 */
typedef struct tuna_stack {
    int  ik;      /**< Index of the most recently selected chunk. */
    char ts[16];  /**< Stores clock_gettime(2) within tuna_pre(). */
} tuna_stack;

/**
 * Invoke the currently selected autotuning algorithm.  Only \ref
 * tuna_post_cost() may be invoked after the selected chunk completes
 * executing.
 *
 * \param[inout] si Durable information local to autotuning site.
 * \param[inout] st Stack-based information stored from this invocation.
 * \param[in   ] ks Durable information tracking for \c nk alternatives.
 * \param[in   ] nk How many alternatives are under consideration?
 *
 * \return The zero-based index of the chunk which should be selected.
 */
int
tuna_pre_cost(tuna_site* si,
              tuna_stack* st,
              const tuna_chunk ks[],
              const int nk);

/**
 * Record the last autotuned chunk invocation using a user-provided \c cost
 * metric.  Either \ref tuna_pre() or \ref tuna_pre_cost() should have been
 * invoked beforehand.
 *
 * \param[in   ] st   Stack-based information from one autotuned invocation.
 * \param[inout] ks   Updated with information learned from chosen alternative.
 * \param[in   ] cost User-provided measure of the employed chunk's cost.
 *                    This may be elapsed time or some sophisticated measure.
 *                    It should be strictly positive.  Lower means better.
 */
void
tuna_post_cost(const tuna_stack* st,
               tuna_chunk ks[],
               const double cost);

/**
 * Invoke the currently selected autotuning algorithm.  Either \ref tuna_post()
 * or \ref tuna_post_cost() may be invoked after the selected chunk completes
 * executing.
 *
 * \param[inout] si Durable information local to one autotuning site.
 * \param[inout] st Stack-based information stored from this invocation.
 * \param[in   ] ks Durable information tracking for \c nk alternatives.
 * \param[in   ] nk How many alternatives are under consideration?
 *
 * \return The zero-based index of the chunk which should be selected.
 */
int
tuna_pre(tuna_site* si,
         tuna_stack* st,
         const tuna_chunk ks[],
         const int nk);

/**
 * Record the results from the last autotuned chunk invocation using
 * internally-managed elapsed time via an appropriate clock.  Method \ref
 * tuna_pre() should have been invoked before the chunk began executing.
 *
 * \param[in   ] st Stack-based information from one autotuned invocation.
 * \param[inout] ks Updated with information learned from chosen alternative.
 *
 * \return The inclusive process time in seconds.
 */
double
tuna_post(const tuna_stack* st,
          tuna_chunk ks[]);

/** @} */

/**
 * Human-readable output for summarizing runtime behavior.
 *
 * The simpler, non-<code>varargs</code> <code>*_fprint</code> routines are
 * intended for non-C purposes while the more flexible but less typesafe
 * <code>*_fprintf</code> variants are intended for C usage.  The \c stream
 * parameters to these methods are \e not flushed.
 * @{
 */

/**
 * Output a single status line about \ref tuna_chunk k to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * @param stream <code>FILE*</code> on which output is produced.
 * @param k      Chunk for which information is output.
 * @param prefix A string used to prefix the output.
 *
 * @return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_chunk_fprint(void *stream,
                  const tuna_chunk* k,
                  const char *prefix);

/**
 * Output a single status line about \ref tuna_site st to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * @param stream <code>FILE*</code> on which output is produced.
 * @param si     Site for which information is output.
 * @param prefix A string used to prefix the output.
 *
 * @return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_site_fprint(void *stream,
                 const tuna_site* si,
                 const char *prefix);

/**
 * Output <code>nk+1</code> status line(s) about \ref tuna_site st and
 * associated \ref tuna_chunk <code>ks[0]</code>, ..., <code>ks[nk-1]</code> to
 * <code>FILE*</code> stream prefixed by \c prefix.
 *
 * @param stream <code>FILE*</code> on which output is produced.
 * @param si     Site for which information is output.
 * @param ks     Chunks about which information is output.
 * @param nk     Number of contiguous chunks in \c ks.
 * @param prefix A string used to prefix the output.
 * @param labels If non-NULL and non-trivial, <code>labels[ik]</code> labels
 *               the <code>ik</code>th chunk in the output.  Otherwise,
 *               each chunk is labelled with its index <code>ik</code>.
 *
 * @return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_fprint(void *stream,
            const tuna_site* si,
            const tuna_chunk ks[],
            const int nk,
            const char *prefix,
            const char *labels[]);

/**
 * \copybrief tuna_chunk_fprint
 *
 * @param stream <code>FILE*</code> on which output is produced.
 * @param k      Chunk for which information is output.
 * @param format A <code>printf</code>-style specifying prefixing
 *               the output.  Any subsequent arguments are consumed
 *               as if <code>fprintf(stream, format, ...)</code> had
 *               been called.
 *
 * @return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_chunk_fprintf(void *stream,
                   const tuna_chunk* k,
                   const char *format,
                   ...);

/**
 * \copybrief tuna_site_fprint
 *
 * @param stream <code>FILE*</code> on which output is produced.
 * @param si     Site for which information is output.
 * @param format A <code>printf</code>-style specifying prefixing
 *               the output.  Any subsequent arguments are consumed
 *               as if <code>fprintf(stream, format, ...)</code> had
 *               been called.
 *
 * @return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_site_fprintf(void *stream,
                  const tuna_site* si,
                  const char *format,
                  ...);

/** @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
