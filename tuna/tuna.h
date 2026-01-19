/*
 * Copyright (C) 2011, 2012, 2013, 2026 Rhys Ulerich
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
 * The implementation is not thread-safe but it can be used in thread-safe
 * fashion provided that some external lock guards access to struct instances.
 */

#ifndef TUNA_H
#define TUNA_H

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>

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
 * Provides statistical accumulators for special cases of interest.
 * @{
 */

/* TODO Defend against overflowing the counter tuna_stats.n */

/**
 * Accumulates running mean and variance details from a data stream.
 * Fill storage with zeros, e.g. from POD zero initialization,
 * to construct or reset an instance.
 *
 * Adapted from <a
 * href="https://red.ices.utexas.edu/projects/suzerain/wiki">Suzerain</a>'s
 * <code>suzerain::running_statistics</code> class which adapted it from
 * http://www.johndcook.com/standard_deviation.html.  Rewritten in C and
 * extended to permit merging statistics from multiple instances.  Storage
 * overhead reduced relative to Cook's presentation.
 */
typedef struct tuna_stats {
    double        m;
    double        s;
    size_t        n;
} tuna_stats;

/** Obtain the running number of samples provided thus far. */
size_t
tuna_stats_cnt(const tuna_stats* stats);

/** Obtain the running mean. */
double
tuna_stats_avg(const tuna_stats* stats);

/** Obtain the running sample variance. */
double
tuna_stats_var(const tuna_stats* stats);

/** Obtain the running sample standard deviation. */
double
tuna_stats_std(const tuna_stats* stats);

/** Obtain the running sum. */
double
tuna_stats_sum(const tuna_stats* stats);

/**
 * Contains the result of tuna_stats_mom().
 */
typedef struct tuna_stats_mom_result {
    size_t n;    /**< Number of samples */
    double avg;  /**< Mean */
    double var;  /**< Sample variance */
} tuna_stats_mom_result;

/**
 * Obtain all running moments at once.
 * \param[in] stats Instance of interest
 * \return A struct containing the number of samples, mean, and variance.
 */
tuna_stats_mom_result
tuna_stats_mom(const tuna_stats* stats);

/** Accumulate a new observation \c x into statistics \c stats. */
void
tuna_stats_obs(tuna_stats* stats,
               const double x);

/**
 * Accumulate \c N distinct observations <code>x[0]</code>, ...,
 * <code>x[N-1]</code> into statistics \c stats.
 */
void
tuna_stats_nobs(tuna_stats* stats,
                const double* x,
                size_t N);

/** Incorporate running information from another instance. */
void
tuna_stats_merge(tuna_stats* dst,
                 const tuna_stats* src);

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
 * Record a new cost observation \c cost about chunk \c chunk.  If cost \c cost is
 * identically zero, no observation is recorded.  Cost might be elapsed time,
 * but it might also be some other performance metric.  Regardless of what is
 * chosen, smaller should mean better.
 */
void
tuna_chunk_obs(tuna_chunk* chunk,
               double cost);

/**
 * Incorporate all cost information recorded about chunk \c chunk into \c stats,
 * including any outliers otherwise discarded from consideration.
 */
void
tuna_chunk_merge(tuna_stats* stats,
                 const tuna_chunk* chunk);

/** @} */

/**
 * Pseudo-random number generation.
 * @{
 */

/**
 * The state required so that tuna_u01() and tuna_n01() can be re-entrant safe
 * and orthogonal to all other pseudo-random generators that might be in use.
 */
typedef unsigned int tuna_state;

/**
 * Retrieves a default state value.  If the whitespace-trimmed environment
 * variable <code>TUNA_SEED</code> can be parsed as a seed, that seed is used
 * to initialize the state.  Otherwise, a time-based seed is used.
 */
tuna_state
tuna_state_default(void);

/** Generate a uniform draw from <tt>[0, 1]</tt>. */
double
tuna_rand_u01(tuna_state* state);

/** Generate a draw from <tt>N(0, 1)</tt>. */
double
tuna_rand_n01(tuna_state* state);

/** @} */

/**
 * Provides variants of <a
 * href="http://en.wikipedia.org/wiki/Welch's_t_test">Welch's t test</a>.
 * @{
 */

/**
 * Contains the result of tuna_welch().
 */
typedef struct tuna_welch_result {
    double t;   /**< Welch's t-statistic */
    double nu;  /**< Number of degrees of freedom per Welch--Satterthwaite equation */
} tuna_welch_result;

/**
 * Compute the Welch t-statistic \c t and degrees of freedom \c nu given
 * samples \c A and \c B.  When only the t-statistic is required, prefer method
 * \ref tuna_welch_t to this one.
 *
 * \param xA  Mean of \c A
 * \param sA2 Variance of \c A
 * \param nA  Number of observations of A
 * \param xB  Mean of \c B
 * \param sB2 Variance of \c B
 * \param nB  Number of observations of B
 *
 * \return A struct containing the t-statistic and degrees of freedom.
 */
tuna_welch_result
tuna_welch(double xA, double sA2, size_t nA,
           double xB, double sB2, size_t nB);

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
 * Type signature for all \e deterministic autotuning algorithms.
 * That is, all randomness occurs outside such algorithmic calls.
 *
 * \param[in] nchunk How many alternatives are under consideration?
 * \param[in] chunks Tracks information about \c nchunk alternatives.
 *                   Must be stored contiguously in memory.
 * \param[in] u01    Provides \c nchunk uniform random draws on [0,1]
 *                   for consumption by the algorithm.
 *
 * \return The zero-based index of the chunk that has been selected.
 */
typedef size_t (*tuna_algo_fn)(const size_t nchunk,
                               const tuna_chunk *chunks,
                               const double *u01);

/**
 * Autotuning algorithm implementation carrying its own name.
 */
typedef struct tuna_algo_impl {
    const char*   name; /**< The algorithm's name */
    tuna_algo_fn  fn;   /**< The algorithm function */
} tuna_algo_impl;

/**
 * Handle to an autotuning algorithm.
 */
typedef const tuna_algo_impl* tuna_algo;

/**
 * Retrieve the name of the algorithm.
 */
const char*
tuna_algo_name(tuna_algo algo);

/**
 * Retrieve a default algorithm when one is left unspecified.  If <code>nchunk <
 * 2</code>, the "zero" algorithm is returned.  If the whitespace-trimmed
 * environment variable <code>TUNA_ALGO</code> case-insensitively names an
 * algorithm (e.g., "welch1", "welch1_nuinf", "zero"), that algorithm will be
 * used.  Otherwise, a sensible default which may or may not take into account
 * \c nchunk is chosen.
 */
tuna_algo
tuna_algo_default(const size_t nchunk);

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
    tuna_algo     al; /**< The chosen tuning algorithm.           */
    tuna_state    st; /**< Random number generator state.         */
} tuna_site;

/**
 * Chunk-independent state maintained for \e every autotuned invocation.
 * Members are managed by Tuna but a non-opaque type is used so the compiler
 * may compute this POD type's size to permit stack-allocated instances.  Stack
 * allocation is recommended to allow reentrant usage of the library.
 */
typedef struct tuna_stack {
    size_t  ik;   /**< Index of the most recently selected chunk. */
    char ts[16];  /**< Stores clock_gettime(2) within tuna_pre(). */
} tuna_stack;

/**
 * Invoke the currently selected autotuning algorithm.  Only \ref
 * tuna_post_cost() may be invoked after the selected chunk completes
 * executing.
 *
 * \param[inout] site   Durable information local to autotuning site.
 * \param[inout] stack  Stack-based information stored from this invocation.
 * \param[in   ] chunks Durable information tracking for \c nchunk alternatives.
 * \param[in   ] nchunk How many alternatives are under consideration?
 *
 * \return The zero-based index of the chunk which should be selected.
 */
size_t
tuna_pre_cost(tuna_site* site,
              tuna_stack* stack,
              const tuna_chunk *chunks,
              const size_t nchunk);

/**
 * Record the last autotuned chunk invocation using a user-provided \c cost
 * metric.  Either \ref tuna_pre() or \ref tuna_pre_cost() should have been
 * invoked beforehand.
 *
 * \param[in   ] stack  Stack-based information from one autotuned invocation.
 * \param[inout] chunks Updated with information learned from chosen alternative.
 * \param[in   ] cost   User-provided measure of the employed chunk's cost.
 *                      This may be elapsed time or some sophisticated measure.
 *                      It should be strictly positive.  Lower means better.
 */
void
tuna_post_cost(const tuna_stack* stack,
               tuna_chunk *chunks,
               const double cost);

/**
 * Invoke the currently selected autotuning algorithm.  Either \ref tuna_post()
 * or \ref tuna_post_cost() may be invoked after the selected chunk completes
 * executing.
 *
 * \param[inout] site   Durable information local to one autotuning site.
 * \param[inout] stack  Stack-based information stored from this invocation.
 * \param[in   ] chunks Durable information tracking for \c nchunk alternatives.
 * \param[in   ] nchunk How many alternatives are under consideration?
 *
 * \return The zero-based index of the chunk which should be selected.
 */
size_t
tuna_pre(tuna_site* site,
         tuna_stack* stack,
         const tuna_chunk *chunks,
         const size_t nchunk);

/**
 * Record the results from the last autotuned chunk invocation using
 * internally-managed elapsed time via an appropriate clock.  Method \ref
 * tuna_pre() should have been invoked before the chunk began executing.
 *
 * \param[in   ] stack  Stack-based information from one autotuned invocation.
 * \param[inout] chunks Updated with information learned from chosen alternative.
 *
 * \return The inclusive process time in seconds.
 */
double
tuna_post(const tuna_stack* stack,
          tuna_chunk *chunks);

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
 * Output <code>nchunk+1</code> status line(s) about \ref tuna_site site and
 * associated \ref tuna_chunk <code>chunks[0]</code>, ..., <code>chunks[nchunk-1]</code> to
 * <code>FILE*</code> stream prefixed by \c prefix.
 *
 * \param stream <code>FILE*</code> on which output is produced.
 * \param site   Site for which information is output.
 * \param chunks Chunks about which information is output.
 * \param nchunk Number of contiguous chunks in \c chunks.
 * \param prefix A string used to prefix the output.
 * \param labels If non-NULL and non-trivial, <code>labels[ik]</code> labels
 *               the <code>ik</code>th chunk in the output.  Otherwise,
 *               each chunk is labelled with its index <code>ik</code>.
 *
 * \return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_fprint(FILE* stream,
            const tuna_site* site,
            const tuna_chunk *chunks,
            const size_t nchunk,
            const char* prefix,
            const char **labels);

/**
 * Output a single status line about \ref tuna_chunk chunk to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * \param stream <code>FILE*</code> on which output is produced.
 * \param chunk  Chunk for which information is output.
 * \param format A <code>printf</code>-style specifying prefixing
 *               the output.
 * \param ap     Variable argument list as if
 *               <code>vfprintf(stream, format, ap)</code> had
 *               been called.
 *
 * \return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_chunk_vfprintf(FILE* stream,
                    const tuna_chunk* chunk,
                    const char* format,
                    va_list ap);

/**
 * Output a single status line about \ref tuna_chunk chunk to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * \param stream <code>FILE*</code> on which output is produced.
 * \param chunk  Chunk for which information is output.
 * \param format A <code>printf</code>-style specifying prefixing
 *               the output.  Any subsequent arguments are consumed
 *               as if <code>fprintf(stream, format, ...)</code> had
 *               been called.
 *
 * \return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_chunk_fprintf(FILE* stream,
                   const tuna_chunk* chunk,
                   const char* format,
                   ...);

/**
 * Output a single status line about \ref tuna_site site to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * \param stream <code>FILE*</code> on which output is produced.
 * \param site   Site for which information is output.
 * \param format A <code>printf</code>-style specifying prefixing
 *               the output.
 * \param ap     Variable argument list as if
 *               <code>vfprintf(stream, format, ap)</code> had
 *               been called.
 *
 * \return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_site_vfprintf(FILE* stream,
                   const tuna_site* site,
                   const char* format,
                   va_list ap);

/**
 * Output a single status line about \ref tuna_site site to <code>FILE*</code>
 * stream prefixed by \c prefix.
 *
 * \param stream <code>FILE*</code> on which output is produced.
 * \param site   Site for which information is output.
 * \param format A <code>printf</code>-style specifying prefixing
 *               the output.  Any subsequent arguments are consumed
 *               as if <code>fprintf(stream, format, ...)</code> had
 *               been called.
 *
 * \return The number of characters output on success.
 *         On error, a negative value is returned.
 */
int
tuna_site_fprintf(FILE* stream,
                  const tuna_site* site,
                  const char* format,
                  ...);

/** @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
