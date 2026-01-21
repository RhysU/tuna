/*
 * Copyright (C) 2011, 2012, 2013, 2026 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * Tuna library implementation as a single ANSI C89 compilation unit.
 *
 * As library avoids features beyond C89 and the necessary evil of POSIX
 * 2001, having a single compilation unit may be important to permit
 * effective interprocedural optimization.  It also simplifies adoption by
 * other projects as it may simply be copied into an existing build tree.
 */

#if (_POSIX_C_SOURCE < 200112L)
#define _POSIX_C_SOURCE (200112L)
#endif

#include <tuna.h>

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#ifndef NAN
/** C89 mechanism to obtain not-a-number. */
#define NAN (sqrt(-1))
#endif

#ifndef M_PI
/** Constant \f$\pi\f$ */
#define M_PI (3.14159265358979323844)
#endif

#ifndef M_1_PI
/** Constant \f$1/\pi\f$ */
#define M_1_PI (.31830988618379067154)
#endif

#ifndef M_SQRT1_2
/** Constant \f$1/\sqrt{2}\f$. */
#define M_SQRT1_2 (.70710678118654752440)
#endif

/** The clock_gettime(2) used for internally-managed timing. */
#define TUNA_CLOCK CLOCK_PROCESS_CPUTIME_ID
#ifndef CLOCK_PROCESS_CPUTIME_ID
# error "CLOCK_PROCESS_CPUTIME_ID unavailable"
#endif

size_t
tuna_stats_cnt(const tuna_stats* const stats)
{
    return stats->n;
}

double
tuna_stats_avg(const tuna_stats* const stats)
{
    return stats->n ? stats->m : NAN;
}

double
tuna_stats_var(const tuna_stats* const stats)
{
    return stats->n ? (stats->n > 1 ? stats->s / (stats->n - 1) : 0) : NAN;
}

double
tuna_stats_std(const tuna_stats* const stats)
{
    return sqrt(tuna_stats_var(stats));
}

double
tuna_stats_sum(const tuna_stats* const stats)
{
    return stats->n * stats->m;
}

tuna_stats_mom_result
tuna_stats_mom(const tuna_stats* const stats)
{
    tuna_stats_mom_result result;
    result.n = stats->n;
    switch (stats->n) {
    case 0:
        result.avg = NAN;
        result.var = NAN;
        break;
    case 1:
        result.avg = stats->m;
        result.var = 0;
        break;
    default:
        result.avg = stats->m;
        result.var = stats->s / (stats->n - 1);
        break;
    }
    return result;
}

void
tuna_stats_obs(tuna_stats* const stats,
               const double x)
{
    /* Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.    */
    /* Knuth shows better behavior than Welford 1962 on test data. */
    const size_t n = ++(stats->n);
    if (n > 1) {  /* Second and subsequent invocation */
        double d  = x - stats->m;
        stats->m     += d / n;
        stats->s     += d * (x - stats->m);
    } else {      /* First invocation requires special treatment */
        stats->m = x;
        stats->s = 0;
    }
}

void
tuna_stats_nobs(tuna_stats* const stats,
                const double* x,
                size_t N)
{
    size_t i;
    for (i = N; i -- > 0 ;) {
        tuna_stats_obs(stats, *x++);
    }
}

void
tuna_stats_merge(tuna_stats* const dst,
                 const tuna_stats* const src)
{
    if (dst == src) {          /* merging with self is a NOP */
        /* NOP */
    } else if (src->n == 0) {  /* src contains no data */
        /* NOP */
    } else if (dst->n == 0) {  /* dst contains no data */
        *dst = *src;
    } else {                   /* merge src into dst */
        size_t total = dst->n + src->n;
        double dM    = dst->m - src->m;  /* Cancellation issues? */
        dst->m       = (dst->n * dst->m + src->n * src->m) / total;
        dst->s       = (dst->n == 1 ? 0 : dst->s)
                       + (src->n == 1 ? 0 : src->s)
                       + ((dM * dM) * (dst->n * src->n)) / total;
        dst->n       = total;
    }
}

/**
 * Enforces <code>*a < *b<code> as a postcondition.
 * If doing so required swapping *a and *b, return 1.  Otherwise 0.
 */
static
int
enforce_lt(double* const a, double* const b)
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

void
tuna_chunk_obs(tuna_chunk* const chunk,
               double cost)
{
    /* First, find smallest observation among set {cost, chunk->outliers[0], ... } */
    /* placing it into cost while maintaining sorted-ness of chunk->outliers.      */
    /* The loop is one sort pass with possibility of short-circuiting.      */
    if (enforce_lt(&cost, chunk->outliers)) {
        size_t i;
        for (i = 1;
             i < tuna_countof(chunk->outliers)
             && enforce_lt(chunk->outliers - 1 + i, chunk->outliers + i);
             ++i)
            ;
    }

    /* Second, when non-zero, record statistics about the best observation. */
    if (cost) {
        tuna_stats_obs(&chunk->stats, cost);
    }

    /* Together, the above steps cause a zero-initialized tuna_chunk to */
    /* gather tuna_noutliers pieces of information before beginning to  */
    /* track any statistics.  This effectively provides some "start up" */
    /* or "burn in" period in addition to preventing highly improbable  */
    /* observations from unduly inflating the discovered variability.   */
}

void
tuna_chunk_merge(tuna_stats* const stats,
                 const tuna_chunk* const chunk)
{
    size_t i;
    tuna_stats_merge(stats, &chunk->stats);
    for (i = 0; i < tuna_countof(chunk->outliers); ++i) {
        if (chunk->outliers[i]) {
            tuna_stats_nobs(stats, chunk->outliers + i, tuna_countof(chunk->outliers) - i);
            break;
        }
    }
}

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
static
double
ltqnorm(const double p)
{
    /* Coefficients in rational approximations. */
    static const double a1 = -3.969683028665376e+01;
    static const double a2 =  2.209460984245205e+02;
    static const double a3 = -2.759285104469687e+02;
    static const double a4 =  1.383577518672690e+02;
    static const double a5 = -3.066479806614716e+01;
    static const double a6 =  2.506628277459239e+00;

    static const double b1 = -5.447609879822406e+01;
    static const double b2 =  1.615858368580409e+02;
    static const double b3 = -1.556989798598866e+02;
    static const double b4 =  6.680131188771972e+01;
    static const double b5 = -1.328068155288572e+01;

    static const double c1 = -7.784894002430293e-03;
    static const double c2 = -3.223964580411365e-01;
    static const double c3 = -2.400758277161838e+00;
    static const double c4 = -2.549732539343734e+00;
    static const double c5 =  4.374664141464968e+00;
    static const double c6 =  2.938163982698783e+00;

    static const double d1 =  7.784695709041462e-03;
    static const double d2 =  3.224671290700398e-01;
    static const double d3 =  2.445134137142996e+00;
    static const double d4 =  3.754408661907416e+00;

    /* Define break-points. */
    const double p_low  = 0.02425;
    const double p_high = 1 - p_low;

    /* Compute rational approximation for... */
    double x, q, r, e, u;
    if (p < 0) {                      /* ...domain error. */
        errno = EDOM;
        return 0;
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1572)
#endif
    } else if (p == 0) {              /* ...improbable point. */
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
        errno = ERANGE;
#ifdef HUGE_VAL
        return -HUGE_VAL;
#else
        return -DBL_MAX;
#endif
    } else if (p < p_low) {           /* ...lower region. */
        q = sqrt(-2 * log(p));
        x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    } else if (p <= p_high) {         /* ...central region */
        q = p - 0.5;
        r = q * q;
        x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    } else if (p < 1) {               /* ...upper region. */
        q = sqrt(-2 * log(1 - p));
        x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1572)
#endif
    } else if (p == 1) {              /* ...improbable point. */
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
        errno = ERANGE;
#ifdef HUGE_VAL
        return HUGE_VAL;
#else
        return DBL_MAX;
#endif
    } else {                          /* ...domain error. */
        errno = EDOM;
        return 0;
    }

    /* One hit with Halley's rational method improves to machine precision. */
    e = 0.5 * erfc(x * -M_SQRT1_2) - p;
#define SQRT_2PI (2.50662827463100050240)
    u = e * SQRT_2PI * exp(x * x / 2);
#undef SQRT_2PI
    x = x - u / (1 + x * u / 2);

    return x;
}

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
static
double
as3(double t, int nu)
{
    double r, a, b, c, f, s, fk;
    int i, k, ks, im2, ioe;

    assert(nu > 0); /* Use errno with EDOM instead? */

    f = nu;
    a = t / sqrt(f);
    r = t;
    b = f / (f + r * r);
    im2 = nu - 2;
    ioe = nu % 2;
    s = 1;
    c = 1;
    f = 1;
    ks = ioe + 2;
    fk = ks;
    if (im2 < 2) {
        goto L20;
    }
    i = im2;
    for (k = ks; k <= i; k += 2) {
        c = c * b * (fk - 1) / fk;
        s += c;
        if (s == f) {
            goto L20;
        }
        f = s;
        fk += 2;
    }
L20:
    if (ioe == 1) {
        goto L30;
    }
    return 0.5 + 0.5 * a * sqrt(b) * s;
L30:
    if (nu == 1) {
        s = 0;
    }
    return 0.5 + (a * b * s + atan(a)) * M_1_PI;
}

double
tuna_rand_u01(tuna_state* state)
{
    return rand_r(state) / (double) RAND_MAX;
}

double
tuna_rand_n01(tuna_state* state)
{
    return ltqnorm(tuna_rand_u01(state));
}

tuna_state
tuna_state_default()
{
    const char* d = getenv("TUNA_SEED");
    unsigned int retval;
    if (!d || 1 != sscanf(d, " %u", &retval)) {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        retval = ts.tv_sec + ts.tv_nsec;
    }
    return retval;
}

tuna_welch_result
tuna_welch(double xA, double sA2, size_t nA,
           double xB, double sB2, size_t nB)
{
    tuna_welch_result result;
    double sA2_nA = sA2 / nA;
    double sB2_nB = sB2 / nB;
    double t_den2 = sA2_nA + sB2_nB;
    double nu_den = sA2_nA * sA2_nA / (nA - 1)
                    + sB2_nB * sB2_nB / (nB - 1);
    result.t  = (xA - xB) / sqrt(t_den2);
    result.nu = t_den2 * t_den2 / nu_den;
    return result;
}

double
tuna_welch_t(double xA, double sA2, size_t nA,
             double xB, double sB2, size_t nB)
{
    return (xA - xB) / sqrt(sA2 / nA + sB2 / nB);
}

double
tuna_welch1_nuinf(double xA, double sA2, size_t nA,
                  double xB, double sB2, size_t nB)
{
    return 1 - erfc(-tuna_welch_t(xA, sA2, nA, xB, sB2, nB) * M_SQRT1_2) / 2;
}

double
tuna_welch1_approx(double xA, double sA2, size_t nA,
                   double xB, double sB2, size_t nB)
{
    tuna_welch_result welch = tuna_welch(xA, sA2, nA, xB, sB2, nB);
    double t = welch.t;
    if (welch.nu > 2) {
        t *= welch.nu / (welch.nu - 2);
    }
    return 1 - erfc(-t * M_SQRT1_2) / 2;
}

double
tuna_welch1(double xA, double sA2, size_t nA,
            double xB, double sB2, size_t nB)
{
    tuna_welch_result welch = tuna_welch(xA, sA2, nA, xB, sB2, nB);
    return 1 - as3(welch.t, welch.nu);
}

/* http://agentzlerich.blogspot.com/2011/01/c-header-only-unit-testing-with-fctx.html */
static
void
trim(char* const a)
{
    char* p = a, *q = a;
    while (isspace(*q)) {
        ++q;
    }
    while (*q) {
        *p++ = *q++;
    }
    *p = '\0';
    while (p > a && isspace(*--p)) {
        *p = '\0';
    }
}


/* Very similar to tuna_algo_welch1_nuinf(...) */
static size_t
tuna_algo_welch1_nuinf_impl(const size_t nchunk,
                            const tuna_chunk* chunks,
                            const double* u01)
{
    size_t i, j;
    tuna_stats_mom_result istats, jstats;
    double p;
    const size_t ncomparisons = nchunk - 1;

    assert(nchunk > 0);
    i = 0;
    istats = tuna_stats_mom(&chunks[0].stats);
    if (istats.n < 2) {
        return i;
    }
    for (j = 1; j < nchunk; ++j) {
        jstats = tuna_stats_mom(&chunks[j].stats);
        if (jstats.n < 2) {
            return j;
        }
        p = tuna_welch1_nuinf(istats.avg, istats.var, istats.n,
                              jstats.avg, jstats.var, jstats.n);
        /* Apply Bonferroni correction for multiple comparisons writing   */
        /* p * ncomparisons < u01[j] instead of p < u01[j] / ncomparisons */
        if (p * ncomparisons < u01[j]) {
            i = j;
            istats = jstats;
        }
    }
    return i;
}

/* Very similar to tuna_algo_welch1(...) */
static size_t
tuna_algo_welch1_impl(const size_t nchunk,
                      const tuna_chunk* chunks,
                      const double* u01)
{
    size_t i, j;
    tuna_stats_mom_result istats, jstats;
    double p;
    const size_t ncomparisons = nchunk - 1;

    assert(nchunk > 0);
    i = 0;
    istats = tuna_stats_mom(&chunks[0].stats);
    if (istats.n < 2) {
        return i;
    }
    for (j = 1; j < nchunk; ++j) {
        jstats = tuna_stats_mom(&chunks[j].stats);
        if (jstats.n < 2) {
            return j;
        }
        p = tuna_welch1(istats.avg, istats.var, istats.n,
                        jstats.avg, jstats.var, jstats.n);
        /* Apply Bonferroni correction for multiple comparisons writing   */
        /* p * ncomparisons < u01[j] instead of p < u01[j] / ncomparisons */
        if (p * ncomparisons < u01[j]) {
            i = j;
            istats = jstats;
        }
    }
    return i;
}

static size_t
tuna_algo_zero_impl(const size_t nchunk,
                    const tuna_chunk* chunks,
                    const double* u01)
{
    (void) nchunk;
    (void) chunks;
    (void) u01;
    return 0;
}

static size_t
tuna_algo_uniform_impl(const size_t nchunk,
                       const tuna_chunk* chunks,
                       const double* u01)
{
    (void) chunks;
    assert(nchunk > 0);
    return (size_t)(u01[0] * nchunk);
}

/* Algorithm instances carrying their own names */
static const tuna_algo tuna_algo_welch1_nuinf_s = {
    "welch1_nuinf",
    tuna_algo_welch1_nuinf_impl
};

static const tuna_algo tuna_algo_welch1_s = {
    "welch1",
    tuna_algo_welch1_impl
};

static const tuna_algo tuna_algo_zero_s = {
    "zero",
    tuna_algo_zero_impl
};

static const tuna_algo tuna_algo_uniform_s = {
    "uniform",
    tuna_algo_uniform_impl
};

/* Public algorithm handles */
const tuna_algo* const tuna_algo_welch1_nuinf = &tuna_algo_welch1_nuinf_s;
const tuna_algo* const tuna_algo_welch1       = &tuna_algo_welch1_s;
const tuna_algo* const tuna_algo_zero         = &tuna_algo_zero_s;
const tuna_algo* const tuna_algo_uniform      = &tuna_algo_uniform_s;

/* Registry of all known algorithms */
static const tuna_algo* const known_algos[] = {
    &tuna_algo_welch1_s,
    &tuna_algo_welch1_nuinf_s,
    &tuna_algo_zero_s,
    &tuna_algo_uniform_s
};

const char*
tuna_algo_name(const tuna_algo* algo)
{
    return algo ? algo->name : "unknown";
}

const tuna_algo*
tuna_algo_default(const size_t nchunk)
{
    char* d;
    size_t i;

    if (nchunk < 2) {
        return tuna_algo_zero;
    }
    d = getenv("TUNA_ALGO");
    if (d) {
        trim(d);
        for (i = 0; i < tuna_countof(known_algos); ++i) {
            if (!strcasecmp(known_algos[i]->name, d)) {
                return known_algos[i];
            }
        }
    }
    return tuna_algo_welch1; /* Default */
}

/* TODO Do something intelligent with clock_getres(2) information */
size_t
tuna_pre_cost(tuna_site* site,
              tuna_stack* stack,
              const tuna_chunk* chunks,
              const size_t nchunk)
{
    size_t i;
    double* u01;

    /* Ensure a zero-initialize stack argument produces good behavior by... */
    if (!site->algo) {
        /* ...providing a default algorithm when not set, and... */
        site->algo = tuna_algo_default(nchunk);

        if (!site->state) {
            /* ...providing a default seed when not set. */
            site->state = tuna_state_default();
        }
    }

    /* Prepare nchunk random numbers for use by the algorithm.        */
    /* Drawing variates here avoids state updates from algorithms. */
    u01 = __builtin_alloca(nchunk * sizeof(double));
    for (i = 0; i < nchunk; ++i) {
        u01[i] = tuna_rand_u01(&site->state);
    }

    /* Invoke chosen algorithm saving selected index for tuna_post_cost(). */
    stack->ik = site->algo->function(nchunk, chunks, u01);

    return stack->ik;
}

void
tuna_post_cost(const tuna_stack* stack,
               tuna_chunk* chunks,
               const double cost)
{
    tuna_chunk_obs(chunks + stack->ik, cost);
}

/**
 * A minimal, standard-compliant <code>struct timespec</code> to provide sane
 * alignment-aware <code>sizeof</code> behavior for the <code>memcpy(3)</code>
 * operation in \ref tuna_pre() and \ref tuna_post().  Required as <a
 * href="http://pubs.opengroup.org/onlinepubs/009695299/basedefs/time.h.html">man
 * time.h</a> states "The <time.h> header shall declare the structure timespec,
 * which has at least the following members..." so <code>sizeof(struct
 * timespec)</code> may vary in undesired ways on conforming systems.
 */
struct tuna_timespec_minimal {
    time_t tv_sec;
    long   tv_nsec;
};

size_t
tuna_pre(tuna_site* site,
         tuna_stack* stack,
         const tuna_chunk* chunks,
         const size_t nchunk)
{
    struct timespec ts;

    /* Delegate chunk selection to tuna_pre_cost */
    tuna_pre_cost(site, stack, chunks, nchunk);

    /* Glimpse at the clock so we may compute elapsed time in tuna_post(). */
    /* Done after algorithm performed to avoid accumulating its runtime.   */
    clock_gettime(TUNA_CLOCK, &ts);

    /* Stash the time in the opaque stack->ts buffer            */
    /* (after being certain the buffer is adequately sized). */
    tuna_assert_static(sizeof(stack->ts) >= sizeof(struct tuna_timespec_minimal));
    memcpy(stack->ts, &ts, sizeof(struct tuna_timespec_minimal));

    return stack->ik;
}

double
tuna_post(const tuna_stack* stack,
          tuna_chunk* chunks)
{
    struct timespec ts, te;
    double elapsed;

    /* Glimpse at the clock before bookkeeping to ignore its overhead. */
    clock_gettime(TUNA_CLOCK, &te);

    /* Retrieve tuna_pre time from the opaque stack->ts buffer  */
    /* (after being certain the buffer is adequately sized). */
    tuna_assert_static(sizeof(stack->ts) >= sizeof(struct tuna_timespec_minimal));
    memcpy(&ts, stack->ts, sizeof(struct tuna_timespec_minimal));

    /* Compute double-valued elapsed time */
    elapsed  = te.tv_nsec - ts.tv_nsec;  /* Nanoseconds...  */
    elapsed *= 1e-9;                     /* ...to seconds   */
    elapsed += te.tv_sec  - ts.tv_sec;   /* ...plus seconds */

    /* Delegate recording the observation to tuna_post_cost() */
    tuna_post_cost(stack, chunks, elapsed);

    return elapsed;
}

int
tuna_fprint(FILE* stream,
            const tuna_site* site,
            const tuna_chunk* chunks,
            const size_t nchunk,
            const char* prefix,
            const char** labels)
{
    size_t ik;
    int nwritten, namelen, status;

    /* Output a bash-like "TUNA$" prompt identifying this tuning site. */
    nwritten = tuna_site_fprintf(stream, site, "TUNA$ %s", prefix);

    /* Find the maximum label length so post-label outputs may be aligned */
    namelen = sizeof("chunk") + 5;
    if (labels) {
        for (ik = 0; ik < nchunk; ++ik) {
            if (labels[ik] && *labels[ik]) {
                int len = strlen(labels[ik]);
                if (len > namelen) {
                    namelen = len;
                }
            }
        }
    }

    /* Output continuation-like "TUNA>", the site, labels, and statistics. */
    for (ik = 0; ik < nchunk && nwritten >= 0; ++ik) {
        if (labels && labels[ik] && *labels[ik]) {
            status = tuna_chunk_fprintf(stream, chunks + ik, "TUNA> %s %-*s",
                                        prefix, namelen, labels[ik]);
        } else {
            status = tuna_chunk_fprintf(stream, chunks + ik, "TUNA> %s chunk%0*zu",
                                        prefix, namelen - sizeof("chunk"), ik);
        }
        nwritten = status >= 0
                   ? nwritten + status
                   : status;
    }
    return nwritten;
}

int
tuna_chunk_vfprintf(FILE* stream,
                    const tuna_chunk* chunk,
                    const char* format,
                    va_list ap)
{
    int nwritten;
    nwritten = vfprintf(stream, format, ap);
    if (nwritten >= 0) {
        int status;
        tuna_stats o;
        tuna_stats_mom_result stats;
        memset(&o, 0, sizeof(o));
        tuna_chunk_merge(&o, chunk);
        stats = tuna_stats_mom(&o);
        status = fprintf(stream,
                         "%s%*lu  %-#*.*g +/- %-#*.*g\n",
                         format[0] ? " " : "",
                         FLT_DIG + 4, (long unsigned) stats.n,
                         FLT_DIG + 4, FLT_DIG, stats.avg,
                         FLT_DIG + 4, FLT_DIG, sqrt(stats.var));
        nwritten = status >= 0
                   ? nwritten + status
                   : status;
    }
    return nwritten;
}

int
tuna_chunk_fprintf(FILE* stream,
                   const tuna_chunk* chunk,
                   const char* format,
                   ...)
{
    int nwritten;
    va_list ap;
    va_start(ap, format);
    nwritten = tuna_chunk_vfprintf(stream, chunk, format, ap);
    va_end(ap);
    return nwritten;
}

int
tuna_site_vfprintf(FILE* stream,
                   const tuna_site* site,
                   const char* format,
                   va_list ap)
{
    const char* name;
    int nwritten;
    nwritten = vfprintf(stream, format, ap);
    name = tuna_algo_name(site->algo);
    if (nwritten >= 0) {
        int status = fprintf(stream,
                             "%s"
                             "%s\n",
                             format[0] ? " " : "",
                             name);
        nwritten = status >= 0
                   ? nwritten + status
                   : status;
    }
    return nwritten;
}

int
tuna_site_fprintf(FILE* stream,
                  const tuna_site* site,
                  const char* format,
                  ...)
{
    int nwritten;
    va_list ap;
    va_start(ap, format);
    nwritten = tuna_site_vfprintf(stream, site, format, ap);
    va_end(ap);
    return nwritten;
}

/******************************************************/
/* Below here is scratch storage for half-baked ideas */
/******************************************************/

/* TODO Document.                                   */
/* TODO Note that nodes own their \c id strings.    */
/* TODO Add enough additional members to be useful. */
typedef struct tuna_registry {
    struct tuna_registry* left;  /**< Left subtree which may be \c NULL.    */
    struct tuna_registry* right; /**< Right subtree which may be \c NULL.   */
    const char id[1];            /**< Sorting key owned via flexible struct */
} tuna_registry;

/* TODO Document. */
tuna_registry*
tuna_registry_alloc(const char* id)
{
    /* Struct hack using id[1] already includes space for NULL terminator. */
    /* Using calloc(3) sets left = right = NULL and enforces termination.  */
    const size_t n = strlen(id);
    tuna_registry* p = calloc(n + sizeof(tuna_registry), 1);
    if (p) {
        memcpy((void*) p->id, (void*) id, n);
    }
    return p;
}

/* TODO Document. */
void
tuna_registry_free(tuna_registry* n)
{
    /* Postorder in which struct hack for id implies just one free(3) call. */
    if (n) {
        tuna_registry_free(n->left);
        tuna_registry_free(n->right);
        free(n);
    }
}

/* TODO Document. */
tuna_registry*
tuna_registry_insert(tuna_registry** n, const char* id)
{
    if (*n) {
        const int cmp = strcmp(id, (*n)->id);
        if (cmp < 0) {
            return tuna_registry_insert(&(*n)->left, id);
        } else if (cmp == 0) {
            return *n;
        } else {
            return tuna_registry_insert(&(*n)->right, id);
        }
    } else {
        return *n = tuna_registry_alloc(id);
    }
}

/* TODO Document. */
tuna_registry*
tuna_registry_find(tuna_registry* n, const char* id)
{
    if (n) {
        const int cmp = strcmp(id, n->id);
        if (cmp < 0) {
            return tuna_registry_find(n->left, id);
        } else if (cmp == 0) {
            return n;
        } else {
            return tuna_registry_find(n->right, id);
        }
    } else {
        return NULL;
    }
}

/** Generate a preorder accumulator producing \c acc named \c func. */
#define TRAVERSE_PREORDER(func, node, acc)                           \
    acc func(node* n, acc (*v)(node*, acc), acc a)                   \
    { return !n ? a : func(n->right, v, func(n->left, v, v(n, a))); }

/** Generate an inorder accumulator producing \c acc named \c func. */
#define TRAVERSE_INORDER(func, node, acc)                            \
    acc func(node* n, acc (*v)(node*, acc), acc a)                   \
    { return !n ? a : func(n->right, v, v(n, func(n->left, v, a))); }

/** Generate a postorder accumulator producing \c acc named \c func. */
#define TRAVERSE_POSTORDER(func, node, acc)                           \
    acc func(node* n, acc (*v)(node*, acc), acc a)                    \
    { return !n ? a : v(n, func(n->right, v, func(n->left, v, a))); }

/* TODO Document. */

TRAVERSE_PREORDER(tuna_registry_preorder_double,   tuna_registry, double)
TRAVERSE_PREORDER(tuna_registry_preorder_int,      tuna_registry, int)
TRAVERSE_PREORDER(tuna_registry_preorder_size_t,   tuna_registry, size_t)
TRAVERSE_PREORDER(tuna_registry_preorder_voidp,    tuna_registry, void*)
TRAVERSE_INORDER(tuna_registry_inorder_double,     tuna_registry, double)
TRAVERSE_INORDER(tuna_registry_inorder_int,        tuna_registry, int)
TRAVERSE_INORDER(tuna_registry_inorder_size_t,     tuna_registry, size_t)
TRAVERSE_INORDER(tuna_registry_inorder_voidp,      tuna_registry, void*)
TRAVERSE_POSTORDER(tuna_registry_postorder_double, tuna_registry, double)
TRAVERSE_POSTORDER(tuna_registry_postorder_int,    tuna_registry, int)
TRAVERSE_POSTORDER(tuna_registry_postorder_size_t, tuna_registry, size_t)
TRAVERSE_POSTORDER(tuna_registry_postorder_voidp,  tuna_registry, void*)
