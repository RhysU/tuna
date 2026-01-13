/*
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
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
#include <tuna.h>

#if (_POSIX_C_SOURCE < 200112L)
#define _POSIX_C_SOURCE (200112L)
#endif

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

#define STRINGIFY(s) STRINGIFY_HELPER(s)
#define STRINGIFY_HELPER(s) #s

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
tuna_stats_cnt(const tuna_stats* const t)
{
    return t->n;
}

double
tuna_stats_avg(const tuna_stats* const t)
{
    return t->n ? t->m : NAN;
}

double
tuna_stats_var(const tuna_stats* const t)
{
    return t->n ? (t->n > 1 ? t->s / (t->n - 1) : 0) : NAN;
}

double
tuna_stats_std(const tuna_stats* const t)
{
    return sqrt(tuna_stats_var(t));
}

double
tuna_stats_sum(const tuna_stats* const t)
{
    return t->n * t->m;
}

size_t
tuna_stats_mom(const tuna_stats* const t,
               double* const avg,
               double* const var)
{
    switch (t->n) {
    case 0:
        *avg = NAN;
        *var = NAN;
        break;
    case 1:
        *avg = t->m;
        *var = 0;
        break;
    default:
        *avg = t->m;
        *var = t->s / (t->n - 1);
        break;
    }
    return t->n;
}

tuna_stats*
tuna_stats_obs(tuna_stats* const t,
               const double x)
{
    /* Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.    */
    /* Knuth shows better behavior than Welford 1962 on test data. */
    const size_t n = ++(t->n);
    if (n > 1) {  /* Second and subsequent invocation */
        double d  = x - t->m;
        t->m     += d / n;
        t->s     += d * (x - t->m);
    } else {      /* First invocation requires special treatment */
        t->m = x;
        t->s = 0;
    }
    return t;
}

tuna_stats*
tuna_stats_nobs(tuna_stats* const t,
                const double* x,
                size_t N)
{
    size_t i;
    for (i = N; i -- > 0 ;) {
        tuna_stats_obs(t, *x++);
    }
    return t;
}

tuna_stats*
tuna_stats_merge(tuna_stats* const dst,
                 const tuna_stats* const src)
{
    if (src->n == 0) {         /* snp contains no data */
        /* NOP */
    } else if (dst->n == 0) {  /* dst contains no data */
        *dst = *src;
    } else {                   /* merge snp into dst */
        size_t total = dst->n + src->n;
        double dM    = dst->m - src->m;  /* Cancellation issues? */
        dst->m       = (dst->n * dst->m + src->n * src->m) / total;
        dst->s       = (dst->n == 1 ? 0 : dst->s)
                       + (src->n == 1 ? 0 : src->s)
                       + ((dM * dM) * (dst->n * src->n)) / total;
        dst->n       = total;
    }
    return dst;
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

tuna_chunk*
tuna_chunk_obs(tuna_chunk* const k,
               double t)
{
    /* First, find smallest observation among set {t, k->outliers[0], ... } */
    /* placing it into t while maintaining sorted-ness of k->outliers.      */
    /* The loop is one sort pass with possibility of short-circuiting.      */
    if (enforce_lt(&t, k->outliers)) {
        size_t i;
        for (i = 1;
             i < sizeof(k->outliers) / sizeof(k->outliers[0])
             && enforce_lt(k->outliers - 1 + i, k->outliers + i);
             ++i)
            ;
    }

    /* Second, when non-zero, record statistics about the best observation. */
    /* Use the internal method to avoid deadlock as we already hold lock.   */
    if (t) {
        tuna_stats_obs(&k->stats, t);
    }

    /* Together, the above steps cause a zero-initialized tuna_chunk to */
    /* gather tuna_noutliers pieces of information before beginning to  */
    /* track any statistics.  This effectively provides some "start up" */
    /* or "burn in" period in addition to preventing highly improbable  */
    /* observations from unduly inflating the discovered variability.   */

    return k;
}

tuna_stats*
tuna_chunk_merge(tuna_stats* const s,
                 const tuna_chunk* const k)
{
    size_t i;
    tuna_stats_merge(s, &k->stats);
    for (i = 0; i < tuna_countof(k->outliers); ++i) {
        if (k->outliers[i]) {
            tuna_stats_nobs(s, k->outliers + i, tuna_countof(k->outliers) - i);
            break;
        }
    }
    return s;
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
tuna_rand_u01(tuna_seed* sd)
{
    return rand_r(sd) / (double) RAND_MAX;
}

double
tuna_rand_n01(tuna_seed* sd)
{
    return ltqnorm(tuna_rand_u01(sd));
}

tuna_seed
tuna_seed_default()
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

void
tuna_welch(double xA, double sA2, size_t nA,
           double xB, double sB2, size_t nB,
           double* const t, double* const nu)
{
    double sA2_nA = sA2 / nA;
    double sB2_nB = sB2 / nB;
    double t_den2 = sA2_nA + sB2_nB;
    double nu_den = sA2_nA * sA2_nA / (nA - 1)
                    + sB2_nB * sB2_nB / (nB - 1);
    *t            = (xA - xB) / sqrt(t_den2);
    *nu           = t_den2 * t_den2 / nu_den;
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
    double t, nu;
    tuna_welch(xA, sA2, nA, xB, sB2, nB, &t, &nu);
    if (nu > 2) {
        t *= nu / (nu - 2);
    }
    return 1 - erfc(-t * M_SQRT1_2) / 2;
}

double
tuna_welch1(double xA, double sA2, size_t nA,
            double xB, double sB2, size_t nB)
{
    double t, nu;
    tuna_welch(xA, sA2, nA, xB, sB2, nB, &t, &nu);
    return 1 - as3(t, nu);
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

/* TODO Should algo be a struct with a pointer to its own name? */

/* Should be kept in sync with tuna_algo_default just below */
const char*
tuna_algo_name(tuna_algo al)
{
    if (al == &tuna_algo_welch1) {
        return "welch1";
    } else if (al == &tuna_algo_welch1_nuinf) {
        return "welch1_nuinf";
    } else if (al == &tuna_algo_zero) {
        return "zero";
    } else {
        return "unknown";
    }
}

/* Should be kept in sync with tuna_algo_name just above */
tuna_algo
tuna_algo_default(const size_t nk)
{
    char* d;
    if (nk < 2) {
        return &tuna_algo_zero;
    }
    d = getenv("TUNA_ALGO");
    if (d) {
        trim(d);
        if (!strcasecmp("welch1", d)) {
            return &tuna_algo_welch1;
        } else if (!strcasecmp("welch1_nuinf", d)) {
            return &tuna_algo_welch1_nuinf;
        } else if (!strcasecmp("zero", d)) {
            return &tuna_algo_zero;
        }
    }
    return &tuna_algo_welch1; /* Default */
}

size_t
tuna_algo_welch1_nuinf(const size_t nk,
                       const tuna_chunk ks[],
                       const double u01[])
{
    size_t i, j;
    size_t icnt, jcnt;
    double iavg, ivar, javg, jvar, p;

    assert(nk > 0);
    i = 0;
    icnt = tuna_stats_mom(&ks[0].stats, &iavg, &ivar);
    if (icnt < 2) {
        return i;
    }
    for (j = 1; j < nk; ++j) {
        jcnt = tuna_stats_mom(&ks[j].stats, &javg, &jvar);
        if (jcnt < 2) {
            return j;
        }
        p    = tuna_welch1_nuinf(iavg, ivar, icnt, javg, jvar, jcnt);
        if (p < u01[j]) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

size_t
tuna_algo_welch1(const size_t nk,
                 const tuna_chunk ks[],
                 const double u01[])
{
    size_t i, j;
    size_t icnt, jcnt;
    double iavg, ivar, javg, jvar, p;

    assert(nk > 0);
    i = 0;
    icnt = tuna_stats_mom(&ks[0].stats, &iavg, &ivar);
    if (icnt < 2) {
        return i;
    }
    for (j = 1; j < nk; ++j) {
        jcnt = tuna_stats_mom(&ks[j].stats, &javg, &jvar);
        if (jcnt < 2) {
            return j;
        }
        p    = tuna_welch1(iavg, ivar, icnt, javg, jvar, jcnt);
        if (p < u01[j]) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

size_t
tuna_algo_zero(const size_t nk,
               const tuna_chunk ks[],
               const double u01[])
{
    (void) nk;
    (void) ks;
    (void) u01;
    return 0;
}

/* TODO Do something intelligent with clock_getres(2) information */
size_t
tuna_pre_cost(tuna_site* si,
              tuna_stack* st,
              const tuna_chunk ks[],
              const size_t nk)
{
    size_t i;
    double* u01;

    /* Ensure a zero-initialize st argument produces good behavior by... */
    if (!si->al) {
        /* ...providing a default algorithm when not set, and... */
        si->al = tuna_algo_default(nk);

        if (!si->sd) {
            /* ...providing a default seed when not set. */
            si->sd = tuna_seed_default();
        }
    }

    /* Prepare nk random numbers for use by the algorithm.        */
    /* Drawing variates here avoids seed updates from algorithms. */
    /* That permits a short critical section on the tuna_site.    */
    u01 = __builtin_alloca(nk * sizeof(double));
    for (i = 0; i < nk; ++i) {
        u01[i] = tuna_rand_u01(&si->sd);
    }

    /* Invoke chosen algorithm saving selected index for tuna_post_cost(). */
    st->ik = si->al(nk, ks, u01);

    return st->ik;
}

void
tuna_post_cost(const tuna_stack* st,
               tuna_chunk ks[],
               const double cost)
{
    tuna_chunk_obs(ks + st->ik, cost);
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
tuna_pre(tuna_site* si,
         tuna_stack* st,
         const tuna_chunk ks[],
         const size_t nk)
{
    struct timespec ts;

    /* Delegate chunk selection to tuna_pre_cost */
    tuna_pre_cost(si, st, ks, nk);

    /* Glimpse at the clock so we may compute elapsed time in tuna_post(). */
    /* Done after algorithm performed to avoid accumulating its runtime.   */
    clock_gettime(TUNA_CLOCK, &ts);

    /* Stash the time in the opaque st->ts buffer            */
    /* (after being certain the buffer is adequately sized). */
    tuna_assert_static(sizeof(st->ts) >= sizeof(struct tuna_timespec_minimal));
    memcpy(st->ts, &ts, sizeof(struct tuna_timespec_minimal));

    return st->ik;
}

double
tuna_post(const tuna_stack* st,
          tuna_chunk ks[])
{
    struct timespec ts, te;
    double elapsed;

    /* Glimpse at the clock before bookkeeping to ignore its overhead. */
    clock_gettime(TUNA_CLOCK, &te);

    /* Retrieve tuna_pre time from the opaque st->ts buffer  */
    /* (after being certain the buffer is adequately sized). */
    tuna_assert_static(sizeof(st->ts) >= sizeof(struct tuna_timespec_minimal));
    memcpy(&ts, st->ts, sizeof(struct tuna_timespec_minimal));

    /* Compute double-valued elapsed time */
    elapsed  = te.tv_nsec - ts.tv_nsec;  /* Nanoseconds...  */
    elapsed *= 1e-9;                     /* ...to seconds   */
    elapsed += te.tv_sec  - ts.tv_sec;   /* ...plus seconds */

    /* Delegate recording the observation to tuna_post_cost() */
    tuna_post_cost(st, ks, elapsed);

    return elapsed;
}

int
tuna_chunk_fprint(void* stream,
                  const tuna_chunk* k,
                  const char* prefix)
{
    return tuna_chunk_fprintf(stream, k, prefix);
}

int
tuna_site_fprint(void* stream,
                 const tuna_site* si,
                 const char* prefix)
{
    return tuna_site_fprintf(stream, si, prefix);
}

int
tuna_fprint(void* stream,
            const tuna_site* si,
            const tuna_chunk ks[],
            const size_t nk,
            const char* prefix,
            const char* labels[])
{
    size_t ik;
    int nwritten, namelen, status;

    /* Output a bash-like "TUNA$" prompt identifying this tuning site. */
    nwritten = tuna_site_fprintf(stream, si, "TUNA$ %s", prefix);

    /* Find the maximum label length so post-label outputs may be aligned */
    namelen = sizeof("chunk") + 5;
    if (labels) {
        for (ik = 0; ik < nk; ++ik) {
            if (labels[ik] && *labels[ik]) {
                int len = strlen(labels[ik]);
                if (len > namelen) {
                    namelen = len;
                }
            }
        }
    }

    /* Output continuation-like "TUNA>", the site, labels, and statistics. */
    for (ik = 0; ik < nk && nwritten >= 0; ++ik) {
        if (labels && labels[ik] && *labels[ik]) {
            status = tuna_chunk_fprintf(stream, ks + ik, "TUNA> %s %-*s",
                                        prefix, namelen, labels[ik]);
        } else {
            status = tuna_chunk_fprintf(stream, ks + ik, "TUNA> %s chunk%0*zu",
                                        prefix, namelen - sizeof("chunk"), ik);
        }
        nwritten = status >= 0
                   ? nwritten + status
                   : status;
    }
    return nwritten;
}

int
tuna_chunk_fprintf(void* stream,
                   const tuna_chunk* k,
                   const char* format,
                   ...)
{
    int nwritten;
    va_list ap;
    va_start(ap, format);
    nwritten = vfprintf(stream, format, ap);
    va_end(ap);
    if (nwritten >= 0) {
        int status;
        size_t cnt;
        double avg, var;
        tuna_stats o;
        memset(&o, 0, sizeof(o));
        tuna_chunk_merge(&o, k);
        cnt = tuna_stats_mom(&o, &avg, &var);
        status = fprintf(stream,
                         "%s"
                         "%"  STRINGIFY(DBL_DIG) "lu "
                         "%#" STRINGIFY(DBL_DIG) "g +/- "
                         "%#" STRINGIFY(DBL_DIG) "g\n",
                         format[0] ? " " : "",
                         (long unsigned) cnt, avg, sqrt(var));
        nwritten = status >= 0
                   ? nwritten + status
                   : status;
    }
    return nwritten;
}

int
tuna_site_fprintf(void* stream,
                  const tuna_site* si,
                  const char* format,
                  ...)
{
    const char* name;
    int nwritten;
    va_list ap;
    va_start(ap, format);
    nwritten = vfprintf(stream, format, ap);
    va_end(ap);
    name = tuna_algo_name(si->al);
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
tuna_registry_alloc(const char id[])
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
tuna_registry_insert(tuna_registry** n, const char id[])
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
tuna_registry_find(tuna_registry* n, const char id[])
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
