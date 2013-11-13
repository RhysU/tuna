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
 * As library eschews features beyond C89, having a single compilation unit may
 * be important to permit effective interprocedural optimization.  It also
 * simplifies adoption by other projects as it may simply be copied into an
 * existing build tree.
 */

#include <tuna.h>

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
tuna_stats_fastavg(const tuna_stats* const t)
{
    assert(tuna_stats_cnt(t) > 0);
    return t->m;
}

double
tuna_stats_var(const tuna_stats* const t)
{
    return t->n ? (t->n > 1 ? t->s / (t->n - 1) : 0) : NAN;
}

double
tuna_stats_fastvar(const tuna_stats* const t)
{
    assert(tuna_stats_cnt(t) > 1);
    return t->s / (t->n - 1);
}

double
tuna_stats_std(const tuna_stats* const t)
{
    return sqrt(tuna_stats_var(t));
}

double
tuna_stats_sum(const tuna_stats* const t)
{
    return tuna_stats_cnt(t) * tuna_stats_avg(t);
}

double
tuna_stats_faststd(const tuna_stats* const t)
{
    assert(tuna_stats_cnt(t) > 1);
    return sqrt(tuna_stats_fastvar(t));
}

tuna_stats*
tuna_stats_fastobs(tuna_stats* const t,
                   const double x)
{
    // Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.
    // Knuth shows better behavior than Welford 1962 on test data.
    assert(tuna_stats_cnt(t) > 0);
    size_t n  = ++(t->n);
    double d  = x - t->m;
    t->m     += d / n;
    t->s     += d * (x - t->m);
    return t;
}

tuna_stats*
tuna_stats_obs(tuna_stats* const t,
               const double x)
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

tuna_stats*
tuna_stats_nobs(tuna_stats* const t,
                const double* x,
                size_t N)
{
    size_t i;
    if (N) {                               // NOP on degenerate input
        tuna_stats_obs(t, *x++);           // Delegate possible n == 1
        for (i = --N; i -- > 0 ;) {        // Henceforth, certainly n > 1
            tuna_stats_fastobs(t, *x++);
        }
    }
    return t;
}

tuna_stats*
tuna_stats_merge(tuna_stats* const dst,
                 const tuna_stats* const src)
{
    if (src->n == 0) {         // src contains no data
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

/**
 * Enforces <code>*a < *b<code> as a postcondition.
 * If doing so required swapping *a and *b, return 1.  Otherwise 0.
 */
static inline
int enforce_lt(double* const a, double* const b)
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
    // First, find smallest observation among set {t, k->outliers[0], ... }
    // placing it into storage t while maintaining sorted-ness of k->outliers.
    // The loop is one bubble sort pass with possibility of short-circuiting.
    if (enforce_lt(&t, k->outliers)) {
        size_t i;
        for (i = 1;
             i < sizeof(k->outliers) / sizeof(k->outliers[0])
             && enforce_lt(k->outliers - 1 + i, k->outliers + i);
             ++i)
            ;
    }

    // Second, when non-zero, record statistics about the best observation.
    if (t) {
        tuna_stats_obs(&k->stats, t);
    }

    // Together, these two steps cause a zero-initialized tuna_chunk to
    // gather tuna_noutliers pieces of information before beginning to
    // track any statistics.  This effectively provides some "start up"
    // or "burn in" period in addition to preventing highly improbable
    // observations from unduly inflating the discovered variability.

    return k;
}

tuna_stats
tuna_chunk_stats(tuna_chunk* const k)
{
    tuna_stats s = k->stats;
    return s;
}

tuna_stats*
tuna_chunk_merge(tuna_stats* const s,
                  const tuna_chunk* const k)
{
    size_t i;
    tuna_stats_merge(s, &k->stats);
    for (i = 0; i < tuna_countof(k->outliers); ++i) {
        if (k->outliers[i]) {
            return tuna_stats_nobs(s, k->outliers + i,
                                   tuna_countof(k->outliers) - i);
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
    // Coefficients in rational approximations.
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

    // Define break-points.
    const double p_low  = 0.02425;
    const double p_high = 1 - p_low;

    // Compute rational approximation for...
    double x, q, r;
    if (p < 0) {                      // ...domain error.
        errno = EDOM;
        return 0;
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1572)
#endif
    } else if (p == 0) {              // ...improbable point.
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
        errno = ERANGE;
#ifdef HUGE_VAL
        return -HUGE_VAL;
#else
        return -DBL_MAX;
#endif
    } else if (p < p_low) {           // ...lower region.
        q = sqrt(-2 * log(p));
        x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    } else if (p <= p_high) {         // ...central region
        q = p - 0.5;
        r = q * q;
        x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    } else if (p < 1) {               // ...upper region.
        q = sqrt(-2 * log(1 - p));
        x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1572)
#endif
    } else if (p == 1) {              // ...improbable point.
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
        errno = ERANGE;
#ifdef HUGE_VAL
        return HUGE_VAL;
#else
        return DBL_MAX;
#endif
    } else {                          // ...domain error.
        errno = EDOM;
        return 0;
    }

    // One iteration of Halley's rational method improves to machine precision.
    const double sqrt_2pi = sqrt(2 * M_PI);
    double e, u;
    e = 0.5 * erfc(x * -M_SQRT1_2) - p;
    u = e * sqrt_2pi * exp(x * x / 2);
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
    assert(nu > 0); // Use errno with EDOM instead?

    double r, a, b, c, f, s, fk;
    int i, k, ks, im2, ioe;

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
    *t            = (xA - xB) / sqrt(t_den2);
    double nu_den = sA2_nA * sA2_nA / (nA - 1)
                    + sB2_nB * sB2_nB / (nB - 1);
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

// http://agentzlerich.blogspot.com/2011/01/c-header-only-unit-testing-with-fctx.html
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

tuna_algo
tuna_algo_default(const int nk)
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
    return &tuna_algo_welch1; // Default
}

int
tuna_algo_welch1_nuinf(const int nk,
                       const tuna_chunk* ks,
                       tuna_seed* seed)
{
    int i, j;
    double icnt, iavg, ivar, jcnt, javg, jvar, p;

    assert(nk > 0);
    i = 0;
    icnt = tuna_stats_cnt(&ks[0].stats);
    if (icnt < 2) {
        return i;
    }
    iavg = tuna_stats_fastavg(&ks[0].stats);
    ivar = tuna_stats_fastvar(&ks[0].stats);
    for (j = 1; j < nk; ++j) {
        jcnt = tuna_stats_cnt(&ks[j].stats);
        if (jcnt < 2) {
            return j;
        }
        javg = tuna_stats_fastavg(&ks[j].stats);
        jvar = tuna_stats_fastvar(&ks[j].stats);
        p    = tuna_welch1_nuinf(iavg, ivar, icnt, javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int
tuna_algo_welch1(const int nk,
                 const tuna_chunk* ks,
                 tuna_seed* seed)
{
    int i, j;
    double icnt, iavg, ivar, jcnt, javg, jvar, p;

    assert(nk > 0);
    i = 0;
    icnt = tuna_stats_cnt(&ks[0].stats);
    if (icnt < 2) {
        return i;
    }
    iavg = tuna_stats_fastavg(&ks[0].stats);
    ivar = tuna_stats_fastvar(&ks[0].stats);
    for (j = 1; j < nk; ++j) {
        jcnt = tuna_stats_cnt(&ks[j].stats);
        if (jcnt < 2) {
            return j;
        }
        javg = tuna_stats_fastavg(&ks[j].stats);
        jvar = tuna_stats_fastvar(&ks[j].stats);
        p    = tuna_welch1(iavg, ivar, icnt, javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int
tuna_algo_zero(const int nk,
               const tuna_chunk* ks,
               tuna_seed* seed)
{
    (void) nk;
    (void) ks;
    (void) seed;
    return 0;
}

// TODO Do something intelligent with clock_getres(2) information

int
tuna_pre_cost(tuna_site* st,
              const tuna_chunk* ks,
              const int nk)
{
    // Ensure a zero-initialize st argument produces good behavior by...
    if (!st->al) {
        // ...providing a default algorithm when not set, and
        st->al = tuna_algo_default(nk);

        if (!st->sd) {
            // ...providing a default seed when not set.
            st->sd = tuna_seed_default();
        }
    }

    // Invoke chosen algorithm recording selection index for tuna_post_cost().
    st->ik = st->al(nk, ks, &st->sd);

    return st->ik;
}

void
tuna_post_cost(tuna_site*  st,
               tuna_chunk* ks,
               const double cost)
{
    tuna_chunk_obs(ks + st->ik, cost);
}

int
tuna_pre(tuna_site* st,
         const tuna_chunk* ks,
         const int nk)
{
    // Delegate chunk selection to tuna_pre_cost
    tuna_pre_cost(st, ks, nk);

    // Glimpse at the clock so we may compute elapsed time in tuna_post()
    clock_gettime(TUNA_CLOCK, &st->ts);

    return st->ik;
}

double
tuna_post(tuna_site*  st,
          tuna_chunk* ks)
{
    // Glimpse at the clock and compute double-valued elapsed time
    struct timespec te;
    clock_gettime(TUNA_CLOCK, &te);
    double elapsed = te.tv_nsec - st->ts.tv_nsec; // Nanoseconds...
    elapsed *= 1e-9;                              // ...to seconds
    elapsed += te.tv_sec - st->ts.tv_sec;         // ...plus seconds

    // Delegate recording the observation to tuna_post_cost()
    tuna_post_cost(st, ks, elapsed);

    return elapsed;
}
