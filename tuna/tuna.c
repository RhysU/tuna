/*
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * Tuna implementation.
 */

#include <tuna.h>

// C99 extern declarations for inlined accessors for tuna_stats_*
extern size_t tuna_stats_cnt    (const tuna_stats* const t);
extern double tuna_stats_avg    (const tuna_stats* const t);
extern double tuna_stats_fastavg(const tuna_stats* const t);
extern double tuna_stats_var    (const tuna_stats* const t);
extern double tuna_stats_fastvar(const tuna_stats* const t);
extern double tuna_stats_std    (const tuna_stats* const t);
extern double tuna_stats_faststd(const tuna_stats* const t);

tuna_stats* tuna_stats_fastobs(tuna_stats* const t,
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

tuna_stats* tuna_stats_obs(tuna_stats* const t,
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

tuna_stats* tuna_stats_nobs(tuna_stats* const t,
                            const double* x,
                            size_t N)
{
    if (N) {                               // NOP on degenerate input
        tuna_stats_obs(t, *x++);           // Delegate possible n == 1
        for (size_t i = --N; i -- > 0 ;) { // Henceforth, certainly n > 1
            tuna_stats_fastobs(t, *x++);
        }
    }
    return t;
}

tuna_stats* tuna_stats_merge(tuna_stats* const dst,
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

tuna_kernel* tuna_kernel_obs(tuna_kernel* const k, double t)
{
    // First, find smallest observation among set {t, k->outliers[0], ... }
    // placing it into storage t while maintaining sorted-ness of k->outliers.
    // The loop is one bubble sort pass with possibility of short-circuiting.
    if (enforce_lt(&t, k->outliers)) {
        for (size_t i = 1;
             i < sizeof(k->outliers) / sizeof(k->outliers[0])
             && enforce_lt(k->outliers - 1 + i, k->outliers + i);
             ++i)
            ;
    }

    // Second, when non-zero, record statistics about the best observation.
    if (t) {
        tuna_stats_obs(&k->stats, t);
    }

    // Together, these two steps cause a zero-initialized tuna_kernel to
    // gather tuna_noutliers pieces of information before beginning to
    // track any statistics.  This effectively provides some "start up"
    // or "burn in" period in addition to preventing highly improbable
    // observations from unduly inflating the discovered variability.

    return k;
}

tuna_stats tuna_kernel_stats(tuna_kernel* const k)
{
    tuna_stats s = k->stats;
    return s;
}

tuna_stats* tuna_kernel_merge(tuna_stats*   const s,
                              const tuna_kernel* const k)
{
    tuna_stats_merge(s, &k->stats);
    for (size_t i = 0; i < tuna_countof(k->outliers); ++i) {
        if (k->outliers[i]) {
            return tuna_stats_nobs(s, k->outliers + i,
                                   tuna_countof(k->outliers) - i);
        }
    }
    return s;
}

double tuna_ltqnorm(const double p)
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

// See header for source attribution
double tuna_as3(double t, int nu)
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

// C99 extern declarations for tuna_rand_*
extern double tuna_rand_u01(tuna_seed* sd);
extern double tuna_rand_n01(tuna_seed* sd);

tuna_seed tuna_seed_default()
{
    const char * d = getenv("TUNA_SEED");
    unsigned int retval;
    if (!d || 1 != sscanf(d, " %u", &retval)) {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        retval = ts.tv_sec + ts.tv_nsec;
    }
    return retval;
}

// C99 extern declarations for inlined tuna_welch_*
extern void   tuna_welch(double xA, double sA2, size_t nA,
                         double xB, double sB2, size_t nB,
                         double* const t, double* const nu);
extern double tuna_welch_t (double xA, double sA2, size_t nA,
                            double xB, double sB2, size_t nB);

double tuna_welch1_nuinf(double xA, double sA2, size_t nA,
                         double xB, double sB2, size_t nB)
{
    return 1 - erfc(-tuna_welch_t(xA, sA2, nA, xB, sB2, nB) * M_SQRT1_2) / 2;
}

double tuna_welch1_approx(double xA, double sA2, size_t nA,
                          double xB, double sB2, size_t nB)
{
    double t, nu;
    tuna_welch(xA, sA2, nA, xB, sB2, nB, &t, &nu);
    if (nu > 2) {
        t *= nu / (nu - 2);
    }
    return 1 - erfc(-t * M_SQRT1_2) / 2;
}

double tuna_welch1(double xA, double sA2, size_t nA,
                   double xB, double sB2, size_t nB)
{
    double t, nu;
    tuna_welch(xA, sA2, nA, xB, sB2, nB, &t, &nu);
    return 1 - tuna_as3(t, nu);
}

// http://agentzlerich.blogspot.com/2011/01/c-header-only-unit-testing-with-fctx.html
static inline
void trim(char * const a)
{
    char *p = a, *q = a;
    while (isspace(*q))            ++q;
    while (*q)                     *p++ = *q++;
    *p = '\0';
    while (p > a && isspace(*--p)) *p = '\0';
}

// TODO Is welch1_nuinf so much better for runtime that it should be default?
//      The behavioral difference can be seen in, e.g., ./examples/smallsort.
tuna_algo tuna_algo_default(void)
{
    char * d = getenv("TUNA_ALGO");
    if (d) {
        trim(d);
        if (!strncasecmp(d, "welch1", sizeof("welch1"))) {
            return &tuna_algo_welch1;
        } else if (!strncasecmp(d, "welch1_nuinf", sizeof("welch1_nuinf"))) {
            return &tuna_algo_welch1_nuinf;
        } else if (!strncasecmp(d, "zero", sizeof("zero"))) {
            return &tuna_algo_zero;
        }
    }
    return &tuna_algo_welch1_nuinf; // Default
}

int tuna_algo_welch1_nuinf(const int nk,
                           const tuna_kernel* ks,
                           tuna_seed* seed)
{
    assert(nk > 0);
    int i = 0;
    double icnt = tuna_stats_cnt    (&ks[0].stats);
    if (icnt < 2) return i;
    double iavg = tuna_stats_fastavg(&ks[0].stats);
    double ivar = tuna_stats_fastvar(&ks[0].stats);
    for (int j = 1; j < nk; ++j) {
        const double jcnt = tuna_stats_cnt    (&ks[j].stats);
        if (jcnt < 2) return j;
        const double javg = tuna_stats_fastavg(&ks[j].stats);
        const double jvar = tuna_stats_fastvar(&ks[j].stats);
        const double p    = tuna_welch1_nuinf (iavg, ivar, icnt,
                                               javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int tuna_algo_welch1(const int nk,
                     const tuna_kernel* ks,
                     tuna_seed* seed)
{
    assert(nk > 0);
    int i = 0;
    double icnt = tuna_stats_cnt    (&ks[0].stats);
    if (icnt < 2) return i;
    double iavg = tuna_stats_fastavg(&ks[0].stats);
    double ivar = tuna_stats_fastvar(&ks[0].stats);
    for (int j = 1; j < nk; ++j) {
        const double jcnt = tuna_stats_cnt    (&ks[j].stats);
        if (jcnt < 2) return j;
        const double javg = tuna_stats_fastavg(&ks[j].stats);
        const double jvar = tuna_stats_fastvar(&ks[j].stats);
        const double p    = tuna_welch1       (iavg, ivar, icnt,
                                               javg, jvar, jcnt);
        if (p < tuna_rand_u01(seed)) {
            i    = j;
            iavg = javg;
            ivar = jvar;
            icnt = jcnt;
        }
    }
    return i;
}

int tuna_algo_zero(const int nk,
                   const tuna_kernel* ks,
                   tuna_seed* seed)
{
    (void) nk;
    (void) ks;
    (void) seed;
    return 0;
}

// TODO Do something intelligent with clock_getres(2) information

int tuna_pre(tuna_site* st,
             const tuna_kernel* ks,
             const int nk)
{
    // Ensure a zero-initialize st argument produces good behavior by...
    if (!st->al) {
        // ...providing a default algorithm when not set, and
        st->al = tuna_algo_default();

        if (!st->sd) {
            // ...providing a default seed when not set.
            st->sd = tuna_seed_default();
        }
    }

    // Invoke chosen algorithm recording selection index for tuna_post_cost().
    st->ik = st->al(nk, ks, &st->sd);

    // Glimpse at the clock so we may compute elapsed time in tuna_post()
    clock_gettime(TUNA_CLOCK, &st->ts);

    return st->ik;
}

inline
void tuna_post_cost(tuna_site*  st,
                    tuna_kernel* ks,
                    const double cost)
{
    tuna_kernel_obs(ks + st->ik, cost);
}

double tuna_post(tuna_site*  st,
                 tuna_kernel* ks)
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
