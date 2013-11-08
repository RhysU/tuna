// See source code attribution in ltqnorm.h

/** \file
 * \copydoc ltqnorm.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <float.h>
#include <math.h>

#include "ltqnorm.h"

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
