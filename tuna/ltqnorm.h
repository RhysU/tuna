/** \file
 * Provides \ref tuna_ltqnorm.
 */

#ifndef TUNA_LTQNORM_H
#define TUNA_LTQNORM_H

#ifdef __cplusplus
extern "C" {
#endif

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
double tuna_ltqnorm(const double p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_LTQNORM_H */
