#ifndef TUNA_AS3_H
#define TUNA_AS3_H

/** \file
 * Provides \ref tuna_as3.
 */

#ifdef __cplusplus
extern "C" {
#endif

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
double tuna_as3(double t, int nu);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
