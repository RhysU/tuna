/*
 * Originally from CDFLIB: http://www.netlib.org/random/dcdflib.c.tar.gz
 * Obtained from http://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html
 *
 * CDFLIB was produced by
 *    Barry Brown, James Lovato, Kathy Russell,
 *    Department of Biomathematics,
 *    University of Texas,
 *    Houston, Texas.
 */

#ifndef TUNA_CDFLIB_H
#define TUNA_CDFLIB_H

/** \file
 * A modified version of CDFLIB by
 * Barry Brown, James Lovato, and Kathy Russell.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief
 *   ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
 *
 * <pre>
 * Discussion:
 *
 *   In this algorithm, DEL(X) is the function defined by
 *
 *     ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
 *                     + DEL(X).
 *
 * Parameters:
 *
 *   Input, double *A, *B, define the arguments.
 *
 *   Output, double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
 * </pre>
 */
double
cdflib_algdiv(
    const double* a,
    const double* b);

/**
 * \brief
 *   ALNREL evaluates the function ln ( 1 + A ).
 *
 * <pre>
 * Modified:
 *
 *   17 November 2006
 *
 * Reference:
 *
 *   Armido DiDinato, Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios,
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, 1993, pages 360-373.
 *
 * Parameters:
 *
 *   Input, double *A, the argument.
 *
 *   Output, double ALNREL, the value of ln ( 1 + A ).
 * </pre>
 */
double
cdflib_alnrel(
    const double* a);

/**
 * \brief
 *   APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
 *
 * <pre>
 * Discussion:
 *
 *   APSER is used only for cases where
 *
 *     A <= min ( EPS, EPS * B ),
 *     B * X <= 1, and
 *     X <= 0.5.
 *
 * Parameters:
 *
 *   Input, double *A, *B, *X, the parameters of teh
 *   incomplete beta ratio.
 *
 *   Input, double *EPS, a tolerance.
 *
 *   Output, double APSER, the computed value of the
 *   incomplete beta ratio.
 * </pre>
 */
double
cdflib_apser(
    const double* a,
    const double* b,
    const double* x,
    const double* eps);


/**
 * \brief
 *   BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
 *
 * <pre>
 * Discussion:
 *
 *   The function DEL(A) is a remainder term that is used in the expression:
 *
 *     ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
 *       - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
 *
 *   or, in other words, DEL ( A ) is defined as:
 *
 *     DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
 *       + A + 0.5 * ln ( 2 * PI ).
 *
 * Parameters:
 *
 *   Input, double *A0, *B0, the arguments.
 *   It is assumed that 8 <= A0 and 8 <= B0.
 *
 *   Output, double *BCORR, the value of the function.
 * </pre>
 */
double
cdflib_bcorr(
    const double* a0,
    const double* b0);

/**
 * \brief
 *   BETA evaluates the beta function.
 *
 * <pre>
 * Modified:
 *
 *   03 December 1999
 *
 * Author:
 *
 *   John Burkardt
 *
 * Parameters:
 *
 *   Input, double A, B, the arguments of the beta function.
 *
 *   Output, double BETA, the value of the beta function.
 * </pre>
 */
double
cdflib_beta(
    double a,
    double b);

/**
 * \brief
 *   BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the function.
 *   A and B should be nonnegative.  It is assumed that both A and B
 *   are greater than or equal to 15.
 *
 *   Input, double *LAMBDA, the value of ( A + B ) * Y - B.
 *   It is assumed that 0 <= LAMBDA.
 *
 *   Input, double *EPS, the tolerance.
 * </pre>
 */
double
cdflib_beta_asym(
    const double* a,
    const double* b,
    const double* lambda,
    const double* eps);

/**
 * \brief
 *   BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the function.
 *   A and B should be nonnegative.  It is assumed that both A and
 *   B are greater than 1.
 *
 *   Input, double *X, *Y.  X is the argument of the
 *   function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
 *
 *   Input, double *LAMBDA, the value of ( A + B ) * Y - B.
 *
 *   Input, double *EPS, a tolerance.
 *
 *   Output, double BETA_FRAC, the value of the continued
 *   fraction approximation for IX(A,B).
 * </pre>
 */
double
cdflib_beta_frac(
    const double* a,
    const double* b,
    const double* x,
    const double* y,
    const double* lambda,
    const double* eps);

/**
 * \brief
 *   BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the function.
 *   A and B should be nonnegative.  It is assumed that 15 <= A
 *   and B <= 1, and that B is less than A.
 *
 *   Input, double *X, *Y.  X is the argument of the
 *   function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
 *
 *   Input/output, double *W, a quantity to which the
 *   result of the computation is to be added on output.
 *
 *   Input, double *EPS, a tolerance.
 *
 *   Output, int *IERR, an error flag, which is 0 if no error
 *   was detected.
 * </pre>
 */
void
cdflib_beta_grat(
    const double* a,
    const double* b,
    const double* x,
    const double* y,
    double* w,
    const double* eps,
    int* ierr);

/**
 * \brief
 *   BETA_INC evaluates the incomplete beta function IX(A,B).
 *
 * <pre>
 * Author:
 *
 *   Alfred H Morris, Jr,
 *   Naval Surface Weapons Center,
 *   Dahlgren, Virginia.
 *
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the function.
 *   A and B should be nonnegative.
 *
 *   Input, double *X, *Y.  X is the argument of the
 *   function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
 *
 *   Output, double *W, *W1, the values of IX(A,B) and
 *   1-IX(A,B).
 *
 *   Output, int *IERR, the error flag.
 *   0, no error was detected.
 *   1, A or B is negative;
 *   2, A = B = 0;
 *   3, X < 0 or 1 < X;
 *   4, Y < 0 or 1 < Y;
 *   5, X + Y /= 1;
 *   6, X = A = 0;
 *   7, Y = B = 0.
 * </pre>
 */
void
cdflib_beta_inc(
    const double* a,
    const double* b,
    const double* x,
    const double* y,
    double* w,
    double* w1,
    int* ierr);

/**
 * \brief
 *   BETA_INC_VALUES returns some values of the incomplete Beta function.
 *
 * <pre>
 * Discussion:
 *
 *   The incomplete Beta function may be written
 *
 *     BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
 *                     / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
 *
 *   Thus,
 *
 *     BETA_INC(A,B,0.0) = 0.0
 *     BETA_INC(A,B,1.0) = 1.0
 *
 *   Note that in Mathematica, the expressions:
 *
 *     BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
 *     BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
 *
 *   and thus, to evaluate the incomplete Beta function requires:
 *
 *     BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
 *
 * Modified:
 *
 *   09 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Karl Pearson,
 *   Tables of the Incomplete Beta Function,
 *   Cambridge University Press, 1968.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *A, *B, the parameters of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_beta_inc_values(
    int* n_data,
    double* a,
    double* b,
    double* x,
    double* fx);

/**
 * \brief
 *   BETA_LOG evaluates the logarithm of the beta function.
 *
 * <pre>
 * Reference:
 *
 *   Armido DiDinato and Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios,
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, 1993, pages 360-373.
 *
 * Parameters:
 *
 *   Input, double *A0, *B0, the parameters of the function.
 *   A0 and B0 should be nonnegative.
 *
 *   Output, double *BETA_LOG, the value of the logarithm
 *   of the Beta function.
 * </pre>
 */
double
cdflib_beta_log(
    const double* a0,
    const double* b0);

/**
 * \brief
 *   BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
 *
 * <pre>
 * Discussion:
 *
 *   BETA_PSER is used when B <= 1 or B*X <= 0.7.
 *
 * Parameters:
 *
 *   Input, double *A, *B, the parameters.
 *
 *   Input, double *X, the point where the function
 *   is to be evaluated.
 *
 *   Input, double *EPS, the tolerance.
 *
 *   Output, double BETA_PSER, the approximate value of IX(A,B)(X).
 * </pre>
 */
double
cdflib_beta_pser(
    const double* a,
    const double* b,
    const double* x,
    const double* eps);

/**
 * \brief
 *   BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the Beta function.
 *   A and B should be nonnegative.
 *
 *   Input, double *X, *Y, define the numerator of the fraction.
 *
 *   Output, double BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
 * </pre>
 */
double
cdflib_beta_cdflib_rcomp(
    const double* a,
    const double* b,
    const double* x,
    const double* y);

/**
 * \brief
 *   BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
 *
 * <pre>
 * Parameters:
 *
 *   Input, int MU, ?
 *
 *   Input, double A, B, the parameters of the Beta function.
 *   A and B should be nonnegative.
 *
 *   Input, double X, Y, ?
 *
 *   Output, double BETA_RCOMP1, the value of
 *   exp(MU) * X**A * Y**B / Beta(A,B).
 * </pre>
 */
double
cdflib_beta_rcomp1(
    const int* mu,
    const double* a,
    const double* b,
    const double* x,
    const double* y);

/**
 * \brief
 *   BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the function.
 *   A and B should be nonnegative.
 *
 *   Input, double *X, *Y, ?
 *
 *   Input, int *N, ?
 *
 *   Input, double *EPS, the tolerance.
 *
 *   Output, double BETA_UP, the value of IX(A,B) - IX(A+N,B).
 * </pre>
 */
double
cdflib_beta_up(
    const double* a,
    const double* b,
    const double* x,
    const double* y,
    const int* n,
    const double* eps);

/**
 * \brief
 *   BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
 *
 * <pre>
 * Discussion:
 *
 *   CDF(X)(A,B) is the probability of at most X successes in A trials,
 *   given that the probability of success on a single trial is B.
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Daniel Zwillinger,
 *   CRC Standard Mathematical Tables and Formulae,
 *   30th Edition, CRC Press, 1996, pages 651-652.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *A, double *B, the parameters of the function.
 *
 *   Output, int *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_binomial_cdf_values(
    int* n_data,
    int* a,
    double* b,
    int* x,
    double* fx);

/**
 * \brief
 *   CDFBET evaluates the CDF of the Beta Distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the beta distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly by code associated with the reference.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The beta density is proportional to t^(A-1) * (1-t)^(B-1).
 *
 * Modified:
 *
 *   09 June 2004
 *
 * Reference:
 *
 *   Armido DiDinato and Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios,
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, 1993, pages 360-373.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which of the next four argument
 *   values is to be calculated from the others.
 *   1: Calculate P and Q from X, Y, A and B;
 *   2: Calculate X and Y from P, Q, A and B;
 *   3: Calculate A from P, Q, X, Y and B;
 *   4: Calculate B from P, Q, X, Y and A.
 *
 *   Input/output, double *P, the integral from 0 to X of the
 *   chi-square distribution.  Input range: [0, 1].
 *
 *   Input/output, double *Q, equals 1-P.  Input range: [0, 1].
 *
 *   Input/output, double *X, the upper limit of integration
 *   of the beta density.  If it is an input value, it should lie in
 *   the range [0,1].  If it is an output value, it will be searched for
 *   in the range [0,1].
 *
 *   Input/output, double *Y, equal to 1-X.  If it is an input
 *   value, it should lie in the range [0,1].  If it is an output value,
 *   it will be searched for in the range [0,1].
 *
 *   Input/output, double *A, the first parameter of the beta
 *   density.  If it is an input value, it should lie in the range
 *   (0, +infinity).  If it is an output value, it will be searched
 *   for in the range [1D-300,1D300].
 *
 *   Input/output, double *B, the second parameter of the beta
 *   density.  If it is an input value, it should lie in the range
 *   (0, +infinity).  If it is an output value, it will be searched
 *   for in the range [1D-300,1D300].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1;
 *   +4, if X + Y /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfbet(
    const int* which,
    double* p,
    double* q,
    double* x,
    double* y,
    double* a,
    double* b,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFBIN evaluates the CDF of the Binomial distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the binomial distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   P is the probability of S or fewer successes in XN binomial trials,
 *   each trial having an individual probability of success of PR.
 *
 * Modified:
 *
 *   09 June 2004
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.24.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which of argument values is to
 *   be calculated from the others.
 *   1: Calculate P and Q from S, XN, PR and OMPR;
 *   2: Calculate S from P, Q, XN, PR and OMPR;
 *   3: Calculate XN from P, Q, S, PR and OMPR;
 *   4: Calculate PR and OMPR from P, Q, S and XN.
 *
 *   Input/output, double *P, the cumulation, from 0 to S,
 *   of the binomial distribution.  If P is an input value, it should
 *   lie in the range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *S, the number of successes observed.
 *   Whether this is an input or output value, it should lie in the
 *   range [0,XN].
 *
 *   Input/output, double *XN, the number of binomial trials.
 *   If this is an input value it should lie in the range: (0, +infinity).
 *   If it is an output value it will be searched for in the
 *   range [1.0D-300, 1.0D+300].
 *
 *   Input/output, double *PR, the probability of success in each
 *   binomial trial.  Whether this is an input or output value, it should
 *   lie in the range: [0,1].
 *
 *   Input/output, double *OMPR, equal to 1-PR.  Whether this is an
 *   input or output value, it should lie in the range [0,1].  Also, it should
 *   be the case that PR + OMPR = 1.
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1;
 *   +4, if PR + OMPR /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfbin(
    const int* which,
    double* p,
    double* q,
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFCHI evaluates the CDF of the chi square distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the chi square distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The CDF of the chi square distribution can be evaluated
 *   within Mathematica by commands such as:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF [ ChiSquareDistribution [ DF ], X ]
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.4.19.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from X and DF;
 *   2: Calculate X from P, Q and DF;
 *   3: Calculate DF from P, Q and X.
 *
 *   Input/output, double *P, the integral from 0 to X of
 *   the chi-square distribution.  If this is an input value, it should
 *   lie in the range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *X, the upper limit of integration
 *   of the chi-square distribution.  If this is an input
 *   value, it should lie in the range: [0, +infinity).  If it is an output
 *   value, it will be searched for in the range: [0,1.0D+300].
 *
 *   Input/output, double *DF, the degrees of freedom of the
 *   chi-square distribution.  If this is an input value, it should lie
 *   in the range: (0, +infinity).  If it is an output value, it will be
 *   searched for in the range: [ 1.0D-300, 1.0D+300].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1;
 *   +10, an error was returned from CUMGAM.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfchi(
    const int* which,
    double* p,
    double* q,
    double* x,
    double* df,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFCHN evaluates the CDF of the Noncentral Chi-Square.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the noncentral chi-square
 *   distribution given values for the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The computation time required for this routine is proportional
 *   to the noncentrality parameter (PNONC).  Very large values of
 *   this parameter can consume immense computer resources.  This is
 *   why the search range is bounded by 10,000.
 *
 *   The CDF of the noncentral chi square distribution can be evaluated
 *   within Mathematica by commands such as:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.25.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from X, DF and PNONC;
 *   2: Calculate X from P, DF and PNONC;
 *   3: Calculate DF from P, X and PNONC;
 *   4: Calculate PNONC from P, X and DF.
 *
 *   Input/output, double *P, the integral from 0 to X of
 *   the noncentral chi-square distribution.  If this is an input
 *   value, it should lie in the range: [0, 1.0-1.0D-16).
 *
 *   Input/output, double *Q, is generally not used by this
 *   subroutine and is only included for similarity with other routines.
 *   However, if P is to be computed, then a value will also be computed
 *   for Q.
 *
 *   Input, double *X, the upper limit of integration of the
 *   noncentral chi-square distribution.  If this is an input value, it
 *   should lie in the range: [0, +infinity).  If it is an output value,
 *   it will be sought in the range: [0,1.0D+300].
 *
 *   Input/output, double *DF, the number of degrees of freedom
 *   of the noncentral chi-square distribution.  If this is an input value,
 *   it should lie in the range: (0, +infinity).  If it is an output value,
 *   it will be searched for in the range: [ 1.0D-300, 1.0D+300].
 *
 *   Input/output, double *PNONC, the noncentrality parameter of
 *   the noncentral chi-square distribution.  If this is an input value, it
 *   should lie in the range: [0, +infinity).  If it is an output value,
 *   it will be searched for in the range: [0,1.0D+4]
 *
 *   Output, int *STATUS, reports on the calculation.
 *   0, if calculation completed correctly;
 *   -I, if input parameter number I is out of range;
 *   1, if the answer appears to be lower than the lowest search bound;
 *   2, if the answer appears to be higher than the greatest search bound.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfchn(
    const int* which,
    double* p,
    double* q,
    double* x,
    double* df,
    double* pnonc,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFF evaluates the CDF of the F distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the F distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The value of the cumulative F distribution is not necessarily
 *   monotone in either degree of freedom.  There thus may be two
 *   values that provide a given CDF value.  This routine assumes
 *   monotonicity and will find an arbitrary one of the two values.
 *
 * Modified:
 *
 *   14 April 2007
 *
 * Reference:
 *
 *   Milton Abramowitz, Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.6.2.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from F, DFN and DFD;
 *   2: Calculate F from P, Q, DFN and DFD;
 *   3: Calculate DFN from P, Q, F and DFD;
 *   4: Calculate DFD from P, Q, F and DFN.
 *
 *   Input/output, double *P, the integral from 0 to F of
 *   the F-density.  If it is an input value, it should lie in the
 *   range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *F, the upper limit of integration
 *   of the F-density.  If this is an input value, it should lie in the
 *   range [0, +infinity).  If it is an output value, it will be searched
 *   for in the range [0,1.0D+300].
 *
 *   Input/output, double *DFN, the number of degrees of
 *   freedom of the numerator sum of squares.  If this is an input value,
 *   it should lie in the range: (0, +infinity).  If it is an output value,
 *   it will be searched for in the range: [ 1.0D-300, 1.0D+300].
 *
 *   Input/output, double *DFD, the number of degrees of freedom
 *   of the denominator sum of squares.  If this is an input value, it should
 *   lie in the range: (0, +infinity).  If it is an output value, it will
 *   be searched for in the  range: [ 1.0D-300, 1.0D+300].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdff(
    const int* which,
    double* p,
    double* q,
    double* f,
    double* dfn,
    double* dfd,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFFNC evaluates the CDF of the Noncentral F distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine originally used 1.0E+300 as the upper bound for the
 *   interval in which many of the missing parameters are to be sought.
 *   Since the underlying rootfinder routine needs to evaluate the
 *   function at this point, it is no surprise that the program was
 *   experiencing overflows.  A less extravagant upper bound
 *   is being tried for now!
 *
 *
 *   This routine calculates any one parameter of the Noncentral F distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The computation time required for this routine is proportional
 *   to the noncentrality parameter PNONC.  Very large values of
 *   this parameter can consume immense computer resources.  This is
 *   why the search range is bounded by 10,000.
 *
 *   The value of the cumulative noncentral F distribution is not
 *   necessarily monotone in either degree of freedom.  There thus
 *   may be two values that provide a given CDF value.  This routine
 *   assumes monotonicity and will find an arbitrary one of the two
 *   values.
 *
 *   The CDF of the noncentral F distribution can be evaluated
 *   within Mathematica by commands such as:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
 *
 * Modified:
 *
 *   15 June 2004
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.6.20.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from F, DFN, DFD and PNONC;
 *   2: Calculate F from P, Q, DFN, DFD and PNONC;
 *   3: Calculate DFN from P, Q, F, DFD and PNONC;
 *   4: Calculate DFD from P, Q, F, DFN and PNONC;
 *   5: Calculate PNONC from P, Q, F, DFN and DFD.
 *
 *   Input/output, double *P, the integral from 0 to F of
 *   the noncentral F-density.  If P is an input value it should
 *   lie in the range [0,1) (Not including 1!).
 *
 *   Dummy, double *Q, is not used by this subroutine,
 *   and is only included for similarity with the other routines.
 *   Its input value is not checked.  If P is to be computed, the
 *   Q is set to 1 - P.
 *
 *   Input/output, double *F, the upper limit of integration
 *   of the noncentral F-density.  If this is an input value, it should
 *   lie in the range: [0, +infinity).  If it is an output value, it
 *   will be searched for in the range: [0,1.0D+30].
 *
 *   Input/output, double *DFN, the number of degrees of freedom
 *   of the numerator sum of squares.  If this is an input value, it should
 *   lie in the range: (0, +infinity).  If it is an output value, it will
 *   be searched for in the range: [ 1.0, 1.0D+30].
 *
 *   Input/output, double *DFD, the number of degrees of freedom
 *   of the denominator sum of squares.  If this is an input value, it should
 *   be in range: (0, +infinity).  If it is an output value, it will be
 *   searched for in the range [1.0, 1.0D+30].
 *
 *   Input/output, double *PNONC, the noncentrality parameter
 *   If this is an input value, it should be nonnegative.
 *   If it is an output value, it will be searched for in the range: [0,1.0D+4].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdffnc(
    const int* which,
    double* p,
    double* q,
    double* f,
    double* dfn,
    double* dfd,
    double* phonc,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFGAM evaluates the CDF of the Gamma Distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the Gamma distribution
 *   given the others.
 *
 *   The cumulative distribution function P is calculated directly.
 *
 *   Computation of the other parameters involves a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The gamma density is proportional to T**(SHAPE - 1) * EXP(- SCALE * T)
 *
 * Reference:
 *
 *   Armido DiDinato and Alfred Morris,
 *   Computation of the incomplete gamma function ratios and their inverse,
 *   ACM Transactions on Mathematical Software,
 *   Volume 12, 1986, pages 377-393.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from X, SHAPE and SCALE;
 *   2: Calculate X from P, Q, SHAPE and SCALE;
 *   3: Calculate SHAPE from P, Q, X and SCALE;
 *   4: Calculate SCALE from P, Q, X and SHAPE.
 *
 *   Input/output, double *P, the integral from 0 to X of the
 *   Gamma density.  If this is an input value, it should lie in the
 *   range: [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *X, the upper limit of integration of
 *   the Gamma density.  If this is an input value, it should lie in the
 *   range: [0, +infinity).  If it is an output value, it will lie in
 *   the range: [0,1E300].
 *
 *   Input/output, double *SHAPE, the shape parameter of the
 *   Gamma density.  If this is an input value, it should lie in the range:
 *   (0, +infinity).  If it is an output value, it will be searched for
 *   in the range: [1.0D-300,1.0D+300].
 *
 *   Input/output, double *SCALE, the scale parameter of the
 *   Gamma density.  If this is an input value, it should lie in the range
 *   (0, +infinity).  If it is an output value, it will be searched for
 *   in the range: (1.0D-300,1.0D+300].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1;
 *   +10, if the Gamma or inverse Gamma routine cannot compute the answer.
 *   This usually happens only for X and SHAPE very large (more than 1.0D+10.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfgam(
    const int* which,
    double* p,
    double* q,
    double* x,
    double* shape,
    double* scale,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFNBN evaluates the CDF of the Negative Binomial distribution
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the negative binomial
 *   distribution given values for the others.
 *
 *   The cumulative negative binomial distribution returns the
 *   probability that there will be F or fewer failures before the
 *   S-th success in binomial trials each of which has probability of
 *   success PR.
 *
 *   The individual term of the negative binomial is the probability of
 *   F failures before S successes and is
 *   Choose( F, S+F-1 ) * PR^(S) * (1-PR)^F
 *
 *   Computation of other parameters involve a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.26.
 *
 * Parameters:
 *
 *   Input, int WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from F, S, PR and OMPR;
 *   2: Calculate F from P, Q, S, PR and OMPR;
 *   3: Calculate S from P, Q, F, PR and OMPR;
 *   4: Calculate PR and OMPR from P, Q, F and S.
 *
 *   Input/output, double P, the cumulation from 0 to F of
 *   the negative binomial distribution.  If P is an input value, it
 *   should lie in the range [0,1].
 *
 *   Input/output, double Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double F, the upper limit of cumulation of
 *   the binomial distribution.  There are F or fewer failures before
 *   the S-th success.  If this is an input value, it may lie in the
 *   range [0,+infinity), and if it is an output value, it will be searched
 *   for in the range [0,1.0D+300].
 *
 *   Input/output, double S, the number of successes.
 *   If this is an input value, it should lie in the range: [0, +infinity).
 *   If it is an output value, it will be searched for in the range:
 *   [0, 1.0D+300].
 *
 *   Input/output, double PR, the probability of success in each
 *   binomial trial.  Whether an input or output value, it should lie in the
 *   range [0,1].
 *
 *   Input/output, double OMPR, the value of (1-PR).  Whether an
 *   input or output value, it should lie in the range [0,1].
 *
 *   Output, int STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1;
 *   +4, if PR + OMPR /= 1.
 *
 *   Output, double BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfnbn(
    const int* which,
    double* p,
    double* q,
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFNOR evaluates the CDF of the Normal distribution.
 *
 * <pre>
 * Discussion:
 *
 *   A slightly modified version of ANORM from SPECFUN
 *   is used to calculate the cumulative standard normal distribution.
 *
 *   The rational functions from pages 90-95 of Kennedy and Gentle
 *   are used as starting values to Newton's Iterations which
 *   compute the inverse standard normal.  Therefore no searches are
 *   necessary for any parameter.
 *
 *   For X < -15, the asymptotic expansion for the normal is used  as
 *   the starting value in finding the inverse standard normal.
 *
 *   The normal density is proportional to
 *   exp( - 0.5D+00 * (( X - MEAN)/SD)**2)
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.2.12.
 *
 *   William Cody,
 *   Algorithm 715: SPECFUN - A Portable FORTRAN Package of
 *     Special Function Routines and Test Drivers,
 *   ACM Transactions on Mathematical Software,
 *   Volume 19, pages 22-32, 1993.
 *
 *   Kennedy and Gentle,
 *   Statistical Computing,
 *   Marcel Dekker, NY, 1980,
 *   QA276.4  K46
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from X, MEAN and SD;
 *   2: Calculate X from P, Q, MEAN and SD;
 *   3: Calculate MEAN from P, Q, X and SD;
 *   4: Calculate SD from P, Q, X and MEAN.
 *
 *   Input/output, double *P, the integral from -infinity to X
 *   of the Normal density.  If this is an input or output value, it will
 *   lie in the range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *X, the upper limit of integration of
 *   the Normal density.
 *
 *   Input/output, double *MEAN, the mean of the Normal density.
 *
 *   Input/output, double *SD, the standard deviation of the
 *   Normal density.  If this is an input value, it should lie in the
 *   range (0,+infinity).
 *
 *   Output, int *STATUS, the status of the calculation.
 *   0, if calculation completed correctly;
 *   -I, if input parameter number I is out of range;
 *   1, if answer appears to be lower than lowest search bound;
 *   2, if answer appears to be higher than greatest search bound;
 *   3, if P + Q /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfnor(
    const int* which,
    double* p,
    double* q,
    double* x,
    double* mean,
    double* sd,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFPOI evaluates the CDF of the Poisson distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the Poisson distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of other parameters involve a search for a value that
 *   produces the desired value of P.  The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.4.21.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1: Calculate P and Q from S and XLAM;
 *   2: Calculate A from P, Q and XLAM;
 *   3: Calculate XLAM from P, Q and S.
 *
 *   Input/output, double *P, the cumulation from 0 to S of the
 *   Poisson density.  Whether this is an input or output value, it will
 *   lie in the range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *S, the upper limit of cumulation of
 *   the Poisson CDF.  If this is an input value, it should lie in
 *   the range: [0, +infinity).  If it is an output value, it will be
 *   searched for in the range: [0,1.0D+300].
 *
 *   Input/output, double *XLAM, the mean of the Poisson
 *   distribution.  If this is an input value, it should lie in the range
 *   [0, +infinity).  If it is an output value, it will be searched for
 *   in the range: [0,1E300].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdfpoi(
    const int* which,
    double* p,
    double* q,
    double* s,
    double* xlam,
    int* status,
    double* bound);

/**
 * \brief
 *   CDFT evaluates the CDF of the T distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates any one parameter of the T distribution
 *   given the others.
 *
 *   The value P of the cumulative distribution function is calculated
 *   directly.
 *
 *   Computation of other parameters involve a search for a value that
 *   produces the desired value of P.   The search relies on the
 *   monotonicity of P with respect to the other parameters.
 *
 *   The original version of this routine allowed the search interval
 *   to extend from -1.0E+300 to +1.0E+300, which is fine until you
 *   try to evaluate a function at such a point!
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.27.
 *
 * Parameters:
 *
 *   Input, int *WHICH, indicates which argument is to be calculated
 *   from the others.
 *   1 : Calculate P and Q from T and DF;
 *   2 : Calculate T from P, Q and DF;
 *   3 : Calculate DF from P, Q and T.
 *
 *   Input/output, double *P, the integral from -infinity to T of
 *   the T-density.  Whether an input or output value, this will lie in the
 *   range [0,1].
 *
 *   Input/output, double *Q, equal to 1-P.  If Q is an input
 *   value, it should lie in the range [0,1].  If Q is an output value,
 *   it will lie in the range [0,1].
 *
 *   Input/output, double *T, the upper limit of integration of
 *   the T-density.  If this is an input value, it may have any value.
 *   It it is an output value, it will be searched for in the range
 *   [ -1.0D+30, 1.0D+30 ].
 *
 *   Input/output, double *DF, the number of degrees of freedom
 *   of the T distribution.  If this is an input value, it should lie
 *   in the range: (0 , +infinity).  If it is an output value, it will be
 *   searched for in the range: [1, 1.0D+10].
 *
 *   Output, int *STATUS, reports the status of the computation.
 *    0, if the calculation completed correctly;
 *   -I, if the input parameter number I is out of range;
 *   +1, if the answer appears to be lower than lowest search bound;
 *   +2, if the answer appears to be higher than greatest search bound;
 *   +3, if P + Q /= 1.
 *
 *   Output, double *BOUND, is only defined if STATUS is nonzero.
 *   If STATUS is negative, then this is the value exceeded by parameter I.
 *   if STATUS is 1 or 2, this is the search bound that was exceeded.
 * </pre>
 */
void
cdflib_cdft(
    const int* which,
    double* p,
    double* q,
    double* t,
    double* df,
    int* status,
    double* bound);

/**
 * \brief
 *   CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
 *
 * <pre>
 * Discussion:
 *
 *   The CDF of the noncentral chi square distribution can be evaluated
 *   within Mathematica by commands such as:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
 *
 * Modified:
 *
 *   12 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *LAMBDA, the noncentrality parameter.
 *
 *   Output, int *DF, the number of degrees of freedom.
 *
 *   Output, double *CDF, the noncentral chi CDF.
 * </pre>
 */
void
cdflib_chi_noncentral_cdf_values(
    int* n_data,
    double* x,
    double* lambda,
    int* df,
    double* cdf);

/**
 * \brief
 *   CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
 *
 * <pre>
 * Discussion:
 *
 *   The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
 *   commands like:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF[ChiSquareDistribution[DF], X ]
 *
 * Modified:
 *
 *   11 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *A, the parameter of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_chi_square_cdf_values(
    int* n_data,
    int* a,
    double* x,
    double* fx);

/**
 * \brief
 *   CUMBET evaluates the cumulative incomplete beta distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine calculates the CDF to X of the incomplete beta distribution
 *   with parameters A and B.  This is the integral from 0 to x
 *   of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
 *
 * Modified:
 *
 *   14 March 2006
 *
 * Reference:
 *
 *   A R Didonato and Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios.
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, Number 3, September 1992, pages 360-373.
 *
 * Parameters:
 *
 *   Input, double *X, the upper limit of integration.
 *
 *   Input, double *Y, the value of 1-X.
 *
 *   Input, double *A, *B, the parameters of the distribution.
 *
 *   Output, double *CUM, *CCUM, the values of the cumulative
 *   density function and complementary cumulative density function.
 * </pre>
 */
void
cdflib_cumbet(
    const double* x,
    const double* y,
    const double* a,
    const double* b,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMBIN evaluates the cumulative binomial distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine returns the probability of 0 to S successes in XN binomial
 *   trials, each of which has a probability of success, PR.
 *
 * Modified:
 *
 *   14 March 2006
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.24.
 *
 * Parameters:
 *
 *   Input, double *S, the upper limit of summation.
 *
 *   Input, double *XN, the number of trials.
 *
 *   Input, double *PR, the probability of success in one trial.
 *
 *   Input, double *OMPR, equals ( 1 - PR ).
 *
 *   Output, double *CUM, the cumulative binomial distribution.
 *
 *   Output, double *CCUM, the complement of the cumulative
 *   binomial distribution.
 * </pre>
 */
void
cdflib_cumbin(
    const double* s,
    const double* xn,
    const double* pr,
    const double* ompr,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMCHI evaluates the cumulative chi-square distribution.
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *X, the upper limit of integration.
 *
 *   Input, double *DF, the degrees of freedom of the
 *   chi-square distribution.
 *
 *   Output, double *CUM, the cumulative chi-square distribution.
 *
 *   Output, double *CCUM, the complement of the cumulative
 *   chi-square distribution.
 * </pre>
 */
void
cdflib_cumchi(
    const double* x,
    const double* df,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMCHN evaluates the cumulative noncentral chi-square distribution.
 *
 * <pre>
 * Discussion:
 *
 *   Calculates the cumulative noncentral chi-square
 *   distribution, i.e., the probability that a random variable
 *   which follows the noncentral chi-square distribution, with
 *   noncentrality parameter PNONC and continuous degrees of
 *   freedom DF, is less than or equal to X.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.4.25.
 *
 * Parameters:
 *
 *   Input, double *X, the upper limit of integration.
 *
 *   Input, double *DF, the number of degrees of freedom.
 *
 *   Input, double *PNONC, the noncentrality parameter of
 *   the noncentral chi-square distribution.
 *
 *   Output, double *CUM, *CCUM, the CDF and complementary
 *   CDF of the noncentral chi-square distribution.
 *
 * Local Parameters:
 *
 *   Local, double EPS, the convergence criterion.  The sum
 *   stops when a term is less than EPS*SUM.
 *
 *   Local, int NTIRED, the maximum number of terms to be evaluated
 *   in each sum.
 *
 *   Local, bool QCONV, is TRUE if convergence was achieved, that is,
 *   the program did not stop on NTIRED criterion.
 * </pre>
 */
void
cdflib_cumchn(
    const double* x,
    const double* df,
    const double* pnonc,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMF evaluates the cumulative F distribution.
 *
 * <pre>
 * Discussion:
 *
 *   CUMF computes the integral from 0 to F of the F density with DFN
 *   numerator and DFD denominator degrees of freedom.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.28.
 *
 * Parameters:
 *
 *   Input, double *F, the upper limit of integration.
 *
 *   Input, double *DFN, *DFD, the number of degrees of
 *   freedom for the numerator and denominator.
 *
 *   Output, double *CUM, *CCUM, the value of the F CDF and
 *   the complementary F CDF.
 * </pre>
 */
void
cdflib_cumf(
    const double* f,
    const double* dfn,
    const double* dfd,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMFNC evaluates the cumulative noncentral F distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine computes the noncentral F distribution with DFN and DFD
 *   degrees of freedom and noncentrality parameter PNONC.
 *
 *   The series is calculated backward and forward from J = LAMBDA/2
 *   (this is the term with the largest Poisson weight) until
 *   the convergence criterion is met.
 *
 *   The sum continues until a succeeding term is less than EPS
 *   times the sum (or the sum is less than 1.0e-20).  EPS is
 *   set to 1.0e-4 in a data statement which can be changed.
 *
 *
 *   The original version of this routine allowed the input values
 *   of DFN and DFD to be negative (nonsensical) or zero (which
 *   caused numerical overflow.)  I have forced both these values
 *   to be at least 1.
 *
 * Modified:
 *
 *   15 June 2004
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
 *
 * Parameters:
 *
 *   Input, double *F, the upper limit of integration.
 *
 *   Input, double *DFN, *DFD, the number of degrees of freedom
 *   in the numerator and denominator.  Both DFN and DFD must be positive,
 *   and normally would be integers.  This routine requires that they
 *   be no less than 1.
 *
 *   Input, double *PNONC, the noncentrality parameter.
 *
 *   Output, double *CUM, *CCUM, the noncentral F CDF and
 *   complementary CDF.
 * </pre>
 */
void
cdflib_cumfnc(
    const double* f,
    const double* dfn,
    const double* dfd,
    const double* pnonc,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMGAM evaluates the cumulative incomplete gamma distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine computes the cumulative distribution function of the
 *   incomplete gamma distribution, i.e., the integral from 0 to X of
 *
 *     (1/GAM(A))*EXP(-T)*T**(A-1) DT
 *
 *   where GAM(A) is the complete gamma function of A, i.e.,
 *
 *     GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
 *
 * Parameters:
 *
 *   Input, double *X, the upper limit of integration.
 *
 *   Input, double *A, the shape parameter of the incomplete
 *   Gamma distribution.
 *
 *   Output, double *CUM, *CCUM, the incomplete Gamma CDF and
 *   complementary CDF.
 * </pre>
 */
void
cdflib_cumgam(
    const double* x,
    const double* a,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMNBN evaluates the cumulative negative binomial distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This routine returns the probability that there will be F or
 *   fewer failures before there are S successes, with each binomial
 *   trial having a probability of success PR.
 *
 *   Prob(# failures = F | S successes, PR)  =
 *                       ( S + F - 1 )
 *                       (            ) * PR^S * (1-PR)^F
 *                       (      F     )
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.5.26.
 *
 * Parameters:
 *
 *   Input, double *F, the number of failures.
 *
 *   Input, double *S, the number of successes.
 *
 *   Input, double *PR, *OMPR, the probability of success on
 *   each binomial trial, and the value of (1-PR).
 *
 *   Output, double *CUM, *CCUM, the negative binomial CDF,
 *   and the complementary CDF.
 * </pre>
 */
void
cdflib_cumnbn(
    const double* s,
    const double* xn,
    const double* pr,
    const double* ompr,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMNOR computes the cumulative normal distribution.
 *
 * <pre>
 * Discussion:
 *
 *   This function evaluates the normal distribution function:
 *
 *                             / x
 *                    1       |       -t*t/2
 *         P(x) = ----------- |      e       dt
 *                sqrt(2 pi)  |
 *                            /-oo
 *
 *   This transportable program uses rational functions that
 *   theoretically approximate the normal distribution function to
 *   at least 18 significant decimal digits.  The accuracy achieved
 *   depends on the arithmetic system, the compiler, the intrinsic
 *   functions, and proper selection of the machine-dependent
 *   constants.
 *
 * Author:
 *
 *   William Cody
 *   Mathematics and Computer Science Division
 *   Argonne National Laboratory
 *   Argonne, IL 60439
 *
 * Reference:
 *
 *   William Cody,
 *   Rational Chebyshev approximations for the error function,
 *   Mathematics of Computation,
 *   1969, pages 631-637.
 *
 *   William Cody,
 *   Algorithm 715:
 *   SPECFUN - A Portable FORTRAN Package of Special Function Routines
 *     and Test Drivers,
 *   ACM Transactions on Mathematical Software,
 *   Volume 19, 1993, pages 22-32.
 *
 * Parameters:
 *
 *   Input, double *ARG, the upper limit of integration.
 *
 *   Output, double *CUM, *CCUM, the Normal density CDF and
 *   complementary CDF.
 *
 * Local Parameters:
 *
 *   Local, double EPS, the argument below which anorm(x)
 *   may be represented by 0.5D+00 and above which  x*x  will not underflow.
 *   A conservative value is the largest machine number X
 *   such that   1.0D+00 + X = 1.0D+00   to machine precision.
 * </pre>
 */
void
cdflib_cumnor(
    const double* arg,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMPOI evaluates the cumulative Poisson distribution.
 *
 * <pre>
 * Discussion:
 *
 *   CUMPOI returns the probability of S or fewer events in a Poisson
 *   distribution with mean XLAM.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   Formula 26.4.21.
 *
 * Parameters:
 *
 *   Input, double *S, the upper limit of cumulation of the
 *   Poisson density function.
 *
 *   Input, double *XLAM, the mean of the Poisson distribution.
 *
 *   Output, double *CUM, *CCUM, the Poisson density CDF and
 *   complementary CDF.
 * </pre>
 */
void
cdflib_cumpoi(
    double* s,
    double* xlam,
    double* cum,
    double* ccum);

/**
 * \brief
 *   CUMT evaluates the cumulative T distribution.
 *
 * <pre>
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   Formula 26.5.27.
 *
 * Parameters:
 *
 *   Input, double *T, the upper limit of integration.
 *
 *   Input, double *DF, the number of degrees of freedom of
 *   the T distribution.
 *
 *   Output, double *CUM, *CCUM, the T distribution CDF and
 *   complementary CDF.
 * </pre>
 */
void
cdflib_cumt(
    double* t,
    double* df,
    double* cum,
    double* ccum);

/**
 * \brief
 *   DBETRM computes the Sterling remainder for the complete beta function.
 *
 * <pre>
 * Discussion:
 *
 *   Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
 *   where Lgamma is the log of the (complete) gamma function
 *
 *   Let ZZ be approximation obtained if each log gamma is approximated
 *   by Sterling's formula, i.e.,
 *   Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5D+00 ) * LOG( Z ) - Z
 *
 *   The Sterling remainder is Log(Beta(A,B)) - ZZ.
 *
 * Parameters:
 *
 *   Input, double *A, *B, the parameters of the Beta function.
 *
 *   Output, double DBETRM, the Sterling remainder.
 * </pre>
 */
double
cdflib_dbetrm(
    double* a,
    double* b);

/**
 * \brief
 *   DEXPM1 evaluates the function EXP(X) - 1.
 *
 * <pre>
 * Reference:
 *
 *   Armido DiDinato and Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios,
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, 1993, pages 360-373.
 *
 * Parameters:
 *
 *   Input, double *X, the value at which exp(X)-1 is desired.
 *
 *   Output, double DEXPM1, the value of exp(X)-1.
 * </pre>
 */
double
cdflib_dexpm1(
    double* x);

/**
 * \brief
 *   DINVNR computes the inverse of the normal distribution.
 *
 * <pre>
 * Discussion:
 *
 *   Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
 *   infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
 *
 *   The rational function on page 95 of Kennedy and Gentle is used as a start
 *   value for the Newton method of finding roots.
 *
 * Reference:
 *
 *   Kennedy and Gentle,
 *   Statistical Computing,
 *   Marcel Dekker, NY, 1980,
 *   QA276.4  K46
 *
 * Parameters:
 *
 *   Input, double *P, *Q, the probability, and the complementary
 *   probability.
 *
 *   Output, double DINVNR, the argument X for which the
 *   Normal CDF has the value P.
 * </pre>
 */
double
cdflib_dinvnr(
    double* p,
    double* q);

/**
 * \brief
 *   DINVR bounds the zero of the function and invokes DZROR.
 *
 * <pre>
 * Discussion:
 *
 *   This routine seeks to find bounds on a root of the function and
 *   invokes ZROR to perform the zero finding.  STINVR must have been
 *   called before this routine in order to set its parameters.
 *
 * Reference:
 *
 *   J C P Bus and T J Dekker,
 *   Two Efficient Algorithms with Guaranteed Convergence for
 *     Finding a Zero of a Function,
 *   ACM Transactions on Mathematical Software,
 *   Volume 1, Number 4, pages 330-345, 1975.
 *
 * Parameters:
 *
 *   Input/output, integer STATUS.  At the beginning of a zero finding
 *   problem, STATUS should be set to 0 and INVR invoked.  The value
 *   of parameters other than X will be ignored on this call.
 *   If INVR needs the function to be evaluated, it will set STATUS to 1
 *   and return.  The value of the function should be set in FX and INVR
 *   again called without changing any of its other parameters.
 *   If INVR finishes without error, it returns with STATUS 0, and X an
 *   approximate root of F(X).
 *   If INVR cannot bound the function, it returns a negative STATUS and
 *   sets QLEFT and QHI.
 *
 *   Output, double precision X, the value at which F(X) is to be evaluated.
 *
 *   Input, double precision FX, the value of F(X) calculated by the user
 *   on the previous call, when INVR returned with STATUS = 1.
 *
 *   Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
 *   case, QLEFT is TRUE if the stepping search terminated unsucessfully
 *   at SMALL, and FALSE if the search terminated unsucessfully at BIG.
 *
 *   Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
 *   case, it is TRUE if Y < F(X) at the termination of the search and FALSE
 *   if F(X) < Y.
 * </pre>
 */
void
cdflib_dinvr(
    int* status,
    double* x,
    double* fx,
    unsigned long* qleft,
    unsigned long* qhi);

/**
 * \brief
 *   DLANOR evaluates the logarithm of the asymptotic Normal CDF.
 *
 * <pre>
 * Discussion:
 *
 *   This routine computes the logarithm of the cumulative normal distribution
 *   from abs ( x ) to infinity for  5 <= abs ( X ).
 *
 *   The relative error at X = 5 is about 0.5D-5.
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions
 *   1966, Formula 26.2.12.
 *
 * Parameters:
 *
 *   Input, double *X, the value at which the Normal CDF is to be
 *   evaluated.  It is assumed that 5 <= abs ( X ).
 *
 *   Output, double DLANOR, the logarithm of the asymptotic
 *   Normal CDF.
 * </pre>
 */
double
cdflib_dlanor(
    double* x);

/**
 * \brief
 *   DPMPAR provides machine constants for double precision arithmetic.
 *
 * <pre>
 * Discussion:
 *
 *    DPMPAR PROVIDES THE double PRECISION MACHINE CONSTANTS FOR
 *    THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
 *    I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
 *    double PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
 *    ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
 *
 *       DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
 *
 *       DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
 *
 *       DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
 *
 *    WRITTEN BY
 *       ALFRED H. MORRIS, JR.
 *       NAVAL SURFACE WARFARE CENTER
 *       DAHLGREN VIRGINIA
 *
 *    MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
 *    CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
 *    MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
 * </pre>
 */
double
cdflib_dpmpar(
    const int* i);

/**
 * \brief
 *   DSTINV seeks a value X such that F(X) = Y.
 *
 * <pre>
 * Discussion:
 *
 *     Double Precision - SeT INverse finder - Reverse Communication
 *                             Function
 *    Concise Description - Given a monotone function F finds X
 *    such that F(X) = Y.  Uses Reverse communication -- see invr.
 *    This routine sets quantities needed by INVR.
 *         More Precise Description of INVR -
 *    F must be a monotone function, the results of QMFINV are
 *    otherwise undefined.  QINCR must be .TRUE. if F is non-
 *    decreasing and .FALSE. if F is non-increasing.
 *    QMFINV will return .TRUE. if and only if F(SMALL) and
 *    F(BIG) bracket Y, i. e.,
 *         QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
 *         QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
 *    if QMFINV returns .TRUE., then the X returned satisfies
 *    the following condition.  let
 *              TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
 *    then if QINCR is .TRUE.,
 *         F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
 *    and if QINCR is .FALSE.
 *         F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
 *                             Arguments
 *    SMALL --> The left endpoint of the interval to be
 *         searched for a solution.
 *                   SMALL is DOUBLE PRECISION
 *    BIG --> The right endpoint of the interval to be
 *         searched for a solution.
 *                   BIG is DOUBLE PRECISION
 *    ABSSTP, RELSTP --> The initial step size in the search
 *         is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
 *                   ABSSTP is DOUBLE PRECISION
 *                   RELSTP is DOUBLE PRECISION
 *    STPMUL --> When a step doesn't bound the zero, the step
 *               size is multiplied by STPMUL and another step
 *               taken.  A popular value is 2.0
 *                   DOUBLE PRECISION STPMUL
 *    ABSTOL, RELTOL --> Two numbers that determine the accuracy
 *         of the solution.  See function for a precise definition.
 *                   ABSTOL is DOUBLE PRECISION
 *                   RELTOL is DOUBLE PRECISION
 *                             Method
 *    Compares F(X) with Y for the input value of X then uses QINCR
 *    to determine whether to step left or right to bound the
 *    desired x.  the initial step size is
 *         MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
 *    Iteratively steps right or left until it bounds X.
 *    At each step which doesn't bound X, the step size is doubled.
 *    The routine is careful never to step beyond SMALL or BIG.  If
 *    it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
 *    after setting QLEFT and QHI.
 *    If X is successfully bounded then Algorithm R of the paper
 *    'Two Efficient Algorithms with Guaranteed Convergence for
 *    Finding a Zero of a Function' by J. C. P. Bus and
 *    T. J. Dekker in ACM Transactions on Mathematical
 *    Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
 *    to find the zero of the function F(X)-Y. This is routine
 *    QRZERO.
 * </pre>
 */
void
cdflib_dstinv(
    double* zsmall,
    double* zbig,
    double* zabsst,
    double* zrelst,
    double* zstpmu,
    double* zabsto,
    double* zrelto);

/**
 * \brief
 *   DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
 *
 * <pre>
 * Discussion:
 *
 *   This routine returns
 *
 *     ln ( Gamma ( Z ) ) - Sterling ( Z )
 *
 *   where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
 *
 *   Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
 *
 *   If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
 *   with values calculated using Maple.
 *
 *   Otherwise, the difference is computed explicitly.
 *
 * Modified:
 *
 *   14 June 2004
 *
 * Parameters:
 *
 *   Input, double *Z, the value at which the Sterling
 *   remainder is to be calculated.  Z must be positive.
 *
 *   Output, double DSTREM, the Sterling remainder.
 * </pre>
 */
double
cdflib_dstrem(
    double* z);

/**
 * \brief
 *   DSTXR sets quantities needed by the zero finder.
 *
 * <pre>
 * Discussion:
 *
 *    Double precision SeT ZeRo finder - Reverse communication version
 *                             Function
 *    Sets quantities needed by ZROR.  The function of ZROR
 *    and the quantities set is given here.
 *    Concise Description - Given a function F
 *    find XLO such that F(XLO) = 0.
 *         More Precise Description -
 *    Input condition. F is a double function of a single
 *    double argument and XLO and XHI are such that
 *         F(XLO)*F(XHI)  .LE.  0.0
 *    If the input condition is met, QRZERO returns .TRUE.
 *    and output values of XLO and XHI satisfy the following
 *         F(XLO)*F(XHI)  .LE. 0.
 *         ABS(F(XLO)  .LE. ABS(F(XHI)
 *         ABS(XLO-XHI)  .LE. TOL(X)
 *    where
 *         TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
 *    If this algorithm does not find XLO and XHI satisfying
 *    these conditions then QRZERO returns .FALSE.  This
 *    implies that the input condition was not met.
 *                             Arguments
 *    XLO --> The left endpoint of the interval to be
 *          searched for a solution.
 *                   XLO is DOUBLE PRECISION
 *    XHI --> The right endpoint of the interval to be
 *          for a solution.
 *                   XHI is DOUBLE PRECISION
 *    ABSTOL, RELTOL --> Two numbers that determine the accuracy
 *                     of the solution.  See function for a
 *                     precise definition.
 *                   ABSTOL is DOUBLE PRECISION
 *                   RELTOL is DOUBLE PRECISION
 *                             Method
 *    Algorithm R of the paper 'Two Efficient Algorithms with
 *    Guaranteed Convergence for Finding a Zero of a Function'
 *    by J. C. P. Bus and T. J. Dekker in ACM Transactions on
 *    Mathematical Software, Volume 1, no. 4 page 330
 *    (Dec. '75) is employed to find the zero of F(X)-Y.
 * </pre>
 */
void
cdflib_dstzr(
    double* zxlo,
    double* zxhi,
    double* zabstl,
    double* zreltl);

/**
 * \brief
 *   DT1 computes an approximate inverse of the cumulative T distribution.
 *
 * <pre>
 * Discussion:
 *
 *   Returns the inverse of the T distribution function, i.e.,
 *   the integral from 0 to INVT of the T density is P. This is an
 *   initial approximation.
 *
 * Parameters:
 *
 *   Input, double *P, *Q, the value whose inverse from the
 *   T distribution CDF is desired, and the value (1-P).
 *
 *   Input, double *DF, the number of degrees of freedom of the
 *   T distribution.
 *
 *   Output, double DT1, the approximate value of X for which
 *   the T density CDF with DF degrees of freedom has value P.
 * </pre>
 */
double
cdflib_dt1(
    double* p,
    double* q,
    double* df);

/**
 * \brief
 *   DZROR seeks the zero of a function using reverse communication.
 *
 * <pre>
 * Discussion:
 *
 *    Performs the zero finding.  STZROR must have been called before
 *    this routine in order to set its parameters.
 *
 *
 *                             Arguments
 *
 *
 *    STATUS <--> At the beginning of a zero finding problem, STATUS
 *                should be set to 0 and ZROR invoked.  (The value
 *                of other parameters will be ignored on this call.)
 *
 *                When ZROR needs the function evaluated, it will set
 *                STATUS to 1 and return.  The value of the function
 *                should be set in FX and ZROR again called without
 *                changing any of its other parameters.
 *
 *                When ZROR has finished without error, it will return
 *                with STATUS 0.  In that case (XLO,XHI) bound the answe
 *
 *                If ZROR finds an error (which implies that F(XLO)-Y an
 *                F(XHI)-Y have the same sign, it returns STATUS -1.  In
 *                this case, XLO and XHI are undefined.
 *                        INTEGER STATUS
 *
 *    X <-- The value of X at which F(X) is to be evaluated.
 *                        DOUBLE PRECISION X
 *
 *    FX --> The value of F(X) calculated when ZROR returns with
 *           STATUS = 1.
 *                        DOUBLE PRECISION FX
 *
 *    XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
 *            inverval in X containing the solution below.
 *                        DOUBLE PRECISION XLO
 *
 *    XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
 *            inverval in X containing the solution above.
 *                        DOUBLE PRECISION XHI
 *
 *    QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
 *               at XLO.  If it is .FALSE. the search terminated
 *               unsucessfully at XHI.
 *                   QLEFT is LOGICAL
 *
 *    QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
 *             search and .FALSE. if F(X) .LT. Y at the
 *             termination of the search.
 *                   QHI is LOGICAL
 *
 * </pre>
 */
void
cdflib_dzror(
    int* status,
    double* x,
    double* fx,
    double* xlo,
    double* xhi,
    unsigned long* qleft,
    unsigned long* qhi);

/**
 * \brief
 *   E0000 is a reverse-communication zero bounder.
 */
void
cdflib_E0000(
    int IENTRY,
    int* status,
    double* x,
    double* fx,
    unsigned long* qleft,
    unsigned long* qhi,
    double* zabsst,
    double* zabsto,
    double* zbig,
    double* zrelst,
    double* zrelto,
    double* zsmall,
    double* zstpmu);

/**
 * \brief
 *   E00001 is a reverse-communication zero finder.
 */
void
cdflib_E0001(
    int IENTRY,
    int* status,
    double* x,
    double* fx,
    double* xlo,
    double* xhi,
    unsigned long* qleft,
    unsigned long* qhi,
    double* zabstl,
    double* zreltl,
    double* zxhi,
    double* zxlo);

/**
 * \brief
 *   ERF_VALUES returns some values of the ERF or "error" function.
 *
 * <pre>
 * Definition:
 *
 *   ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_erf_values(
    int* n_data,
    double* x,
    double* fx);

/**
 * \brief
 *   ERROR_F evaluates the error function ERF.
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *X, the argument.
 *
 *   Output, double ERROR_F, the value of the error function at X.
 * </pre>
 */
double
cdflib_error_f(
    double* x);

/**
 * \brief
 *   ERROR_FC evaluates the complementary error function ERFC.
 *
 * <pre>
 * Modified:
 *
 *   09 December 1999
 *
 * Parameters:
 *
 *   Input, int *IND, chooses the scaling.
 *   If IND is nonzero, then the value returned has been multiplied by
 *   EXP(X*X).
 *
 *   Input, double *X, the argument of the function.
 *
 *   Output, double ERROR_FC, the value of the complementary
 *   error function.
 * </pre>
 */
double
cdflib_error_fc(
    const int* ind,
    double* x);

/**
 * \brief
 *   ESUM evaluates exp ( MU + X ).
 *
 * <pre>
 * Parameters:
 *
 *   Input, int *MU, part of the argument.
 *
 *   Input, double *X, part of the argument.
 *
 *   Output, double ESUM, the value of exp ( MU + X ).
 * </pre>
 */
double
cdflib_esum(
    const int* mu,
    const double* x);

/**
 * \brief
 *   EVAL_POL evaluates a polynomial at X.
 *
 * <pre>
 * Discussion:
 *
 *   EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
 *
 * Modified:
 *
 *   15 December 1999
 *
 * Parameters:
 *
 *   Input, double precision A(0:N), coefficients of the polynomial.
 *
 *   Input, int *N, length of A.
 *
 *   Input, double *X, the point at which the polynomial
 *   is to be evaluated.
 *
 *   Output, double EVAL_POL, the value of the polynomial at X.
 * </pre>
 */
double
cdflib_eval_pol(
    const double a[],
    const int* n,
    double* x);

/**
 * \brief
 *   EXPARG returns the largest or smallest legal argument for EXP.
 *
 * <pre>
 * Discussion:
 *
 *   Only an approximate limit for the argument of EXP is desired.
 *
 * Modified:
 *
 *   09 December 1999
 *
 * Parameters:
 *
 *   Input, int *L, indicates which limit is desired.
 *   If L = 0, then the largest positive argument for EXP is desired.
 *   Otherwise, the largest negative argument for EXP for which the
 *   result is nonzero is desired.
 *
 *   Output, double EXPARG, the desired value.
 * </pre>
 */
double
cdflib_exparg(
    const int* l);

/**
 * \brief
 *   F_CDF_VALUES returns some values of the F CDF test function.
 *
 * <pre>
 * Discussion:
 *
 *   The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
 *   commands like:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF[FRatioDistribution[ DFN, DFD ], X ]
 *
 * Modified:
 *
 *   11 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *A, int *B, the parameters of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_f_cdf_values(
    int* n_data,
    int* a,
    int* b,
    double* x,
    double* fx);

/**
 * \brief
 *   F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
 *
 * <pre>
 * Discussion:
 *
 *   The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
 *   in Mathematica by commands like:
 *
 *     Needs["Statistics`ContinuousDistributions`"]
 *     CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
 *
 * Modified:
 *
 *   12 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Stephen Wolfram,
 *   The Mathematica Book,
 *   Fourth Edition,
 *   Wolfram Media / Cambridge University Press, 1999.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *A, int *B, double *LAMBDA, the
 *   parameters of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_f_noncentral_cdf_values(
    int* n_data,
    int* a,
    int* b,
    double* lambda,
    double* x,
    double* fx);

/**
 * \brief
 *   FIFDINT truncates a double number to an integer.
 *
 * <pre>
 * Parameters:
 *
 *   a     -     number to be truncated
 * </pre>
 */
double
cdflib_fifdint(
    double a);

/**
 * \brief
 *   FIFDMAX1 returns the maximum of two numbers a and b
 *
 * <pre>
 * Parameters:
 *
 *   a     -      first number
 *   b     -      second number
 * </pre>
 */
double
cdflib_fifdmax1(
    double a,
    double b);

/**
 * \brief
 *   FIFDMIN1 returns the minimum of two numbers.
 *
 * <pre>
 * Parameters:
 *
 *   a     -     first number
 *   b     -     second number
 * </pre>
 */
double
cdflib_fifdmin1(
    double a,
    double b);

/**
 * \brief
 *   FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
 *
 * <pre>
 * Parameters:
 *
 *   mag     -     magnitude
 *   sign    -     sign to be transfered
 * </pre>
 */
double
cdflib_fifdsign(
    double mag,
    double sign);

/**
 * \brief
 *   FIFIDINT truncates a double number to a long integer
 *
 * <pre>
 * Parameters:
 *
 *   a - number to be truncated
 * </pre>
 */
long
cdflib_fifidint(
    double a);

/**
 * \brief
 *   FIFMOD returns the modulo of a and b
 *
 * <pre>
 * Parameters:
 *
 *   a - numerator
 *   b - denominator
 * </pre>
 */
long
cdflib_fifmod(
    long a,
    long b);

/**
 * \brief
 *   FPSER evaluates IX(A,B)(X) for very small B.
 *
 * <pre>
 * Discussion:
 *
 *   This routine is appropriate for use when
 *
 *     B < min ( EPS, EPS * A )
 *
 *   and
 *
 *     X <= 0.5.
 *
 * Parameters:
 *
 *   Input, double *A, *B, parameters of the function.
 *
 *   Input, double *X, the point at which the function is to
 *   be evaluated.
 *
 *   Input, double *EPS, a tolerance.
 *
 *   Output, double FPSER, the value of IX(A,B)(X).
 * </pre>
 */
double
cdflib_fpser(
    double* a,
    double* b,
    double* x,
    double* eps);

/**
 * \brief
 *   FTNSTOP prints a message to standard error and then exits.
 *
 * <pre>
 * Parameters:
 *
 *   Input, char *MSG, the message to be printed.
 * </pre>
 */
void
cdflib_ftnstop(
    char* msg);

/**
 * \brief
 *   GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, forms the argument of the Gamma function.
 *
 *   Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
 * </pre>
 */
double
cdflib_gam1(
    const double* a);

/**
 * \brief
 *   GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
 *
 * <pre>
 * Discussion:
 *
 *   This is certified spaghetti code.
 *
 * Author:
 *
 *   Alfred H Morris, Jr,
 *   Naval Surface Weapons Center,
 *   Dahlgren, Virginia.
 *
 * Parameters:
 *
 *   Input, double *A, *X, the arguments of the incomplete
 *   gamma ratio.  A and X must be nonnegative.  A and X cannot
 *   both be zero.
 *
 *   Output, double *ANS, *QANS.  On normal output,
 *   ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
 *   A or X is negative, or both are 0, or when the answer is
 *   computationally indeterminate because A is extremely large
 *   and X is very close to A.
 *
 *   Input, int *IND, indicates the accuracy request:
 *   0, as much accuracy as possible.
 *   1, to within 1 unit of the 6-th significant digit,
 *   otherwise, to within 1 unit of the 3rd significant digit.
 * </pre>
 */
void
cdflib_gamma_inc(
    const double* a,
    const double* x,
    double* ans,
    double* qans,
    const int* ind);

/**
 * \brief
 *   GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
 *
 * <pre>
 * Discussion:
 *
 *   The routine is given positive A, and nonnegative P and Q where P + Q = 1.
 *   The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
 *   Schroder iteration is employed.  The routine attempts to compute X
 *   to 10 significant digits if this is possible for the particular computer
 *   arithmetic being used.
 *
 * Author:
 *
 *   Alfred H Morris, Jr,
 *   Naval Surface Weapons Center,
 *   Dahlgren, Virginia.
 *
 * Parameters:
 *
 *   Input, double *A, the parameter in the incomplete gamma
 *   ratio.  A must be positive.
 *
 *   Output, double *X, the computed point for which the
 *   incomplete gamma functions have the values P and Q.
 *
 *   Input, double *X0, an optional initial approximation
 *   for the solution X.  If the user does not want to supply an
 *   initial approximation, then X0 should be set to 0, or a negative
 *   value.
 *
 *   Input, double *P, *Q, the values of the incomplete gamma
 *   functions, for which the corresponding argument is desired.
 *
 *   Output, int *IERR, error flag.
 *   0, the solution was obtained. Iteration was not used.
 *   0 < K, The solution was obtained. IERR iterations were performed.
 *   -2, A <= 0
 *   -3, No solution was obtained. The ratio Q/A is too large.
 *   -4, P + Q /= 1
 *   -6, 20 iterations were performed. The most recent value obtained
 *       for X is given.  This cannot occur if X0 <= 0.
 *   -7, Iteration failed. No value is given for X.
 *       This may occur when X is approximately 0.
 *   -8, A value for X has been obtained, but the routine is not certain
 *       of its accuracy.  Iteration cannot be performed in this
 *       case. If X0 <= 0, this can occur only when P or Q is
 *       approximately 0. If X0 is positive then this can occur when A is
 *       exceedingly close to X and A is extremely large (say A .GE. 1.E20).
 * </pre>
 */
void
cdflib_gamma_inc_inv(
    double* a,
    double* x,
    double* x0,
    double* p,
    double* q,
    int* ierr);

/**
 * \brief
 *   GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
 *
 * <pre>
 * Discussion:
 *
 *   The (normalized) incomplete Gamma function P(A,X) is defined as:
 *
 *     PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
 *
 *   With this definition, for all A and X,
 *
 *     0 <= PN(A,X) <= 1
 *
 *   and
 *
 *     PN(A,INFINITY) = 1.0
 *
 *   Mathematica can compute this value as
 *
 *     1 - GammaRegularized[A,X]
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *A, the parameter of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_gamma_inc_values(
    int* n_data,
    double* a,
    double* x,
    double* fx);

/**
 * \brief
 *   GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, defines the argument of the function.
 *
 *   Output, double GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
 * </pre>
 */
double
cdflib_gamma_ln1(
    double* a);

/**
 * \brief
 *   GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
 *
 * <pre>
 * Author:
 *
 *   Alfred H Morris, Jr,
 *   Naval Surface Weapons Center,
 *   Dahlgren, Virginia.
 *
 * Reference:
 *
 *   Armido DiDinato and Alfred Morris,
 *   Algorithm 708:
 *   Significant Digit Computation of the Incomplete Beta Function Ratios,
 *   ACM Transactions on Mathematical Software,
 *   Volume 18, 1993, pages 360-373.
 *
 * Parameters:
 *
 *   Input, double *A, the argument of the function.
 *   A should be positive.
 *
 *   Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
 * </pre>
 */
double
cdflib_gamma_log(
    double* a);

/**
 * \brief
 *   GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *X, the parameters of the functions.
 *   It is assumed that A <= 1.
 *
 *   Input, double *R, the value exp(-X) * X**A / Gamma(A).
 *
 *   Output, double *P, *Q, the values of P(A,X) and Q(A,X).
 *
 *   Input, double *EPS, the tolerance.
 * </pre>
 */
void
cdflib_gamma_rat1(
    const double* a,
    const double* x,
    const double* r,
    double* p,
    double* q,
    const double* eps);

/**
 * \brief
 *   GAMMA_VALUES returns some values of the Gamma function.
 *
 * <pre>
 * Definition:
 *
 *   GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
 *
 * Recursion:
 *
 *   GAMMA(X+1) = X*GAMMA(X)
 *
 * Restrictions:
 *
 *   0 < X ( a software restriction).
 *
 * Special values:
 *
 *   GAMMA(0.5) = sqrt(PI)
 *
 *   For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_gamma_values(
    int* n_data,
    double* x,
    double* fx);

/**
 * \brief
 *   GAMMA_X evaluates the gamma function.
 *
 * <pre>
 * Discussion:
 *
 *   This routine was renamed from "GAMMA" to avoid a conflict with the
 *   C/C++ math library routine.
 *
 * Author:
 *
 *   Alfred H Morris, Jr,
 *   Naval Surface Weapons Center,
 *   Dahlgren, Virginia.
 *
 * Parameters:
 *
 *   Input, double *A, the argument of the Gamma function.
 *
 *   Output, double GAMMA_X, the value of the Gamma function.
 * </pre>
 */
double
cdflib_gamma_x(
    const double* a);

/**
 * \brief
 *   GSUMLN evaluates the function ln(Gamma(A + B)).
 *
 * <pre>
 * Discussion:
 *
 *   GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
 *
 * Parameters:
 *
 *   Input, double *A, *B, values whose sum is the argument of
 *   the Gamma function.
 *
 *   Output, double GSUMLN, the value of ln(Gamma(A+B)).
 * </pre>
 */
double
cdflib_gsumln(
    double* a,
    double* b);

/**
 * \brief
 *   IPMPAR returns integer machine constants.
 *
 * <pre>
 * Discussion:
 *
 *   Input arguments 1 through 3 are queries about integer arithmetic.
 *   We assume integers are represented in the N-digit, base-A form
 *
 *     sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
 *
 *   where 0 <= X(0:N-1) < A.
 *
 *   Then:
 *
 *     IPMPAR(1) = A, the base of integer arithmetic;
 *     IPMPAR(2) = N, the number of base A digits;
 *     IPMPAR(3) = A**N - 1, the largest magnitude.
 *
 *   It is assumed that the single and double precision floating
 *   point arithmetics have the same base, say B, and that the
 *   nonzero numbers are represented in the form
 *
 *     sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
 *
 *   where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
 *   EMIN <= E <= EMAX.
 *
 *   Input argument 4 is a query about the base of real arithmetic:
 *
 *     IPMPAR(4) = B, the base of single and double precision arithmetic.
 *
 *   Input arguments 5 through 7 are queries about single precision
 *   floating point arithmetic:
 *
 *    IPMPAR(5) = M, the number of base B digits for single precision.
 *    IPMPAR(6) = EMIN, the smallest exponent E for single precision.
 *    IPMPAR(7) = EMAX, the largest exponent E for single precision.
 *
 *   Input arguments 8 through 10 are queries about double precision
 *   floating point arithmetic:
 *
 *    IPMPAR(8) = M, the number of base B digits for double precision.
 *    IPMPAR(9) = EMIN, the smallest exponent E for double precision.
 *    IPMPAR(10) = EMAX, the largest exponent E for double precision.
 *
 * Reference:
 *
 *   Phyllis Fox, Andrew Hall, and Norman Schryer,
 *   Algorithm 528,
 *   Framework for a Portable FORTRAN Subroutine Library,
 *   ACM Transactions on Mathematical Software,
 *   Volume 4, 1978, pages 176-188.
 *
 * Parameters:
 *
 *   Input, int *I, the index of the desired constant.
 *
 *   Output, int IPMPAR, the value of the desired constant.
 * </pre>
 */
int
cdflib_ipmpar(
    const int* i);

/**
 * \brief
 *   NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
 *
 * <pre>
 * Discussion:
 *
 *   Assume that a coin has a probability P of coming up heads on
 *   any one trial.  Suppose that we plan to flip the coin until we
 *   achieve a total of S heads.  If we let F represent the number of
 *   tails that occur in this process, then the value of F satisfies
 *   a negative binomial PDF:
 *
 *     PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
 *
 *   The negative binomial CDF is the probability that there are F or
 *   fewer failures upon the attainment of the S-th success.  Thus,
 *
 *     CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
 *
 * Modified:
 *
 *   07 June 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   F C Powell,
 *   Statistical Tables for Sociology, Biology and Physical Sciences,
 *   Cambridge University Press, 1982.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *F, the maximum number of failures.
 *
 *   Output, int *S, the number of successes.
 *
 *   Output, double *P, the probability of a success on one trial.
 *
 *   Output, double *CDF, the probability of at most F failures before the
 *   S-th success.
 * </pre>
 */
void
cdflib_negative_binomial_cdf_values(
    int* n_data,
    int* f,
    int* s,
    double* p,
    double* cdf);

/**
 * \brief
 *   NORMAL_CDF_VALUES returns some values of the Normal CDF.
 *
 * <pre>
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output double *FX, the value of the function.
 * </pre>
 */
void
cdflib_normal_cdf_values(
    int* n_data,
    double* x,
    double* fx);

/**
 * \brief
 *   POISSON_CDF_VALUES returns some values of the Poisson CDF.
 *
 * <pre>
 * Discussion:
 *
 *   CDF(X)(A) is the probability of at most X successes in unit time,
 *   given that the expected mean number of successes is A.
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 *   Daniel Zwillinger,
 *   CRC Standard Mathematical Tables and Formulae,
 *   30th Edition, CRC Press, 1996, pages 653-658.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *A, the parameter of the function.
 *
 *   Output, int *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_poisson_cdf_values(
    int* n_data,
    double* a,
    int* x,
    double* fx);

/**
 * \brief
 *   PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
 *
 * <pre>
 * Discussion:
 *
 *   The main computation involves evaluation of rational Chebyshev
 *   approximations.  PSI was written at Argonne National Laboratory
 *   for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
 *
 * Reference:
 *
 *   William Cody, Strecok and Thacher,
 *   Chebyshev Approximations for the Psi Function,
 *   Mathematics of Computation,
 *   Volume 27, 1973, pages 123-127.
 *
 * Parameters:
 *
 *   Input, double *XX, the argument of the psi function.
 *
 *   Output, double PSI, the value of the psi function.  PSI
 *   is assigned the value 0 when the psi function is undefined.
 * </pre>
 */
double
cdflib_psi(
    const double* xx);

/**
 * \brief
 *   PSI_VALUES returns some values of the Psi or Digamma function.
 *
 * <pre>
 * Discussion:
 *
 *   PSI(X) = d LN ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
 *
 *   PSI(1) = - Euler's constant.
 *
 *   PSI(X+1) = PSI(X) + 1 / X.
 *
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_psi_values(
    int* n_data,
    double* x,
    double* fx);

/**
 * \brief
 *   RCOMP evaluates exp(-X) * X**A / Gamma(A).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *A, *X, arguments of the quantity to be computed.
 *
 *   Output, double RCOMP, the value of exp(-X) * X**A / Gamma(A).
 *
 * Local parameters:
 *
 *   RT2PIN = 1/SQRT(2*PI)
 * </pre>
 */
double
cdflib_rcomp(
    double* a,
    double* x);

/**
 * \brief
 *   REXP evaluates the function EXP(X) - 1.
 *
 * <pre>
 * Modified:
 *
 *   09 December 1999
 *
 * Parameters:
 *
 *   Input, double *X, the argument of the function.
 *
 *   Output, double REXP, the value of EXP(X)-1.
 * </pre>
 */
double
cdflib_rexp(
    double* x);

/**
 * \brief
 *   RLOG computes  X - 1 - LN(X).
 *
 * <pre>
 * Modified:
 *
 *   09 December 1999
 *
 * Parameters:
 *
 *   Input, double *X, the argument of the function.
 *
 *   Output, double RLOG, the value of the function.
 * </pre>
 */
double
cdflib_rlog(
    double* x);

/**
 * \brief
 *   RLOG1 evaluates the function X - ln ( 1 + X ).
 *
 * <pre>
 * Parameters:
 *
 *   Input, double *X, the argument.
 *
 *   Output, double RLOG1, the value of X - ln ( 1 + X ).
 * </pre>
 */
double
cdflib_rlog1(
    double* x);

/**
 * \brief
 *   STUDENT_CDF_VALUES returns some values of the Student CDF.
 *
 * <pre>
 * Modified:
 *
 *   31 May 2004
 *
 * Author:
 *
 *   John Burkardt
 *
 * Reference:
 *
 *   Milton Abramowitz and Irene Stegun,
 *   Handbook of Mathematical Functions,
 *   US Department of Commerce, 1964.
 *
 * Parameters:
 *
 *   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
 *   first call.  On each call, the routine increments N_DATA by 1, and
 *   returns the corresponding data; when there is no more data, the
 *   output value of N_DATA will be 0 again.
 *
 *   Output, int *A, the parameter of the function.
 *
 *   Output, double *X, the argument of the function.
 *
 *   Output, double *FX, the value of the function.
 * </pre>
 */
void
cdflib_student_cdf_values(
    int* n_data,
    int* a,
    double* x,
    double* fx);

/**
 * \brief
 *   STVALN provides starting values for the inverse of the normal distribution.
 *
 * <pre>
 * Discussion:
 *
 *   The routine returns X such that
 *     P = CUMNOR(X),
 *   that is,
 *     P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
 *
 * Reference:
 *
 *   Kennedy and Gentle,
 *   Statistical Computing,
 *   Marcel Dekker, NY, 1980, page 95,
 *   QA276.4  K46
 *
 * Parameters:
 *
 *   Input, double *P, the probability whose normal deviate
 *   is sought.
 *
 *   Output, double STVALN, the normal deviate whose probability
 *   is P.
 * </pre>
 */
double
cdflib_stvaln(
    double* p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_CDFLIB_H */
