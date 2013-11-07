/*
 * Modified from CDFLIB: http://www.netlib.org/random/dcdflib.c.tar.gz
 * Obtained from http://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html
 *
 * Barry Brown, James Lovato, Kathy Russell,
 * Department of Biomathematics,
 * University of Texas,
 * Houston, Texas.
 */

#ifndef TUNA_CDFLIB_H
#define TUNA_CDFLIB_H

/** @file
 * A modified version of CDFLIB by
 * Barry Brown, James Lovato, and Kathy Russell.
 */

double cdflib_algdiv ( double *a, double *b );
double cdflib_alnrel ( double *a );
double cdflib_apser ( double *a, double *b, double *x, double *eps );
double cdflib_bcorr ( double *a0, double *b0 );
double cdflib_beta ( double a, double b );
double cdflib_beta_asym ( double *a, double *b, double *lambda, double *eps );
double cdflib_beta_frac ( double *a, double *b, double *x, double *y, double *lambda,
  double *eps );
void cdflib_beta_grat ( double *a, double *b, double *x, double *y, double *w,
  double *eps,int *ierr );
void cdflib_beta_inc ( double *a, double *b, double *x, double *y, double *w,
  double *w1, int *ierr );
void cdflib_beta_inc_values ( int *n_data, double *a, double *b, double *x, double *fx );
double cdflib_beta_log ( double *a0, double *b0 );
double cdflib_beta_pser ( double *a, double *b, double *x, double *eps );
double cdflib_beta_cdflib_rcomp ( double *a, double *b, double *x, double *y );
double cdflib_beta_rcomp1 ( int *mu, double *a, double *b, double *x, double *y );
double cdflib_beta_up ( double *a, double *b, double *x, double *y, int *n, double *eps );
void cdflib_binomial_cdf_values ( int *n_data, int *a, double *b, int *x, double *fx );
void cdflib_cdfbet ( int *which, double *p, double *q, double *x, double *y,
  double *a, double *b, int *status, double *bound );
void cdflib_cdfbin ( int *which, double *p, double *q, double *s, double *xn,
  double *pr, double *ompr, int *status, double *bound );
void cdflib_cdfchi ( int *which, double *p, double *q, double *x, double *df,
  int *status, double *bound );
void cdflib_cdfchn ( int *which, double *p, double *q, double *x, double *df,
  double *pnonc, int *status, double *bound );
void cdflib_cdff ( int *which, double *p, double *q, double *f, double *dfn,
  double *dfd, int *status, double *bound );
void cdflib_cdffnc ( int *which, double *p, double *q, double *f, double *dfn,
  double *dfd, double *phonc, int *status, double *bound );
void cdflib_cdfgam ( int *which, double *p, double *q, double *x, double *shape,
  double *scale, int *status, double *bound );
void cdflib_cdfnbn ( int *which, double *p, double *q, double *s, double *xn,
  double *pr, double *ompr, int *status, double *bound );
void cdflib_cdfnor ( int *which, double *p, double *q, double *x, double *mean,
  double *sd, int *status, double *bound );
void cdflib_cdfpoi ( int *which, double *p, double *q, double *s, double *xlam,
  int *status, double *bound );
void cdflib_cdft ( int *which, double *p, double *q, double *t, double *df,
  int *status, double *bound );
void cdflib_chi_noncentral_cdf_values ( int *n_data, double *x, double *lambda, 
  int *df, double *cdf );
void cdflib_chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx );
void cdflib_cumbet ( double *x, double *y, double *a, double *b, double *cum,
  double *ccum );
void cdflib_cumbin ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum );
void cdflib_cumchi ( double *x, double *df, double *cum, double *ccum );
void cdflib_cumchn ( double *x, double *df, double *pnonc, double *cum,
  double *ccum );
void cdflib_cumf ( double *f, double *dfn, double *dfd, double *cum, double *ccum );
void cdflib_cumfnc ( double *f, double *dfn, double *dfd, double *pnonc,
  double *cum, double *ccum );
void cdflib_cumgam ( double *x, double *a, double *cum, double *ccum );
void cdflib_cumnbn ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum );
void cdflib_cumnor ( double *arg, double *result, double *ccum );
void cdflib_cumpoi ( double *s, double *xlam, double *cum, double *ccum );
void cdflib_cumt ( double *t, double *df, double *cum, double *ccum );
double cdflib_dbetrm ( double *a, double *b );
double cdflib_dexpm1 ( double *x );
double cdflib_dinvnr ( double *p, double *q );
void cdflib_dinvr ( int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi );
double cdflib_dlanor ( double *x );
double cdflib_dpmpar ( int *i );
void cdflib_dstinv ( double *zsmall, double *zbig, double *zabsst,
  double *zrelst, double *zstpmu, double *zabsto, double *zrelto );
double cdflib_dstrem ( double *z );
void cdflib_dstzr ( double *zxlo, double *zxhi, double *zabstl, double *zreltl );
double cdflib_dt1 ( double *p, double *q, double *df );
void cdflib_dzror ( int *status, double *x, double *fx, double *xlo,
  double *xhi, unsigned long *qleft, unsigned long *qhi );
void cdflib_E0000 ( int IENTRY, int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi, double *zabsst,
  double *zabsto, double *zbig, double *zrelst,
  double *zrelto, double *zsmall, double *zstpmu );
void cdflib_E0001 ( int IENTRY, int *status, double *x, double *fx,
  double *xlo, double *xhi, unsigned long *qleft,
  unsigned long *qhi, double *zabstl, double *zreltl,
  double *zxhi, double *zxlo );
void cdflib_erf_values ( int *n_data, double *x, double *fx );
double cdflib_error_f ( double *x );
double cdflib_error_fc ( int *ind, double *x );
double cdflib_esum ( int *mu, double *x );
double cdflib_eval_pol ( double a[], int *n, double *x );
double cdflib_exparg ( int *l );
void cdflib_f_cdf_values ( int *n_data, int *a, int *b, double *x, double *fx );
void cdflib_f_noncentral_cdf_values ( int *n_data, int *a, int *b, double *lambda, 
  double *x, double *fx );
double cdflib_fifdint ( double a );
double cdflib_fifdmax1 ( double a, double b );
double cdflib_fifdmin1 ( double a, double b );
double cdflib_fifdsign ( double mag, double sign );
long cdflib_fifidint ( double a );
long cdflib_fifmod ( long a, long b );
double cdflib_fpser ( double *a, double *b, double *x, double *eps );
void cdflib_ftnstop ( char *msg );
double cdflib_gam1 ( double *a );
void cdflib_gamma_inc ( double *a, double *x, double *ans, double *qans, int *ind );
void cdflib_gamma_inc_inv ( double *a, double *x, double *x0, double *p, double *q,
  int *ierr );
void cdflib_gamma_inc_values ( int *n_data, double *a, double *x, double *fx );
double cdflib_gamma_ln1 ( double *a );
double cdflib_gamma_log ( double *a );
void cdflib_gamma_rat1 ( double *a, double *x, double *r, double *p, double *q,
  double *eps );
void cdflib_gamma_values ( int *n_data, double *x, double *fx );
double cdflib_gamma_x ( double *a );
double cdflib_gsumln ( double *a, double *b );
int cdflib_ipmpar ( int *i );
void negative_cdflib_binomial_cdf_values ( int *n_data, int *f, int *s, double *p, 
  double *cdf );
void cdflib_normal_cdf_values ( int *n_data, double *x, double *fx );
void cdflib_poisson_cdf_values ( int *n_data, double *a, int *x, double *fx );
double cdflib_psi ( double *xx );
void cdflib_psi_values ( int *n_data, double *x, double *fx );
double cdflib_rcomp ( double *a, double *x );
double cdflib_rexp ( double *x );
double cdflib_rlog ( double *x );
double cdflib_rlog1 ( double *x );
void cdflib_student_cdf_values ( int *n_data, int *a, double *x, double *fx );
double cdflib_stvaln ( double *p );
void cdflib_timestamp ( void );

#endif /* TUNA_CDFLIB_H */
