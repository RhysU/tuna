/*
 * Modified from CDFLIB: http://www.netlib.org/random/dcdflib.c.tar.gz
 * Obtained from http://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html
 *
 * Barry Brown, James Lovato, Kathy Russell,
 * Department of Biomathematics,
 * University of Texas,
 * Houston, Texas.
 */

/** \file
 * \copydoc cdflib.h
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "cdflib.h"

double
cdflib_algdiv(
    double* a,
    double* b)
{
  static double algdiv;
  static double c;
  static const double c0 =  0.833333333333333e-01;
  static const double c1 = -0.277777777760991e-02;
  static const double c2 =  0.793650666825390e-03;
  static const double c3 = -0.595202931351870e-03;
  static const double c4 =  0.837308034031215e-03;
  static const double c5 = -0.165322962780713e-02;
  static double d;
  static double h;
  static double s11;
  static double s3;
  static double s5;
  static double s7;
  static double s9;
  static double t;
  static double T1;
  static double u;
  static double v;
  static double w;
  static double x;
  static double x2;

  if ( *b <= *a )
  {
    h = *b / *a;
    c = 1.0e0 / ( 1.0e0 + h );
    x = h / ( 1.0e0 + h );
    d = *a + ( *b - 0.5e0 );
  }
  else
  {
    h = *a / *b;
    c = h / ( 1.0e0 + h );
    x = 1.0e0 / ( 1.0e0 + h );
    d = *b + ( *a - 0.5e0 );
  }
//
//  SET SN = (1 - X**N)/(1 - X)
//
  x2 = x * x;
  s3 = 1.0e0 + ( x + x2 );
  s5 = 1.0e0 + ( x + x2 * s3 );
  s7 = 1.0e0 + ( x + x2 * s5 );
  s9 = 1.0e0 + ( x + x2 * s7 );
  s11 = 1.0e0 + ( x + x2 * s9 );
//
//  SET W = DEL(B) - DEL(A + B)
//
  t = pow ( 1.0e0 / *b, 2.0 );

  w = (((( c5 * s11  * t
         + c4 * s9 ) * t
         + c3 * s7 ) * t
         + c2 * s5 ) * t
         + c1 * s3 ) * t
         + c0;

  w *= ( c / *b );
//
//  Combine the results.
//
  T1 = *a / *b;
  u = d * cdflib_alnrel ( &T1 );
  v = *a * ( log ( *b ) - 1.0e0 );

  if ( v < u )
  {
    algdiv = w - v - u;
  }
  else
  {
    algdiv = w - u - v;
  }
  return algdiv;
}

double
cdflib_alnrel(
    double* a)
{
  double alnrel;
  static const double p1 = -0.129418923021993e+01;
  static const double p2 =  0.405303492862024e+00;
  static const double p3 = -0.178874546012214e-01;
  static const double q1 = -0.162752256355323e+01;
  static const double q2 =  0.747811014037616e+00;
  static const double q3 = -0.845104217945565e-01;
  double t;
  double t2;
  double w;
  double x;

  if ( fabs ( *a ) <= 0.375e0 )
  {
    t = *a / ( *a + 2.0e0 );
    t2 = t * t;
    w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)
      / (((q3*t2+q2)*t2+q1)*t2+1.0e0);
    alnrel = 2.0e0 * t * w;
  }
  else
  {
    x = 1.0e0 + *a;
    alnrel = log ( x );
  }
  return alnrel;
}

double
cdflib_apser(
    double* a,
    double* b,
    double* x,
    double* eps)
{
  static const double g = 0.577215664901533e0;
  static double apser,aj,bx,c,j,s,t,tol;

    bx = *b**x;
    t = *x-bx;
    if(*b**eps > 2.e-2) goto S10;
    c = log(*x)+cdflib_psi(b)+g+t;
    goto S20;
S10:
    c = log(bx)+g+t;
S20:
    tol = 5.0e0**eps*fabs(c);
    j = 1.0e0;
    s = 0.0e0;
S30:
    j = j + 1.0e0;
    t = t * (*x-bx/j);
    aj = t/j;
    s = s + aj;
    if(fabs(aj) > tol) goto S30;
    apser = -(*a*(c+s));
    return apser;
}

double
cdflib_bcorr(
    double* a0,
    double* b0)
{
  static const double c0 =  0.833333333333333e-01;
  static const double c1 = -0.277777777760991e-02;
  static const double c2 =  0.793650666825390e-03;
  static const double c3 = -0.595202931351870e-03;
  static const double c4 =  0.837308034031215e-03;
  static const double c5 = -0.165322962780713e-02;
  static double bcorr,a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;

  a = cdflib_fifdmin1 ( *a0, *b0 );
  b = cdflib_fifdmax1 ( *a0, *b0 );
  h = a / b;
  c = h / ( 1.0e0 + h );
  x = 1.0e0 / ( 1.0e0 + h );
  x2 = x * x;
//
//  SET SN = (1 - X**N)/(1 - X)
//
  s3 = 1.0e0 + ( x + x2 );
  s5 = 1.0e0 + ( x + x2 * s3 );
  s7 = 1.0e0 + ( x + x2 * s5 );
  s9 = 1.0e0 + ( x + x2 * s7 );
  s11 = 1.0e0 + ( x + x2 * s9 );
//
//  SET W = DEL(B) - DEL(A + B)
//
  t = pow ( 1.0e0 / b, 2.0 );

  w = (((( c5 * s11  * t + c4
              * s9 ) * t + c3
              * s7 ) * t + c2
              * s5 ) * t + c1
              * s3 ) * t + c0;
  w *= ( c / b );
//
//  COMPUTE  DEL(A) + W
//
  t = pow ( 1.0e0 / a, 2.0 );

  bcorr = ((((( c5 * t + c4 )
                   * t + c3 )
                   * t + c2 )
                   * t + c1 )
                   * t + c0 ) / a + w;
  return bcorr;
}

double
cdflib_beta(
    double a,
    double b)
{
  return ( exp ( cdflib_beta_log ( &a, &b ) ) );
}

double
cdflib_beta_asym(
    double* a,
    double* b,
    double* lambda,
    double* eps)
{
  static const double e0 = 1.12837916709551e0;
  static const double e1 = .353553390593274e0;
  static int num = 20;
//
//  NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
//            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
//            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
//     E0 = 2/SQRT(PI)
//     E1 = 2**(-3/2)
//
  static int K3 = 1;
  static double value;
  static double bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,u,w,w0,z,z0,
    z2,zn,znm1;
  static int i,im1,imj,j,m,mm1,mmj,n,np1;
  static double a0[21],b0[21],c[21],d[21],T1,T2;

    value = 0.0e0;
    if(*a >= *b) goto S10;
    h = *a/ *b;
    r0 = 1.0e0/(1.0e0+h);
    r1 = (*b-*a)/ *b;
    w0 = 1.0e0/sqrt(*a*(1.0e0+h));
    goto S20;
S10:
    h = *b/ *a;
    r0 = 1.0e0/(1.0e0+h);
    r1 = (*b-*a)/ *a;
    w0 = 1.0e0/sqrt(*b*(1.0e0+h));
S20:
    T1 = -(*lambda/ *a);
    T2 = *lambda/ *b;
    f = *a*cdflib_rlog1(&T1)+*b*cdflib_rlog1(&T2);
    t = exp(-f);
    if(t == 0.0e0) return value;
    z0 = sqrt(f);
    z = 0.5e0*(z0/e1);
    z2 = f+f;
    a0[0] = 2.0e0/3.0e0*r1;
    c[0] = -(0.5e0*a0[0]);
    d[0] = -c[0];
    j0 = 0.5e0/e0 * cdflib_error_fc ( &K3, &z0 );
    j1 = e1;
    sum = j0+d[0]*w0*j1;
    s = 1.0e0;
    h2 = h*h;
    hn = 1.0e0;
    w = w0;
    znm1 = z;
    zn = z2;
    for ( n = 2; n <= num; n += 2 )
    {
        hn = h2*hn;
        a0[n-1] = 2.0e0*r0*(1.0e0+h*hn)/((double)n+2.0e0);
        np1 = n+1;
        s += hn;
        a0[np1-1] = 2.0e0*r1*s/((double)n+3.0e0);
        for ( i = n; i <= np1; i++ )
        {
            r = -(0.5e0*((double)i+1.0e0));
            b0[0] = r*a0[0];
            for ( m = 2; m <= i; m++ )
            {
                bsum = 0.0e0;
                mm1 = m-1;
                for ( j = 1; j <= mm1; j++ )
                {
                    mmj = m-j;
                    bsum += (((double)j*r-(double)mmj)*a0[j-1]*b0[mmj-1]);
                }
                b0[m-1] = r*a0[m-1]+bsum/(double)m;
            }
            c[i-1] = b0[i-1]/((double)i+1.0e0);
            dsum = 0.0e0;
            im1 = i-1;
            for ( j = 1; j <= im1; j++ )
            {
                imj = i-j;
                dsum += (d[imj-1]*c[j-1]);
            }
            d[i-1] = -(dsum+c[i-1]);
        }
        j0 = e1*znm1+((double)n-1.0e0)*j0;
        j1 = e1*zn+(double)n*j1;
        znm1 = z2*znm1;
        zn = z2*zn;
        w = w0*w;
        t0 = d[n-1]*w*j0;
        w = w0*w;
        t1 = d[np1-1]*w*j1;
        sum += (t0+t1);
        if(fabs(t0)+fabs(t1) <= *eps*sum) goto S80;
    }
S80:
    u = exp(-cdflib_bcorr(a,b));
    value = e0*t*u*sum;
    return value;
}

double
cdflib_beta_frac(
    double* a,
    double* b,
    double* x,
    double* y,
    double* lambda,
    double* eps)
{
  static double bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;

  bfrac = cdflib_beta_cdflib_rcomp ( a, b, x, y );

  if ( bfrac == 0.0e0 )
  {
    return bfrac;
  }

  c = 1.0e0+*lambda;
  c0 = *b/ *a;
  c1 = 1.0e0+1.0e0/ *a;
  yp1 = *y+1.0e0;
  n = 0.0e0;
  p = 1.0e0;
  s = *a+1.0e0;
  an = 0.0e0;
  bn = anp1 = 1.0e0;
  bnp1 = c/c1;
  r = c1/c;
//
//  CONTINUED FRACTION CALCULATION
//
S10:
  n = n + 1.0e0;
  t = n/ *a;
  w = n*(*b-n)**x;
  e = *a/s;
  alpha = p*(p+c0)*e*e*(w**x);
  e = (1.0e0+t)/(c1+t+t);
  beta = n+w/s+e*(c+n*yp1);
  p = 1.0e0+t;
  s += 2.0e0;
//
//  UPDATE AN, BN, ANP1, AND BNP1
//
  t = alpha*an+beta*anp1;
  an = anp1;
  anp1 = t;
  t = alpha*bn+beta*bnp1;
  bn = bnp1;
  bnp1 = t;
  r0 = r;
  r = anp1/bnp1;

  if ( fabs(r-r0) <= (*eps) * r )
  {
    goto S20;
  }
//
//  RESCALE AN, BN, ANP1, AND BNP1
//
  an /= bnp1;
  bn /= bnp1;
  anp1 = r;
  bnp1 = 1.0e0;
  goto S10;
//
//  TERMINATION
//
S20:
  bfrac = bfrac * r;
  return bfrac;
}

void
cdflib_beta_grat(
    double* a,
    double* b,
    double* x,
    double* y,
    double* w,
    double* eps,
    int* ierr)
{
  static double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
  static int i,n,nm1;
  static double c[30],d[30],T1;

    bm1 = *b-0.5e0-0.5e0;
    nu = *a+0.5e0*bm1;
    if(*y > 0.375e0) goto S10;
    T1 = -*y;
    lnx = cdflib_alnrel(&T1);
    goto S20;
S10:
    lnx = log(*x);
S20:
    z = -(nu*lnx);
    if(*b*z == 0.0e0) goto S70;
//
//  COMPUTATION OF THE EXPANSION
//  SET R = EXP(-Z)*Z**B/GAMMA(B)
//
    r = *b*(1.0e0+cdflib_gam1(b))*exp(*b*log(z));
    r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
    u = cdflib_algdiv(b,a)+*b*log(nu);
    u = r*exp(-u);
    if(u == 0.0e0) goto S70;
    cdflib_gamma_rat1 ( b, &z, &r, &p, &q, eps );
    v = 0.25e0*pow(1.0e0/nu,2.0);
    t2 = 0.25e0*lnx*lnx;
    l = *w/u;
    j = q/r;
    sum = j;
    t = cn = 1.0e0;
    n2 = 0.0e0;
    for ( n = 1; n <= 30; n++ )
    {
        bp2n = *b+n2;
        j = (bp2n*(bp2n+1.0e0)*j+(z+bp2n+1.0e0)*t)*v;
        n2 = n2 + 2.0e0;
        t *= t2;
        cn /= (n2*(n2+1.0e0));
        c[n-1] = cn;
        s = 0.0e0;
        if(n == 1) goto S40;
        nm1 = n-1;
        coef = *b-(double)n;
        for ( i = 1; i <= nm1; i++ )
        {
            s = s + (coef*c[i-1]*d[n-i-1]);
            coef = coef + *b;
        }
S40:
        d[n-1] = bm1*cn+s/(double)n;
        dj = d[n-1]*j;
        sum = sum + dj;
        if(sum <= 0.0e0) goto S70;
        if(fabs(dj) <= *eps*(sum+l)) goto S60;
    }
S60:
//
//  ADD THE RESULTS TO W
//
    *ierr = 0;
    *w = *w + (u*sum);
    return;
S70:
//
//  THE EXPANSION CANNOT BE COMPUTED
//
    *ierr = 1;
    return;
}

void
cdflib_beta_inc(
    double* a,
    double* b,
    double* x,
    double* y,
    double* w,
    double* w1,
    int* ierr)
{
  static int K1 = 1;
  static double a0,b0,eps,lambda,t,x0,y0,z;
  static int ierr1,ind,n;
  static double T2,T3,T4,T5;
//
//  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
//  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
//
    eps = cdflib_dpmpar ( &K1 );
    *w = *w1 = 0.0e0;
    if(*a < 0.0e0 || *b < 0.0e0) goto S270;
    if(*a == 0.0e0 && *b == 0.0e0) goto S280;
    if(*x < 0.0e0 || *x > 1.0e0) goto S290;
    if(*y < 0.0e0 || *y > 1.0e0) goto S300;
    z = *x+*y-0.5e0-0.5e0;
    if(fabs(z) > 3.0e0*eps) goto S310;
    *ierr = 0;
    if(*x == 0.0e0) goto S210;
    if(*y == 0.0e0) goto S230;
    if(*a == 0.0e0) goto S240;
    if(*b == 0.0e0) goto S220;
    eps = cdflib_fifdmax1(eps,1.e-15);
    if(cdflib_fifdmax1(*a,*b) < 1.e-3*eps) goto S260;
    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if(cdflib_fifdmin1(a0,b0) > 1.0e0) goto S40;
//
//  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
//
    if(*x <= 0.5e0) goto S10;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
S10:
    if(b0 < cdflib_fifdmin1(eps,eps*a0)) goto S90;
    if(a0 < cdflib_fifdmin1(eps,eps*b0) && b0*x0 <= 1.0e0) goto S100;
    if(cdflib_fifdmax1(a0,b0) > 1.0e0) goto S20;
    if(a0 >= cdflib_fifdmin1(0.2e0,b0)) goto S110;
    if(pow(x0,a0) <= 0.9e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    n = 20;
    goto S140;
S20:
    if(b0 <= 1.0e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    if(x0 >= 0.1e0) goto S30;
    if(pow(x0*b0,a0) <= 0.7e0) goto S110;
S30:
    if(b0 > 15.0e0) goto S150;
    n = 20;
    goto S140;
S40:
//
//  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
//
    if(*a > *b) goto S50;
    lambda = *a-(*a+*b)**x;
    goto S60;
S50:
    lambda = (*a+*b)**y-*b;
S60:
    if(lambda >= 0.0e0) goto S70;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    lambda = fabs(lambda);
S70:
    if(b0 < 40.0e0 && b0*x0 <= 0.7e0) goto S110;
    if(b0 < 40.0e0) goto S160;
    if(a0 > b0) goto S80;
    if(a0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*a0) goto S130;
    goto S200;
S80:
    if(b0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*b0) goto S130;
    goto S200;
S90:
//
//  EVALUATION OF THE APPROPRIATE ALGORITHM
//
    *w = cdflib_fpser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S100:
    *w1 = cdflib_apser(&a0,&b0,&x0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S110:
    *w = cdflib_beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S120:
    *w1 = cdflib_beta_pser(&b0,&a0,&y0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S130:
    T2 = 15.0e0*eps;
    *w = cdflib_beta_frac ( &a0,&b0,&x0,&y0,&lambda,&T2 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S140:
    *w1 = cdflib_beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    b0 = b0 + (double)n;
S150:
    T3 = 15.0e0*eps;
    cdflib_beta_grat (&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S160:
    n = ( int ) b0;
    b0 -= (double)n;
    if(b0 != 0.0e0) goto S170;
    n -= 1;
    b0 = 1.0e0;
S170:
    *w = cdflib_beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    if(x0 > 0.7e0) goto S180;
    *w = *w + cdflib_beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S180:
    if(a0 > 15.0e0) goto S190;
    n = 20;
    *w = *w + cdflib_beta_up ( &a0, &b0, &x0, &y0, &n, &eps );
    a0 = a0 + (double)n;
S190:
    T4 = 15.0e0*eps;
    cdflib_beta_grat ( &a0, &b0, &x0, &y0, w, &T4, &ierr1 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S200:
    T5 = 100.0e0*eps;
    *w = cdflib_beta_asym ( &a0, &b0, &lambda, &T5 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S210:
//
//  TERMINATION OF THE PROCEDURE
//
    if(*a == 0.0e0) goto S320;
S220:
    *w = 0.0e0;
    *w1 = 1.0e0;
    return;
S230:
    if(*b == 0.0e0) goto S330;
S240:
    *w = 1.0e0;
    *w1 = 0.0e0;
    return;
S250:
    if(ind == 0) return;
    t = *w;
    *w = *w1;
    *w1 = t;
    return;
S260:
//
//  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
//
    *w = *b/(*a+*b);
    *w1 = *a/(*a+*b);
    return;
S270:
//
//  ERROR RETURN
//
    *ierr = 1;
    return;
S280:
    *ierr = 2;
    return;
S290:
    *ierr = 3;
    return;
S300:
    *ierr = 4;
    return;
S310:
    *ierr = 5;
    return;
S320:
    *ierr = 6;
    return;
S330:
    *ierr = 7;
    return;
}

void
cdflib_beta_inc_values(
    int* n_data,
    double* a,
    double* b,
    double* x,
    double* fx)
{
  enum { N_MAX = 30 };

  static const double a_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     1.0E+00,  1.0E+00,  1.0E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.5E+00, 10.0E+00, 10.0E+00,
    10.0E+00, 10.0E+00, 20.0E+00, 20.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00, 30.0E+00,
    30.0E+00, 40.0E+00 };
  static const double b_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  0.5E+00,
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.0E+00,  0.5E+00,  5.0E+00,
     5.0E+00, 10.0E+00,  5.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00, 20.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00 };
  static const double fx_vec[N_MAX] = {
    0.0637686E+00, 0.2048328E+00, 1.0000000E+00, 0.0E+00,
    0.0050126E+00, 0.0513167E+00, 0.2928932E+00, 0.5000000E+00,
    0.028E+00,     0.104E+00,     0.216E+00,     0.352E+00,
    0.500E+00,     0.648E+00,     0.784E+00,     0.896E+00,
    0.972E+00,     0.4361909E+00, 0.1516409E+00, 0.0897827E+00,
    1.0000000E+00, 0.5000000E+00, 0.4598773E+00, 0.2146816E+00,
    0.9507365E+00, 0.5000000E+00, 0.8979414E+00, 0.2241297E+00,
    0.7586405E+00, 0.7001783E+00 };
  static const double x_vec[N_MAX] = {
    0.01E+00, 0.10E+00, 1.00E+00, 0.0E+00,
    0.01E+00, 0.10E+00, 0.50E+00, 0.50E+00,
    0.1E+00,  0.2E+00,  0.3E+00,  0.4E+00,
    0.5E+00,  0.6E+00,  0.7E+00,  0.8E+00,
    0.9E+00,  0.50E+00, 0.90E+00, 0.50E+00,
    1.00E+00, 0.50E+00, 0.80E+00, 0.60E+00,
    0.80E+00, 0.50E+00, 0.60E+00, 0.70E+00,
    0.80E+00, 0.70E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *b = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_beta_log(
    double* a0,
    double* b0)
{
  static const double e = .918938533204673e0;
  static double value,a,b,c,h,u,v,w,z;
  static int i,n;
  static double T1;

    a = cdflib_fifdmin1(*a0,*b0);
    b = cdflib_fifdmax1(*a0,*b0);
    if(a >= 8.0e0) goto S100;
    if(a >= 1.0e0) goto S20;
//
//  PROCEDURE WHEN A .LT. 1
//
    if(b >= 8.0e0) goto S10;
    T1 = a+b;
    value = cdflib_gamma_log ( &a )+( cdflib_gamma_log ( &b )- cdflib_gamma_log ( &T1 ));
    return value;
S10:
    value = cdflib_gamma_log ( &a )+cdflib_algdiv(&a,&b);
    return value;
S20:
//
//  PROCEDURE WHEN 1 .LE. A .LT. 8
//
    if(a > 2.0e0) goto S40;
    if(b > 2.0e0) goto S30;
    value = cdflib_gamma_log ( &a )+ cdflib_gamma_log ( &b )-cdflib_gsumln(&a,&b);
    return value;
S30:
    w = 0.0e0;
    if(b < 8.0e0) goto S60;
    value = cdflib_gamma_log ( &a )+cdflib_algdiv(&a,&b);
    return value;
S40:
//
//  REDUCTION OF A WHEN B .LE. 1000
//
    if(b > 1000.0e0) goto S80;
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        h = a/b;
        w *= (h/(1.0e0+h));
    }
    w = log(w);
    if(b < 8.0e0) goto S60;
    value = w+ cdflib_gamma_log ( &a )+cdflib_algdiv(&a,&b);
    return value;
S60:
//
//  REDUCTION OF B WHEN B .LT. 8
//
    n = ( int ) ( b - 1.0e0 );
    z = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b -= 1.0e0;
        z *= (b/(a+b));
    }
    value = w+log(z)+( cdflib_gamma_log ( &a )+( cdflib_gamma_log ( &b )-cdflib_gsumln(&a,&b)));
    return value;
S80:
//
//  REDUCTION OF A WHEN B .GT. 1000
//
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        w *= (a/(1.0e0+a/b));
    }
    value = log(w)-(double)n*log(b)+( cdflib_gamma_log ( &a )+cdflib_algdiv(&a,&b));
    return value;
S100:
//
//  PROCEDURE WHEN A .GE. 8
//
    w = cdflib_bcorr(&a,&b);
    h = a/b;
    c = h/(1.0e0+h);
    u = -((a-0.5e0)*log(c));
    v = b*cdflib_alnrel(&h);
    if(u <= v) goto S110;
    value = -(0.5e0*log(b))+e+w-v-u;
    return value;
S110:
    value = -(0.5e0*log(b))+e+w-u-v;
    return value;
}

double
cdflib_beta_pser(
    double* a,
    double* b,
    double* x,
    double* eps)
{
  static double bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
  static int i,m;

    bpser = 0.0e0;
    if(*x == 0.0e0) return bpser;
//
//  COMPUTE THE FACTOR X**A/(A*BETA(A,B))
//
    a0 = cdflib_fifdmin1(*a,*b);
    if(a0 < 1.0e0) goto S10;
    z = *a*log(*x)-cdflib_beta_log(a,b);
    bpser = exp(z)/ *a;
    goto S100;
S10:
    b0 = cdflib_fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S90;
    if(b0 > 1.0e0) goto S40;
//
//  PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
//
    bpser = pow(*x,*a);
    if(bpser == 0.0e0) return bpser;
    apb = *a+*b;
    if(apb > 1.0e0) goto S20;
    z = 1.0e0+cdflib_gam1(&apb);
    goto S30;
S20:
    u = *a+*b-1.e0;
    z = (1.0e0+cdflib_gam1(&u))/apb;
S30:
    c = (1.0e0+cdflib_gam1(a))*(1.0e0+cdflib_gam1(b))/z;
    bpser *= (c*(*b/apb));
    goto S100;
S40:
//
//  PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
//
    u = cdflib_gamma_ln1 ( &a0 );
    m = ( int ) ( b0 - 1.0e0 );
    if(m < 1) goto S60;
    c = 1.0e0;
    for ( i = 1; i <= m; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S60:
    z = *a*log(*x)-u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S70;
    t = 1.0e0+cdflib_gam1(&apb);
    goto S80;
S70:
    u = a0+b0-1.e0;
    t = (1.0e0+cdflib_gam1(&u))/apb;
S80:
    bpser = exp(z)*(a0/ *a)*(1.0e0+cdflib_gam1(&b0))/t;
    goto S100;
S90:
//
//  PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
//
    u = cdflib_gamma_ln1 ( &a0 ) + cdflib_algdiv ( &a0, &b0 );
    z = *a*log(*x)-u;
    bpser = a0/ *a*exp(z);
S100:
    if(bpser == 0.0e0 || *a <= 0.1e0**eps) return bpser;
//
//  COMPUTE THE SERIES
//
    sum = n = 0.0e0;
    c = 1.0e0;
    tol = *eps/ *a;
S110:
    n = n + 1.0e0;
    c *= ((0.5e0+(0.5e0-*b/n))**x);
    w = c/(*a+n);
    sum = sum + w;
    if(fabs(w) > tol) goto S110;
    bpser *= (1.0e0+*a*sum);
    return bpser;
}

double
cdflib_beta_cdflib_rcomp(
    double* a,
    double* b,
    double* x,
    double* y)
{
  static const double Const = .398942280401433e0;
  static double brcomp,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  static int i,n;
//
//  CONST = 1/SQRT(2*PI)
//
  static double T1,T2;

    brcomp = 0.0e0;
    if(*x == 0.0e0 || *y == 0.0e0) return brcomp;
    a0 = cdflib_fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = cdflib_alnrel(&T1);
    goto S30;
S10:
    if(*y > 0.375e0) goto S20;
    T2 = -*y;
    lnx = cdflib_alnrel(&T2);
    lny = log(*y);
    goto S30;
S20:
    lnx = log(*x);
    lny = log(*y);
S30:
    z = *a*lnx+*b*lny;
    if(a0 < 1.0e0) goto S40;
    z -= cdflib_beta_log(a,b);
    brcomp = exp(z);
    return brcomp;
S40:
//
//  PROCEDURE FOR A .LT. 1 OR B .LT. 1
//
    b0 = cdflib_fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S120;
    if(b0 > 1.0e0) goto S70;
//
//  ALGORITHM FOR B0 .LE. 1
//
    brcomp = exp(z);
    if(brcomp == 0.0e0) return brcomp;
    apb = *a+*b;
    if(apb > 1.0e0) goto S50;
    z = 1.0e0+cdflib_gam1(&apb);
    goto S60;
S50:
    u = *a+*b-1.e0;
    z = (1.0e0+cdflib_gam1(&u))/apb;
S60:
    c = (1.0e0+cdflib_gam1(a))*(1.0e0+cdflib_gam1(b))/z;
    brcomp = brcomp*(a0*c)/(1.0e0+a0/b0);
    return brcomp;
S70:
//
//  ALGORITHM FOR 1 .LT. B0 .LT. 8
//
    u = cdflib_gamma_ln1 ( &a0 );
    n = ( int ) ( b0 - 1.0e0 );
    if(n < 1) goto S90;
    c = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S90:
    z -= u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S100;
    t = 1.0e0+cdflib_gam1(&apb);
    goto S110;
S100:
    u = a0+b0-1.e0;
    t = (1.0e0+cdflib_gam1(&u))/apb;
S110:
    brcomp = a0*exp(z)*(1.0e0+cdflib_gam1(&b0))/t;
    return brcomp;
S120:
//
//  ALGORITHM FOR B0 .GE. 8
//
    u = cdflib_gamma_ln1 ( &a0 ) + cdflib_algdiv ( &a0, &b0 );
    brcomp = a0*exp(z-u);
    return brcomp;
S130:
//
//  PROCEDURE FOR A .GE. 8 AND B .GE. 8
//
    if(*a > *b) goto S140;
    h = *a/ *b;
    x0 = h/(1.0e0+h);
    y0 = 1.0e0/(1.0e0+h);
    lambda = *a-(*a+*b)**x;
    goto S150;
S140:
    h = *b/ *a;
    x0 = 1.0e0/(1.0e0+h);
    y0 = h/(1.0e0+h);
    lambda = (*a+*b)**y-*b;
S150:
    e = -(lambda/ *a);
    if(fabs(e) > 0.6e0) goto S160;
    u = cdflib_rlog1(&e);
    goto S170;
S160:
    u = e-log(*x/x0);
S170:
    e = lambda/ *b;
    if(fabs(e) > 0.6e0) goto S180;
    v = cdflib_rlog1(&e);
    goto S190;
S180:
    v = e-log(*y/y0);
S190:
    z = exp(-(*a*u+*b*v));
    brcomp = Const*sqrt(*b*x0)*z*exp(-cdflib_bcorr(a,b));
    return brcomp;
}

double
cdflib_beta_rcomp1(
    int* mu,
    double* a,
    double* b,
    double* x,
    double* y)
{
  static const double Const = .398942280401433e0;
  static double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  static int i,n;
//
//     CONST = 1/SQRT(2*PI)
//
  static double T1,T2,T3,T4;

    a0 = cdflib_fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = cdflib_alnrel(&T1);
    goto S30;
S10:
    if(*y > 0.375e0) goto S20;
    T2 = -*y;
    lnx = cdflib_alnrel(&T2);
    lny = log(*y);
    goto S30;
S20:
    lnx = log(*x);
    lny = log(*y);
S30:
    z = *a*lnx+*b*lny;
    if(a0 < 1.0e0) goto S40;
    z -= cdflib_beta_log(a,b);
    brcmp1 = cdflib_esum(mu,&z);
    return brcmp1;
S40:
//
//   PROCEDURE FOR A .LT. 1 OR B .LT. 1
//
    b0 = cdflib_fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S120;
    if(b0 > 1.0e0) goto S70;
//
//  ALGORITHM FOR B0 .LE. 1
//
    brcmp1 = cdflib_esum(mu,&z);
    if(brcmp1 == 0.0e0) return brcmp1;
    apb = *a+*b;
    if(apb > 1.0e0) goto S50;
    z = 1.0e0+cdflib_gam1(&apb);
    goto S60;
S50:
    u = *a+*b-1.e0;
    z = (1.0e0+cdflib_gam1(&u))/apb;
S60:
    c = (1.0e0+cdflib_gam1(a))*(1.0e0+cdflib_gam1(b))/z;
    brcmp1 = brcmp1*(a0*c)/(1.0e0+a0/b0);
    return brcmp1;
S70:
//
//  ALGORITHM FOR 1 .LT. B0 .LT. 8
//
    u = cdflib_gamma_ln1 ( &a0 );
    n = ( int ) ( b0 - 1.0e0 );
    if(n < 1) goto S90;
    c = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S90:
    z -= u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S100;
    t = 1.0e0+cdflib_gam1(&apb);
    goto S110;
S100:
    u = a0+b0-1.e0;
    t = (1.0e0+cdflib_gam1(&u))/apb;
S110:
    brcmp1 = a0*cdflib_esum(mu,&z)*(1.0e0+cdflib_gam1(&b0))/t;
    return brcmp1;
S120:
//
//  ALGORITHM FOR B0 .GE. 8
//
    u = cdflib_gamma_ln1 ( &a0 ) + cdflib_algdiv ( &a0, &b0 );
    T3 = z-u;
    brcmp1 = a0*cdflib_esum(mu,&T3);
    return brcmp1;
S130:
//
//    PROCEDURE FOR A .GE. 8 AND B .GE. 8
//
    if(*a > *b) goto S140;
    h = *a/ *b;
    x0 = h/(1.0e0+h);
    y0 = 1.0e0/(1.0e0+h);
    lambda = *a-(*a+*b)**x;
    goto S150;
S140:
    h = *b/ *a;
    x0 = 1.0e0/(1.0e0+h);
    y0 = h/(1.0e0+h);
    lambda = (*a+*b)**y-*b;
S150:
    e = -(lambda/ *a);
    if(fabs(e) > 0.6e0) goto S160;
    u = cdflib_rlog1(&e);
    goto S170;
S160:
    u = e-log(*x/x0);
S170:
    e = lambda/ *b;
    if(fabs(e) > 0.6e0) goto S180;
    v = cdflib_rlog1(&e);
    goto S190;
S180:
    v = e-log(*y/y0);
S190:
    T4 = -(*a*u+*b*v);
    z = cdflib_esum(mu,&T4);
    brcmp1 = Const*sqrt(*b*x0)*z*exp(-cdflib_bcorr(a,b));
    return brcmp1;
}

double
cdflib_beta_up(
    double* a,
    double* b,
    double* x,
    double* y,
    int* n,
    double* eps)
{
  static int K1 = 1;
  static int K2 = 0;
  static double bup,ap1,apb,d,l,r,t,w;
  static int i,k,kp1,mu,nm1;
//
//  OBTAIN THE SCALING FACTOR EXP(-MU) AND
//  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
//
    apb = *a+*b;
    ap1 = *a+1.0e0;
    mu = 0;
    d = 1.0e0;
    if(*n == 1 || *a < 1.0e0) goto S10;
    if(apb < 1.1e0*ap1) goto S10;
    mu = ( int ) fabs ( cdflib_exparg(&K1) );
    k = ( int ) cdflib_exparg ( &K2 );
    if(k < mu) mu = k;
    t = mu;
    d = exp(-t);
S10:
    bup = cdflib_beta_rcomp1 ( &mu, a, b, x, y ) / *a;
    if(*n == 1 || bup == 0.0e0) return bup;
    nm1 = *n-1;
    w = d;
//
//  LET K BE THE INDEX OF THE MAXIMUM TERM
//
    k = 0;
    if(*b <= 1.0e0) goto S50;
    if(*y > 1.e-4) goto S20;
    k = nm1;
    goto S30;
S20:
    r = (*b-1.0e0)**x/ *y-*a;
    if(r < 1.0e0) goto S50;
    t = ( double ) nm1;
    k = nm1;
    if ( r < t ) k = ( int ) r;
S30:
//
//          ADD THE INCREASING TERMS OF THE SERIES
//
    for ( i = 1; i <= k; i++ )
    {
        l = i-1;
        d = (apb+l)/(ap1+l)**x*d;
        w = w + d;
    }
    if(k == nm1) goto S70;
S50:
//
//          ADD THE REMAINING TERMS OF THE SERIES
//
    kp1 = k+1;
    for ( i = kp1; i <= nm1; i++ )
    {
        l = i-1;
        d = (apb+l)/(ap1+l)**x*d;
        w = w + d;
        if(d <= *eps*w) goto S70;
    }
S70:
//
//  TERMINATE THE PROCEDURE
//
    bup *= w;
    return bup;
}

void
cdflib_binomial_cdf_values(
    int* n_data,
    int* a,
    double* b,
    int* x,
    double* fx)
{
  enum { N_MAX = 17 };

  static const int a_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  4,  4,  4,
     4, 10, 10, 10,
    10, 10, 10, 10,
    10 };
  static const double b_vec[N_MAX] = {
    0.05E+00, 0.05E+00, 0.05E+00, 0.50E+00,
    0.50E+00, 0.25E+00, 0.25E+00, 0.25E+00,
    0.25E+00, 0.05E+00, 0.10E+00, 0.15E+00,
    0.20E+00, 0.25E+00, 0.30E+00, 0.40E+00,
    0.50E+00 };
  static const double fx_vec[N_MAX] = {
    0.9025E+00, 0.9975E+00, 1.0000E+00, 0.2500E+00,
    0.7500E+00, 0.3164E+00, 0.7383E+00, 0.9492E+00,
    0.9961E+00, 0.9999E+00, 0.9984E+00, 0.9901E+00,
    0.9672E+00, 0.9219E+00, 0.8497E+00, 0.6331E+00,
    0.3770E+00 };
  static const int x_vec[N_MAX] = {
     0, 1, 2, 0,
     1, 0, 1, 2,
     3, 4, 4, 4,
     4, 4, 4, 4,
     4 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0.0E+00;
    *x = 0;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

void
cdflib_cdfbet(
    int* which,
    double* p,
    double* q,
    double* x,
    double* y,
    double* a,
    double* b,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K3 = 1.0e0;
  static double K8 = 0.5e0;
  static double K9 = 5.0e0;
  static double fx,xhi,xlo,cum,ccum,xy,pq;
  static unsigned long qhi,qleft,qporq;
  static double T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q < 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S150;
//
//     X
//
    if(!(*x < 0.0e0 || *x > 1.0e0)) goto S140;
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = -4;
    return;
S150:
S140:
    if(*which == 2) goto S190;
//
//     Y
//
    if(!(*y < 0.0e0 || *y > 1.0e0)) goto S180;
    if(!(*y < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = -5;
    return;
S190:
S180:
    if(*which == 3) goto S210;
//
//     A
//
    if(!(*a <= 0.0e0)) goto S200;
    *bound = 0.0e0;
    *status = -6;
    return;
S210:
S200:
    if(*which == 4) goto S230;
//
//     B
//
    if(!(*b <= 0.0e0)) goto S220;
    *bound = 0.0e0;
    *status = -7;
    return;
S230:
S220:
    if(*which == 1) goto S270;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S260;
    if(!(pq < 0.0e0)) goto S240;
    *bound = 0.0e0;
    goto S250;
S240:
    *bound = 1.0e0;
S250:
    *status = 3;
    return;
S270:
S260:
    if(*which == 2) goto S310;
//
//     X + Y
//
    xy = *x+*y;
    if(!(fabs(xy-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S300;
    if(!(xy < 0.0e0)) goto S280;
    *bound = 0.0e0;
    goto S290;
S280:
    *bound = 1.0e0;
S290:
    *status = 4;
    return;
S310:
S300:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        cdflib_cumbet(x,y,a,b,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating X and Y
//
        T4 = atol;
        T5 = tol;
        cdflib_dstzr(&K2,&K3,&T4,&T5);
        if(!qporq) goto S340;
        *status = 0;
        cdflib_dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
        *y = one-*x;
S320:
        if(!(*status == 1)) goto S330;
        cdflib_cumbet(x,y,a,b,&cum,&ccum);
        fx = cum-*p;
        cdflib_dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
        *y = one-*x;
        goto S320;
S330:
        goto S370;
S340:
        *status = 0;
        cdflib_dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
        *x = one-*y;
S350:
        if(!(*status == 1)) goto S360;
        cdflib_cumbet(x,y,a,b,&cum,&ccum);
        fx = ccum-*q;
        cdflib_dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
        *x = one-*y;
        goto S350;
S370:
S360:
        if(!(*status == -1)) goto S400;
        if(!qleft) goto S380;
        *status = 1;
        *bound = 0.0e0;
        goto S390;
S380:
        *status = 2;
        *bound = 1.0e0;
S400:
S390:
        ;
    }
    else if(3 == *which) {
//
//     Computing A
//
        *a = 5.0e0;
        T6 = zero;
        T7 = inf;
        T10 = atol;
        T11 = tol;
        cdflib_dstinv(&T6,&T7,&K8,&K8,&K9,&T10,&T11);
        *status = 0;
        cdflib_dinvr(status,a,&fx,&qleft,&qhi);
S410:
        if(!(*status == 1)) goto S440;
        cdflib_cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S420;
        fx = cum-*p;
        goto S430;
S420:
        fx = ccum-*q;
S430:
        cdflib_dinvr(status,a,&fx,&qleft,&qhi);
        goto S410;
S440:
        if(!(*status == -1)) goto S470;
        if(!qleft) goto S450;
        *status = 1;
        *bound = zero;
        goto S460;
S450:
        *status = 2;
        *bound = inf;
S470:
S460:
        ;
    }
    else if(4 == *which) {
//
//     Computing B
//
        *b = 5.0e0;
        T12 = zero;
        T13 = inf;
        T14 = atol;
        T15 = tol;
        cdflib_dstinv(&T12,&T13,&K8,&K8,&K9,&T14,&T15);
        *status = 0;
        cdflib_dinvr(status,b,&fx,&qleft,&qhi);
S480:
        if(!(*status == 1)) goto S510;
        cdflib_cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S490;
        fx = cum-*p;
        goto S500;
S490:
        fx = ccum-*q;
S500:
        cdflib_dinvr(status,b,&fx,&qleft,&qhi);
        goto S480;
S510:
        if(!(*status == -1)) goto S540;
        if(!qleft) goto S520;
        *status = 1;
        *bound = zero;
        goto S530;
S520:
        *status = 2;
        *bound = inf;
S530:
        ;
    }
S540:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
# undef one
}

void
cdflib_cdfbin(
    int* which,
    double* p,
    double* q,
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    int* status,
    double* bound)
{
# define atol (1.0e-50)
# define tol (1.0e-8)
# define zero (1.0e-300)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double K11 = 1.0e0;
  static double fx,xhi,xlo,cum,ccum,pq,prompr;
  static unsigned long qhi,qleft,qporq;
  static double T5,T6,T7,T8,T9,T10,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 && *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q < 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 3) goto S130;
//
//     XN
//
    if(!(*xn <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -5;
    return;
S130:
S120:
    if(*which == 2) goto S170;
//
//     S
//
    if(!(*s < 0.0e0 || (*which != 3 && *s > *xn))) goto S160;
    if(!(*s < 0.0e0)) goto S140;
    *bound = 0.0e0;
    goto S150;
S140:
    *bound = *xn;
S150:
    *status = -4;
    return;
S170:
S160:
    if(*which == 4) goto S210;
//
//     PR
//
    if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S200;
    if(!(*pr < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = -6;
    return;
S210:
S200:
    if(*which == 4) goto S250;
//
//     OMPR
//
    if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S240;
    if(!(*ompr < 0.0e0)) goto S220;
    *bound = 0.0e0;
    goto S230;
S220:
    *bound = 1.0e0;
S230:
    *status = -7;
    return;
S250:
S240:
    if(*which == 1) goto S290;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S280;
    if(!(pq < 0.0e0)) goto S260;
    *bound = 0.0e0;
    goto S270;
S260:
    *bound = 1.0e0;
S270:
    *status = 3;
    return;
S290:
S280:
    if(*which == 4) goto S330;
//
//     PR + OMPR
//
    prompr = *pr+*ompr;
    if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S320;
    if(!(prompr < 0.0e0)) goto S300;
    *bound = 0.0e0;
    goto S310;
S300:
    *bound = 1.0e0;
S310:
    *status = 4;
    return;
S330:
S320:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cdflib_cumbin(s,xn,pr,ompr,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T5 = atol;
        T6 = tol;
        cdflib_dstinv(&K2,xn,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
S340:
        if(!(*status == 1)) goto S370;
        cdflib_cumbin(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S350;
        fx = cum-*p;
        goto S360;
S350:
        fx = ccum-*q;
S360:
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
        goto S340;
S370:
        if(!(*status == -1)) goto S400;
        if(!qleft) goto S380;
        *status = 1;
        *bound = 0.0e0;
        goto S390;
S380:
        *status = 2;
        *bound = *xn;
S400:
S390:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XN
//
        *xn = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        cdflib_dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        cdflib_dinvr(status,xn,&fx,&qleft,&qhi);
S410:
        if(!(*status == 1)) goto S440;
        cdflib_cumbin(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S420;
        fx = cum-*p;
        goto S430;
S420:
        fx = ccum-*q;
S430:
        cdflib_dinvr(status,xn,&fx,&qleft,&qhi);
        goto S410;
S440:
        if(!(*status == -1)) goto S470;
        if(!qleft) goto S450;
        *status = 1;
        *bound = zero;
        goto S460;
S450:
        *status = 2;
        *bound = inf;
S470:
S460:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PR and OMPR
//
        T12 = atol;
        T13 = tol;
        cdflib_dstzr(&K2,&K11,&T12,&T13);
        if(!qporq) goto S500;
        *status = 0;
        cdflib_dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
S480:
        if(!(*status == 1)) goto S490;
        cdflib_cumbin(s,xn,pr,ompr,&cum,&ccum);
        fx = cum-*p;
        cdflib_dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
        goto S480;
S490:
        goto S530;
S500:
        *status = 0;
        cdflib_dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
S510:
        if(!(*status == 1)) goto S520;
        cdflib_cumbin(s,xn,pr,ompr,&cum,&ccum);
        fx = ccum-*q;
        cdflib_dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
        goto S510;
S530:
S520:
        if(!(*status == -1)) goto S560;
        if(!qleft) goto S540;
        *status = 1;
        *bound = 0.0e0;
        goto S550;
S540:
        *status = 2;
        *bound = 1.0e0;
S550:
        ;
    }
S560:
    return;
# undef atol
# undef tol
# undef zero
# undef inf
# undef one
}

void
cdflib_cdfchi(
    int* which,
    double* p,
    double* q,
    double* x,
    double* df,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq,porq;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T11;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     X
//
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 1) goto S190;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S180;
    if(!(pq < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = 3;
    return;
S190:
S180:
    if(*which == 1) goto S220;
//
//     Select the minimum of P or Q
//
    qporq = *p <= *q;
    if(!qporq) goto S200;
    porq = *p;
    goto S210;
S200:
    porq = *q;
S220:
S210:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        *status = 0;
        cdflib_cumchi(x,df,p,q);
        if(porq > 1.5e0) {
            *status = 10;
            return;
        }
    }
    else if(2 == *which) {
//
//     Calculating X
//
        *x = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        cdflib_dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        cdflib_dinvr(status,x,&fx,&qleft,&qhi);
S230:
        if(!(*status == 1)) goto S270;
        cdflib_cumchi(x,df,&cum,&ccum);
        if(!qporq) goto S240;
        fx = cum-*p;
        goto S250;
S240:
        fx = ccum-*q;
S250:
        if(!(fx+porq > 1.5e0)) goto S260;
        *status = 10;
        return;
S260:
        cdflib_dinvr(status,x,&fx,&qleft,&qhi);
        goto S230;
S270:
        if(!(*status == -1)) goto S300;
        if(!qleft) goto S280;
        *status = 1;
        *bound = 0.0e0;
        goto S290;
S280:
        *status = 2;
        *bound = inf;
S300:
S290:
        ;
    }
    else if(3 == *which) {
//
//  Calculating DF
//
        *df = 5.0e0;
        T8 = zero;
        T9 = inf;
        T10 = atol;
        T11 = tol;
        cdflib_dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
        *status = 0;
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
S310:
        if(!(*status == 1)) goto S350;
        cdflib_cumchi(x,df,&cum,&ccum);
        if(!qporq) goto S320;
        fx = cum-*p;
        goto S330;
S320:
        fx = ccum-*q;
S330:
        if(!(fx+porq > 1.5e0)) goto S340;
        *status = 10;
        return;
S340:
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
        goto S310;
S350:
        if(!(*status == -1)) goto S380;
        if(!qleft) goto S360;
        *status = 1;
        *bound = zero;
        goto S370;
S360:
        *status = 2;
        *bound = inf;
S370:
        ;
    }
S380:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
}

void
cdflib_cdfchn(
    int* which,
    double* p,
    double* q,
    double* x,
    double* df,
    double* pnonc,
    int* status,
    double* bound)
{
# define tent4 1.0e4
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define one (1.0e0-1.0e-16)
# define inf 1.0e300

  static double K1 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double fx,cum,ccum;
  static unsigned long qhi,qleft;
  static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 2) goto S90;
//
//     X
//
    if(!(*x < 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -4;
    return;
S90:
S80:
    if(*which == 3) goto S110;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S100;
    *bound = 0.0e0;
    *status = -5;
    return;
S110:
S100:
    if(*which == 4) goto S130;
//
//     PNONC
//
    if(!(*pnonc < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -6;
    return;
S130:
S120:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        cdflib_cumchn(x,df,pnonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating X
//
        *x = 5.0e0;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        cdflib_dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        cdflib_dinvr(status,x,&fx,&qleft,&qhi);
S140:
        if(!(*status == 1)) goto S150;
        cdflib_cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,x,&fx,&qleft,&qhi);
        goto S140;
S150:
        if(!(*status == -1)) goto S180;
        if(!qleft) goto S160;
        *status = 1;
        *bound = 0.0e0;
        goto S170;
S160:
        *status = 2;
        *bound = inf;
S180:
S170:
        ;
    }
    else if(3 == *which) {
//
//     Calculating DF
//
        *df = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        cdflib_dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
S190:
        if(!(*status == 1)) goto S200;
        cdflib_cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
        goto S190;
S200:
        if(!(*status == -1)) goto S230;
        if(!qleft) goto S210;
        *status = 1;
        *bound = zero;
        goto S220;
S210:
        *status = 2;
        *bound = inf;
S230:
S220:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PNONC
//
        *pnonc = 5.0e0;
        T11 = tent4;
        T12 = atol;
        T13 = tol;
        cdflib_dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
        *status = 0;
        cdflib_dinvr(status,pnonc,&fx,&qleft,&qhi);
S240:
        if(!(*status == 1)) goto S250;
        cdflib_cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,pnonc,&fx,&qleft,&qhi);
        goto S240;
S250:
        if(!(*status == -1)) goto S280;
        if(!qleft) goto S260;
        *status = 1;
        *bound = zero;
        goto S270;
S260:
        *status = 2;
        *bound = tent4;
S270:
        ;
    }
S280:
    return;
# undef tent4
# undef tol
# undef atol
# undef zero
# undef one
# undef inf
}

void
cdflib_cdff(
    int* which,
    double* p,
    double* q,
    double* f,
    double* dfn,
    double* dfd,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double pq,fx,cum,ccum;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;

  *status = 0;
  *bound = 0.0;
//
//  Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     F
//
    if(!(*f < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     DFN
//
    if(!(*dfn <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     DFD
//
    if(!(*dfd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
    if(*which == 1) goto S210;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * cdflib_dpmpar ( &K1 ) ) ) goto S200;
    if(!(pq < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = 3;
    return;
S210:
S200:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cdflib_cumf(f,dfn,dfd,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating F
//
        *f = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        cdflib_dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        cdflib_dinvr(status,f,&fx,&qleft,&qhi);
S220:
        if(!(*status == 1)) goto S250;
        cdflib_cumf(f,dfn,dfd,&cum,&ccum);
        if(!qporq) goto S230;
        fx = cum-*p;
        goto S240;
S230:
        fx = ccum-*q;
S240:
        cdflib_dinvr(status,f,&fx,&qleft,&qhi);
        goto S220;
S250:
        if(!(*status == -1)) goto S280;
        if(!qleft) goto S260;
        *status = 1;
        *bound = 0.0e0;
        goto S270;
S260:
        *status = 2;
        *bound = inf;
S280:
S270:
        ;
    }
//
//  Calculate DFN.
//
//  Note that, in the original calculation, the lower bound for DFN was 0.
//  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
//  The lower bound was set to the more reasonable value of 1.
//  JVB, 14 April 2007.
//
  else if ( 3 == *which )
  {

    T8 = 1.0;
    T9 = inf;
    T10 = atol;
    T11 = tol;
    cdflib_dstinv ( &T8, &T9, &K4, &K4, &K5, &T10, &T11 );

    *status = 0;
    *dfn = 5.0;
    fx = 0.0;

    cdflib_dinvr ( status, dfn, &fx, &qleft, &qhi );

    while ( *status == 1 )
    {
      cdflib_cumf ( f, dfn, dfd, &cum, &ccum );

      if ( *p <= *q )
      {
        fx = cum - *p;
      }
      else
      {
        fx = ccum - *q;
      }
      cdflib_dinvr ( status, dfn, &fx, &qleft, &qhi );
    }

    if ( *status == -1 )
    {
      if ( qleft )
      {
        *status = 1;
        *bound = 1.0;
      }
      else
      {
        *status = 2;
        *bound = inf;
      }
    }
  }
//
//  Calculate DFD.
//
//  Note that, in the original calculation, the lower bound for DFD was 0.
//  Using DFD = 0 causes an error in CUMF when it calls BETA_INC.
//  The lower bound was set to the more reasonable value of 1.
//  JVB, 14 April 2007.
//
//
  else if ( 4 == *which )
  {

    T12 = 1.0;
    T13 = inf;
    T14 = atol;
    T15 = tol;
    cdflib_dstinv ( &T12, &T13, &K4, &K4, &K5, &T14, &T15 );

    *status = 0;
    *dfd = 5.0;
    fx = 0.0;
    cdflib_dinvr ( status, dfd, &fx, &qleft, &qhi );

    while ( *status == 1 )
    {
      cdflib_cumf ( f, dfn, dfd, &cum, &ccum );

      if ( *p <= *q )
      {
        fx = cum - *p;
      }
      else
      {
        fx = ccum - *q;
      }
      cdflib_dinvr ( status, dfd, &fx, &qleft, &qhi );
    }

    if ( *status == -1 )
    {
      if ( qleft )
      {
        *status = 1;
        *bound = 1.0;
      }
      else
      {
        *status = 2;
        *bound = inf;
      }
    }
  }

  return;
# undef tol
# undef atol
# undef zero
# undef inf
}

void
cdflib_cdffnc(
    int* which,
    double* p,
    double* q,
    double* f,
    double* dfn,
    double* dfd,
    double* phonc,
    int* status,
    double* bound)
{
# define tent4 1.0e4
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define one (1.0e0-1.0e-16)
# define inf 1.0e300

  static double K1 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double fx,cum,ccum;
  static unsigned long qhi,qleft;
  static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 5)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 5.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 2) goto S90;
//
//     F
//
    if(!(*f < 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -4;
    return;
S90:
S80:
    if(*which == 3) goto S110;
//
//     DFN
//
    if(!(*dfn <= 0.0e0)) goto S100;
    *bound = 0.0e0;
    *status = -5;
    return;
S110:
S100:
    if(*which == 4) goto S130;
//
//     DFD
//
    if(!(*dfd <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -6;
    return;
S130:
S120:
    if(*which == 5) goto S150;
//
//     PHONC
//
    if(!(*phonc < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -7;
    return;
S150:
S140:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cdflib_cumfnc(f,dfn,dfd,phonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating F
//
        *f = 5.0e0;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        cdflib_dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        cdflib_dinvr(status,f,&fx,&qleft,&qhi);
S160:
        if(!(*status == 1)) goto S170;
        cdflib_cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,f,&fx,&qleft,&qhi);
        goto S160;
S170:
        if(!(*status == -1)) goto S200;
        if(!qleft) goto S180;
        *status = 1;
        *bound = 0.0e0;
        goto S190;
S180:
        *status = 2;
        *bound = inf;
S200:
S190:
        ;
    }
    else if(3 == *which) {
//
//     Calculating DFN
//
        *dfn = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        cdflib_dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        cdflib_dinvr(status,dfn,&fx,&qleft,&qhi);
S210:
        if(!(*status == 1)) goto S220;
        cdflib_cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,dfn,&fx,&qleft,&qhi);
        goto S210;
S220:
        if(!(*status == -1)) goto S250;
        if(!qleft) goto S230;
        *status = 1;
        *bound = zero;
        goto S240;
S230:
        *status = 2;
        *bound = inf;
S250:
S240:
        ;
    }
    else if(4 == *which) {
//
//     Calculating DFD
//
        *dfd = 5.0e0;
        T11 = zero;
        T12 = inf;
        T13 = atol;
        T14 = tol;
        cdflib_dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
        *status = 0;
        cdflib_dinvr(status,dfd,&fx,&qleft,&qhi);
S260:
        if(!(*status == 1)) goto S270;
        cdflib_cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,dfd,&fx,&qleft,&qhi);
        goto S260;
S270:
        if(!(*status == -1)) goto S300;
        if(!qleft) goto S280;
        *status = 1;
        *bound = zero;
        goto S290;
S280:
        *status = 2;
        *bound = inf;
S300:
S290:
        ;
    }
    else if(5 == *which) {
//
//     Calculating PHONC
//
        *phonc = 5.0e0;
        T15 = tent4;
        T16 = atol;
        T17 = tol;
        cdflib_dstinv(&K1,&T15,&K3,&K3,&K4,&T16,&T17);
        *status = 0;
        cdflib_dinvr(status,phonc,&fx,&qleft,&qhi);
S310:
        if(!(*status == 1)) goto S320;
        cdflib_cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        cdflib_dinvr(status,phonc,&fx,&qleft,&qhi);
        goto S310;
S320:
        if(!(*status == -1)) goto S350;
        if(!qleft) goto S330;
        *status = 1;
        *bound = 0.0e0;
        goto S340;
S330:
        *status = 2;
        *bound = tent4;
S340:
        ;
    }
S350:
    return;
# undef tent4
# undef tol
# undef atol
# undef zero
# undef one
# undef inf
}

void
cdflib_cdfgam(
    int* which,
    double* p,
    double* q,
    double* x,
    double* shape,
    double* scale,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K5 = 0.5e0;
  static double K6 = 5.0e0;
  static double xx,fx,xscale,cum,ccum,pq,porq;
  static int ierr;
  static unsigned long qhi,qleft,qporq;
  static double T2,T3,T4,T7,T8,T9;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     X
//
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     SHAPE
//
    if(!(*shape <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     SCALE
//
    if(!(*scale <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
    if(*which == 1) goto S210;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S200;
    if(!(pq < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = 3;
    return;
S210:
S200:
    if(*which == 1) goto S240;
//
//     Select the minimum of P or Q
//
    qporq = *p <= *q;
    if(!qporq) goto S220;
    porq = *p;
    goto S230;
S220:
    porq = *q;
S240:
S230:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        *status = 0;
        xscale = *x**scale;
        cdflib_cumgam(&xscale,shape,p,q);
        if(porq > 1.5e0) *status = 10;
    }
    else if(2 == *which) {
//
//     Computing X
//
        T2 = -1.0e0;
        cdflib_gamma_inc_inv ( shape, &xx, &T2, p, q, &ierr );
        if(ierr < 0.0e0) {
            *status = 10;
            return;
        }
        else  {
            *x = xx/ *scale;
            *status = 0;
        }
    }
    else if(3 == *which) {
//
//     Computing SHAPE
//
        *shape = 5.0e0;
        xscale = *x**scale;
        T3 = zero;
        T4 = inf;
        T7 = atol;
        T8 = tol;
        cdflib_dstinv(&T3,&T4,&K5,&K5,&K6,&T7,&T8);
        *status = 0;
        cdflib_dinvr(status,shape,&fx,&qleft,&qhi);
S250:
        if(!(*status == 1)) goto S290;
        cdflib_cumgam(&xscale,shape,&cum,&ccum);
        if(!qporq) goto S260;
        fx = cum-*p;
        goto S270;
S260:
        fx = ccum-*q;
S270:
        if(!((qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))) goto S280;
        *status = 10;
        return;
S280:
        cdflib_dinvr(status,shape,&fx,&qleft,&qhi);
        goto S250;
S290:
        if(!(*status == -1)) goto S320;
        if(!qleft) goto S300;
        *status = 1;
        *bound = zero;
        goto S310;
S300:
        *status = 2;
        *bound = inf;
S320:
S310:
        ;
    }
    else if(4 == *which) {
//
//     Computing SCALE
//
        T9 = -1.0e0;
        cdflib_gamma_inc_inv ( shape, &xx, &T9, p, q, &ierr );
        if(ierr < 0.0e0) {
            *status = 10;
            return;
        }
        else  {
            *scale = xx/ *x;
            *status = 0;
        }
    }
    return;
# undef tol
# undef atol
# undef zero
# undef inf
}

void
cdflib_cdfnbn(
    int* which,
    double* p,
    double* q,
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double K11 = 1.0e0;
  static double fx,xhi,xlo,pq,prompr,cum,ccum;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     S
//
    if(!(*s < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     XN
//
    if(!(*xn < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S190;
//
//     PR
//
    if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S180;
    if(!(*pr < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = -6;
    return;
S190:
S180:
    if(*which == 4) goto S230;
//
//     OMPR
//
    if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S220;
    if(!(*ompr < 0.0e0)) goto S200;
    *bound = 0.0e0;
    goto S210;
S200:
    *bound = 1.0e0;
S210:
    *status = -7;
    return;
S230:
S220:
    if(*which == 1) goto S270;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S260;
    if(!(pq < 0.0e0)) goto S240;
    *bound = 0.0e0;
    goto S250;
S240:
    *bound = 1.0e0;
S250:
    *status = 3;
    return;
S270:
S260:
    if(*which == 4) goto S310;
//
//     PR + OMPR
//
    prompr = *pr+*ompr;
    if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S300;
    if(!(prompr < 0.0e0)) goto S280;
    *bound = 0.0e0;
    goto S290;
S280:
    *bound = 1.0e0;
S290:
    *status = 4;
    return;
S310:
S300:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cdflib_cumnbn(s,xn,pr,ompr,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        cdflib_dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
S320:
        if(!(*status == 1)) goto S350;
        cdflib_cumnbn(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S330;
        fx = cum-*p;
        goto S340;
S330:
        fx = ccum-*q;
S340:
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
        goto S320;
S350:
        if(!(*status == -1)) goto S380;
        if(!qleft) goto S360;
        *status = 1;
        *bound = 0.0e0;
        goto S370;
S360:
        *status = 2;
        *bound = inf;
S380:
S370:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XN
//
        *xn = 5.0e0;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        cdflib_dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
        *status = 0;
        cdflib_dinvr(status,xn,&fx,&qleft,&qhi);
S390:
        if(!(*status == 1)) goto S420;
        cdflib_cumnbn(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S400;
        fx = cum-*p;
        goto S410;
S400:
        fx = ccum-*q;
S410:
        cdflib_dinvr(status,xn,&fx,&qleft,&qhi);
        goto S390;
S420:
        if(!(*status == -1)) goto S450;
        if(!qleft) goto S430;
        *status = 1;
        *bound = 0.0e0;
        goto S440;
S430:
        *status = 2;
        *bound = inf;
S450:
S440:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PR and OMPR
//
        T12 = atol;
        T13 = tol;
        cdflib_dstzr(&K2,&K11,&T12,&T13);
        if(!qporq) goto S480;
        *status = 0;
        cdflib_dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
S460:
        if(!(*status == 1)) goto S470;
        cdflib_cumnbn(s,xn,pr,ompr,&cum,&ccum);
        fx = cum-*p;
        cdflib_dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
        goto S460;
S470:
        goto S510;
S480:
        *status = 0;
        cdflib_dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
S490:
        if(!(*status == 1)) goto S500;
        cdflib_cumnbn(s,xn,pr,ompr,&cum,&ccum);
        fx = ccum-*q;
        cdflib_dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
        goto S490;
S510:
S500:
        if(!(*status == -1)) goto S540;
        if(!qleft) goto S520;
        *status = 1;
        *bound = 0.0e0;
        goto S530;
S520:
        *status = 2;
        *bound = 1.0e0;
S530:
        ;
    }
S540:
    return;
# undef tol
# undef atol
# undef inf
# undef one
}

void
cdflib_cdfnor(
    int* which,
    double* p,
    double* q,
    double* x,
    double* mean,
    double* sd,
    int* status,
    double* bound)
{
  static int K1 = 1;
  static double z,pq;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    *status = 0;
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 1) goto S150;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S140;
    if(!(pq < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = 3;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     SD
//
    if(!(*sd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Computing P
//
        z = (*x-*mean)/ *sd;
        cdflib_cumnor(&z,p,q);
    }
    else if(2 == *which) {
//
//     Computing X
//
        z = cdflib_dinvnr(p,q);
        *x = *sd*z+*mean;
    }
    else if(3 == *which) {
//
//     Computing the MEAN
//
        z = cdflib_dinvnr(p,q);
        *mean = *x-*sd*z;
    }
    else if(4 == *which) {
//
//     Computing SD
//
        z = cdflib_dinvnr(p,q);
        *sd = (*x-*mean)/z;
    }
    return;
}

void
cdflib_cdfpoi(
    int* which,
    double* p,
    double* q,
    double* s,
    double* xlam,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     S
//
    if(!(*s < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     XLAM
//
    if(!(*xlam < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 1) goto S190;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S180;
    if(!(pq < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = 3;
    return;
S190:
S180:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cdflib_cumpoi(s,xlam,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        cdflib_dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
S200:
        if(!(*status == 1)) goto S230;
        cdflib_cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S210;
        fx = cum-*p;
        goto S220;
S210:
        fx = ccum-*q;
S220:
        cdflib_dinvr(status,s,&fx,&qleft,&qhi);
        goto S200;
S230:
        if(!(*status == -1)) goto S260;
        if(!qleft) goto S240;
        *status = 1;
        *bound = 0.0e0;
        goto S250;
S240:
        *status = 2;
        *bound = inf;
S260:
S250:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XLAM
//
        *xlam = 5.0e0;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        cdflib_dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
        *status = 0;
        cdflib_dinvr(status,xlam,&fx,&qleft,&qhi);
S270:
        if(!(*status == 1)) goto S300;
        cdflib_cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S280;
        fx = cum-*p;
        goto S290;
S280:
        fx = ccum-*q;
S290:
        cdflib_dinvr(status,xlam,&fx,&qleft,&qhi);
        goto S270;
S300:
        if(!(*status == -1)) goto S330;
        if(!qleft) goto S310;
        *status = 1;
        *bound = 0.0e0;
        goto S320;
S310:
        *status = 2;
        *bound = inf;
S320:
        ;
    }
S330:
    return;
# undef tol
# undef atol
# undef inf
}

void
cdflib_cdft(
    int* which,
    double* p,
    double* q,
    double* t,
    double* df,
    int* status,
    double* bound)
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e30
# define maxdf 1.0e10

  static int K1 = 1;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq;
  static unsigned long qhi,qleft,qporq;
  static double T2,T3,T6,T7,T8,T9,T10,T11;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 3) goto S130;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -5;
    return;
S130:
S120:
    if(*which == 1) goto S170;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*cdflib_dpmpar(&K1))) goto S160;
    if(!(pq < 0.0e0)) goto S140;
    *bound = 0.0e0;
    goto S150;
S140:
    *bound = 1.0e0;
S150:
    *status = 3;
    return;
S170:
S160:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Computing P and Q
//
        cdflib_cumt(t,df,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Computing T
//     .. Get initial approximation for T
//
        *t = cdflib_dt1(p,q,df);
        T2 = -inf;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        cdflib_dstinv(&T2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        cdflib_dinvr(status,t,&fx,&qleft,&qhi);
S180:
        if(!(*status == 1)) goto S210;
        cdflib_cumt(t,df,&cum,&ccum);
        if(!qporq) goto S190;
        fx = cum-*p;
        goto S200;
S190:
        fx = ccum-*q;
S200:
        cdflib_dinvr(status,t,&fx,&qleft,&qhi);
        goto S180;
S210:
        if(!(*status == -1)) goto S240;
        if(!qleft) goto S220;
        *status = 1;
        *bound = -inf;
        goto S230;
S220:
        *status = 2;
        *bound = inf;
S240:
S230:
        ;
    }
    else if(3 == *which) {
//
//     Computing DF
//
        *df = 5.0e0;
        T8 = zero;
        T9 = maxdf;
        T10 = atol;
        T11 = tol;
        cdflib_dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
        *status = 0;
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
S250:
        if(!(*status == 1)) goto S280;
        cdflib_cumt(t,df,&cum,&ccum);
        if(!qporq) goto S260;
        fx = cum-*p;
        goto S270;
S260:
        fx = ccum-*q;
S270:
        cdflib_dinvr(status,df,&fx,&qleft,&qhi);
        goto S250;
S280:
        if(!(*status == -1)) goto S310;
        if(!qleft) goto S290;
        *status = 1;
        *bound = zero;
        goto S300;
S290:
        *status = 2;
        *bound = maxdf;
S300:
        ;
    }
S310:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
# undef maxdf
}

void
cdflib_chi_noncentral_cdf_values(
    int* n_data,
    double* x,
    double* lambda,
    int* df,
    double* cdf)
{
  enum { N_MAX = 27 };

  static const double cdf_vec[N_MAX] = {
    0.839944E+00, 0.695906E+00, 0.535088E+00,
    0.764784E+00, 0.620644E+00, 0.469167E+00,
    0.307088E+00, 0.220382E+00, 0.150025E+00,
    0.307116E-02, 0.176398E-02, 0.981679E-03,
    0.165175E-01, 0.202342E-03, 0.498448E-06,
    0.151325E-01, 0.209041E-02, 0.246502E-03,
    0.263684E-01, 0.185798E-01, 0.130574E-01,
    0.583804E-01, 0.424978E-01, 0.308214E-01,
    0.105788E+00, 0.794084E-01, 0.593201E-01 };
  static const int df_vec[N_MAX] = {
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
     60,  80, 100,
      1,   2,   3,
     10,  10,  10,
     10,  10,  10,
     10,  10,  10 };
  static const double lambda_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,
     1.0E+00,  1.0E+00,  1.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00,
    30.0E+00, 30.0E+00, 30.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00 };
  static const double x_vec[N_MAX] = {
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
    60.000E+00, 60.000E+00, 60.000E+00,
     0.050E+00,  0.050E+00,  0.050E+00,
     4.000E+00,  4.000E+00,  4.000E+00,
     5.000E+00,  5.000E+00,  5.000E+00,
     6.000E+00,  6.000E+00,  6.000E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *lambda = 0.0E+00;
    *df = 0;
    *cdf = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *df = df_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
}

void
cdflib_chi_square_cdf_values(
    int* n_data,
    int* a,
    double* x,
    double* fx)
{
  enum { N_MAX = 21 };

  static const int a_vec[N_MAX] = {
     1,  2,  1,  2,
     1,  2,  3,  4,
     1,  2,  3,  4,
     5,  3,  3,  3,
     3,  3, 10, 10,
    10 };
  static const double fx_vec[N_MAX] = {
    0.0796557E+00, 0.00498752E+00, 0.112463E+00,    0.00995017E+00,
    0.472911E+00,  0.181269E+00,   0.0597575E+00,   0.0175231E+00,
    0.682689E+00,  0.393469E+00,   0.198748E+00,    0.090204E+00,
    0.0374342E+00, 0.427593E+00,   0.608375E+00,    0.738536E+00,
    0.828203E+00,  0.88839E+00,    0.000172116E+00, 0.00365985E+00,
    0.0185759E+00 };
  static const double x_vec[N_MAX] = {
    0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
    0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
    1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
    1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
    5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
    3.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

void
cdflib_cumbet(
    double* x,
    double* y,
    double* a,
    double* b,
    double* cum,
    double* ccum)
{
  static int ierr;

  if ( *x <= 0.0 )
  {
    *cum = 0.0;
    *ccum = 1.0;
  }
  else if ( *y <= 0.0 )
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  else
  {
    cdflib_beta_inc ( a, b, x, y, cum, ccum, &ierr );
  }
  return;
}

void
cdflib_cumbin(
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    double* cum,
    double* ccum)
{
  static double T1,T2;

  if ( *s < *xn )
  {
    T1 = *s + 1.0;
    T2 = *xn - *s;
    cdflib_cumbet ( pr, ompr, &T1, &T2, ccum, cum );
  }
  else
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  return;
}

void
cdflib_cumchi(
    double* x,
    double* df,
    double* cum,
    double* ccum)
{
  static double a;
  static double xx;

  a = *df * 0.5;
  xx = *x * 0.5;
  cdflib_cumgam ( &xx, &a, cum, ccum );
  return;
}

void
cdflib_cumchn(
    double* x,
    double* df,
    double* pnonc,
    double* cum,
    double* ccum)
{
# define dg(i) (*df+2.0e0*(double)(i))
# define qsmall(xx) (int)(sum < 1.0e-20 || (xx) < eps*sum)
# define qtired(i) (int)((i) > ntired)

  static double eps = 1.0e-5;
  static int ntired = 1000;
  static double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
    sumadj,term,wt,xnonc;
  static int i,icent,iterb,iterf;
  static double T1,T2,T3;

    if(!(*x <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return;
S10:
    if(!(*pnonc <= 1.0e-10)) goto S20;
//
//     When non-centrality parameter is (essentially) zero,
//     use cumulative chi-square distribution
//
    cdflib_cumchi(x,df,cum,ccum);
    return;
S20:
    xnonc = *pnonc/2.0e0;
//
//     The following code calculates the weight, chi-square, and
//     adjustment term for the central term in the infinite series.
//     The central term is the one in which the poisson weight is
//     greatest.  The adjustment term is the amount that must
//     be subtracted from the chi-square to move up two degrees
//     of freedom.
//
    icent = cdflib_fifidint(xnonc);
    if(icent == 0) icent = 1;
    chid2 = *x/2.0e0;
//
//     Calculate central weight term
//
    T1 = (double)(icent+1);
    lfact = cdflib_gamma_log ( &T1 );
    lcntwt = -xnonc+(double)icent*log(xnonc)-lfact;
    centwt = exp(lcntwt);
//
//     Calculate central chi-square
//
    T2 = dg(icent);
    cdflib_cumchi(x,&T2,&pcent,ccum);
//
//     Calculate central adjustment term
//
    dfd2 = dg(icent)/2.0e0;
    T3 = 1.0e0+dfd2;
    lfact = cdflib_gamma_log ( &T3 );
    lcntaj = dfd2*log(chid2)-chid2-lfact;
    centaj = exp(lcntaj);
    sum = centwt*pcent;
//
//     Sum backwards from the central term towards zero.
//     Quit whenever either
//     (1) the zero term is reached, or
//     (2) the term gets small relative to the sum, or
//     (3) More than NTIRED terms are totaled.
//
    iterb = 0;
    sumadj = 0.0e0;
    adj = centaj;
    wt = centwt;
    i = icent;
    goto S40;
S30:
    if(qtired(iterb) || qsmall(term) || i == 0) goto S50;
S40:
    dfd2 = dg(i)/2.0e0;
//
//     Adjust chi-square for two fewer degrees of freedom.
//     The adjusted value ends up in PTERM.
//
    adj = adj*dfd2/chid2;
    sumadj = sumadj + adj;
    pterm = pcent+sumadj;
//
//     Adjust poisson weight for J decreased by one
//
    wt *= ((double)i/xnonc);
    term = wt*pterm;
    sum = sum + term;
    i -= 1;
    iterb = iterb + 1;
    goto S30;
S50:
    iterf = 0;
//
//     Now sum forward from the central term towards infinity.
//     Quit when either
//     (1) the term gets small relative to the sum, or
//     (2) More than NTIRED terms are totaled.
//
    sumadj = adj = centaj;
    wt = centwt;
    i = icent;
    goto S70;
S60:
    if(qtired(iterf) || qsmall(term)) goto S80;
S70:
//
//     Update weights for next higher J
//
    wt *= (xnonc/(double)(i+1));
//
//     Calculate PTERM and add term to sum
//
    pterm = pcent-sumadj;
    term = wt*pterm;
    sum = sum + term;
//
//  Update adjustment term for DF for next iteration
//
    i = i + 1;
    dfd2 = dg(i)/2.0e0;
    adj = adj*chid2/dfd2;
    sumadj = sum + adj;
    iterf = iterf + 1;
    goto S60;
S80:
    *cum = sum;
    *ccum = 0.5e0+(0.5e0-*cum);
    return;
# undef dg
# undef qsmall
# undef qtired
}

void
cdflib_cumf(
    double* f,
    double* dfn,
    double* dfd,
    double* cum,
    double* ccum)
{
# define half 0.5e0
# define done 1.0e0

  static double dsum,prod,xx,yy;
  static int ierr;
  static double T1,T2;

  if(!(*f <= 0.0e0)) goto S10;
  *cum = 0.0e0;
  *ccum = 1.0e0;
  return;
S10:
  prod = *dfn**f;
//
//     XX is such that the incomplete beta with parameters
//     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
//     YY is 1 - XX
//     Calculate the smaller of XX and YY accurately
//
  dsum = *dfd+prod;
  xx = *dfd/dsum;

  if ( xx > half )
  {
    yy = prod/dsum;
    xx = done-yy;
  }
  else
  {
    yy = done-xx;
  }

  T1 = *dfd*half;
  T2 = *dfn*half;
  cdflib_beta_inc ( &T1, &T2, &xx, &yy, ccum, cum, &ierr );
  return;
# undef half
# undef done
}

void
cdflib_cumfnc(
    double* f,
    double* dfn,
    double* dfd,
    double* pnonc,
    double* cum,
    double* ccum)
{
# define qsmall(x) (int)(sum < 1.0e-20 || (x) < eps*sum)
# define half 0.5e0
# define done 1.0e0

  static const double eps = 1.0e-4;
  static double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
    upterm,xmult,xnonc;
  static int i,icent,ierr;
  static double T1,T2,T3,T4,T5,T6;

    if(!(*f <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return;
S10:
    if(!(*pnonc < 1.0e-10)) goto S20;
//
//  Handle case in which the non-centrality parameter is
//  (essentially) zero.
//
    cdflib_cumf(f,dfn,dfd,cum,ccum);
    return;
S20:
    xnonc = *pnonc/2.0e0;
//
//  Calculate the central term of the poisson weighting factor.
//
    icent = ( int ) xnonc;
    if(icent == 0) icent = 1;
//
//  Compute central weight term
//
    T1 = (double)(icent+1);
    centwt = exp(-xnonc+(double)icent*log(xnonc)- cdflib_gamma_log ( &T1 ) );
//
//  Compute central incomplete beta term
//  Assure that minimum of arg to beta and 1 - arg is computed
//  accurately.
//
    prod = *dfn**f;
    dsum = *dfd+prod;
    yy = *dfd/dsum;
    if(yy > half) {
        xx = prod/dsum;
        yy = done-xx;
    }
    else  xx = done-yy;
    T2 = *dfn*half+(double)icent;
    T3 = *dfd*half;
    cdflib_beta_inc ( &T2, &T3, &xx, &yy, &betdn, &dummy, &ierr );
    adn = *dfn/2.0e0+(double)icent;
    aup = adn;
    b = *dfd/2.0e0;
    betup = betdn;
    sum = centwt*betdn;
//
//  Now sum terms backward from icent until convergence or all done
//
    xmult = centwt;
    i = icent;
    T4 = adn+b;
    T5 = adn+1.0e0;
    dnterm = exp( cdflib_gamma_log ( &T4 ) - cdflib_gamma_log ( &T5 )
      - cdflib_gamma_log ( &b ) + adn * log ( xx ) + b * log(yy));
S30:
    if(qsmall(xmult*betdn) || i <= 0) goto S40;
    xmult *= ((double)i/xnonc);
    i -= 1;
    adn -= 1.0;
    dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
    betdn += dnterm;
    sum += (xmult*betdn);
    goto S30;
S40:
    i = icent+1;
//
//  Now sum forwards until convergence
//
    xmult = centwt;
    if(aup-1.0+b == 0) upterm = exp(-cdflib_gamma_log ( &aup )
      - cdflib_gamma_log ( &b ) + (aup-1.0)*log(xx)+
      b*log(yy));
    else  {
        T6 = aup-1.0+b;
        upterm = exp( cdflib_gamma_log ( &T6 ) - cdflib_gamma_log ( &aup )
          - cdflib_gamma_log ( &b ) + (aup-1.0)*log(xx)+b*
          log(yy));
    }
    goto S60;
S50:
    if(qsmall(xmult*betup)) goto S70;
S60:
    xmult *= (xnonc/(double)i);
    i += 1;
    aup += 1.0;
    upterm = (aup+b-2.0e0)*xx/(aup-1.0)*upterm;
    betup -= upterm;
    sum += (xmult*betup);
    goto S50;
S70:
    *cum = sum;
    *ccum = 0.5e0+(0.5e0-*cum);
    return;
# undef qsmall
# undef half
# undef done
}

void
cdflib_cumgam(
    double* x,
    double* a,
    double* cum,
    double* ccum)
{
  static int K1 = 0;

  if(!(*x <= 0.0e0)) goto S10;
  *cum = 0.0e0;
  *ccum = 1.0e0;
  return;
S10:
  cdflib_gamma_inc ( a, x, cum, ccum, &K1 );
//
//     Call gratio routine
//
    return;
}

void
cdflib_cumnbn(
    double* s,
    double* xn,
    double* pr,
    double* ompr,
    double* cum,
    double* ccum)
{
  static double T1;

  T1 = *s+1.e0;
  cdflib_cumbet(pr,ompr,xn,&T1,cum,ccum);
  return;
}

void
cdflib_cumnor(
    double* arg,
    double* result,
    double* ccum)
{
  static const double a[5] = {
    2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
  };
  static const double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
  };
  static const double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
  };
  static const double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
  };
  static const double half = 0.5e0;
  static const double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
  };
  static const double one = 1.0e0;
  static const double q[5] = {
    1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
  };
  static const double sixten = 1.60e0;
  static const double sqrpi = 3.9894228040143267794e-1;
  static const double thrsh = 0.66291e0;
  static const double root32 = 5.656854248e0;
  static const double zero = 0.0e0;
  static int K1 = 1;
  static int K2 = 2;
  static int i;
  static double del,eps,temp,x,xden,xnum,y,xsq,min;
//
//  Machine dependent constants
//
    eps = cdflib_dpmpar(&K1)*0.5e0;
    min = cdflib_dpmpar(&K2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh) {
//
//  Evaluate  anorm  for  |X| <= 0.66291
//
        xsq = zero;
        if(y > eps) xsq = x*x;
        xnum = a[4]*xsq;
        xden = xsq;
        for ( i = 0; i < 3; i++ )
        {
            xnum = (xnum+a[i])*xsq;
            xden = (xden+b[i])*xsq;
        }
        *result = x*(xnum+a[3])/(xden+b[3]);
        temp = *result;
        *result = half+temp;
        *ccum = half-temp;
    }
//
//  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
//
    else if(y <= root32) {
        xnum = c[8]*y;
        xden = y;
        for ( i = 0; i < 7; i++ )
        {
            xnum = (xnum+c[i])*y;
            xden = (xden+d[i])*y;
        }
        *result = (xnum+c[7])/(xden+d[7]);
        xsq = cdflib_fifdint(y*sixten)/sixten;
        del = (y-xsq)*(y+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
//
//  Evaluate  anorm  for |X| > sqrt(32)
//
    else  {
        *result = zero;
        xsq = one/(x*x);
        xnum = p[5]*xsq;
        xden = xsq;
        for ( i = 0; i < 4; i++ )
        {
            xnum = (xnum+p[i])*xsq;
            xden = (xden+q[i])*xsq;
        }
        *result = xsq*(xnum+p[4])/(xden+q[4]);
        *result = (sqrpi-*result)/y;
        xsq = cdflib_fifdint(x*sixten)/sixten;
        del = (x-xsq)*(x+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    if(*result < min) *result = 0.0e0;
//
//  Fix up for negative argument, erf, etc.
//
    if(*ccum < min) *ccum = 0.0e0;
}

void
cdflib_cumpoi(
    double* s,
    double* xlam,
    double* cum,
    double* ccum)
{
  static double chi,df;

  df = 2.0e0*(*s+1.0e0);
  chi = 2.0e0**xlam;
  cdflib_cumchi(&chi,&df,ccum,cum);
  return;
}

void
cdflib_cumt(
    double* t,
    double* df,
    double* cum,
    double* ccum)
{
  static double a;
  static double dfptt;
  static double K2 = 0.5e0;
  static double oma;
  static double T1;
  static double tt;
  static double xx;
  static double yy;

  tt = (*t) * (*t);
  dfptt = ( *df ) + tt;
  xx = *df / dfptt;
  yy = tt / dfptt;
  T1 = 0.5e0 * ( *df );
  cdflib_cumbet ( &xx, &yy, &T1, &K2, &a, &oma );

  if ( *t <= 0.0e0 )
  {
    *cum = 0.5e0 * a;
    *ccum = oma + ( *cum );
  }
  else
  {
    *ccum = 0.5e0 * a;
    *cum = oma + ( *ccum );
  }
  return;
}

double
cdflib_dbetrm(
    double* a,
    double* b)
{
  static double dbetrm,T1,T2,T3;
//
//     Try to sum from smallest to largest
//
    T1 = *a+*b;
    dbetrm = -cdflib_dstrem(&T1);
    T2 = cdflib_fifdmax1(*a,*b);
    dbetrm += cdflib_dstrem(&T2);
    T3 = cdflib_fifdmin1(*a,*b);
    dbetrm += cdflib_dstrem(&T3);
    return dbetrm;
}

double
cdflib_dexpm1(
    double* x)
{
  static const double p1 = .914041914819518e-09;
  static const double p2 = .238082361044469e-01;
  static const double q1 = -.499999999085958e+00;
  static const double q2 = .107141568980644e+00;
  static const double q3 = -.119041179760821e-01;
  static const double q4 = .595130811860248e-03;
  static double dexpm1;
  double w;

  if ( fabs(*x) <= 0.15e0 )
  {
    dexpm1 =   *x * ( ( (
        p2   * *x
      + p1 ) * *x
      + 1.0e0 )
      /((((
        q4   * *x
      + q3 ) * *x
      + q2 ) * *x
      + q1 ) * *x
      + 1.0e0 ) );
  }
  else if ( *x <= 0.0e0 )
  {
    w = exp(*x);
    dexpm1 = w-0.5e0-0.5e0;
  }
  else
  {
    w = exp(*x);
    dexpm1 = w*(0.5e0+(0.5e0-1.0e0/w));
  }

  return dexpm1;
}

double
cdflib_dinvnr(
    double* p,
    double* q)
{
# define maxit 100
# define eps (1.0e-13)
# define r2pi 0.3989422804014326e0
# define nhalf (-0.5e0)
# define dennor(x) (r2pi*exp(nhalf*(x)*(x)))

  static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
  static int i;
  static unsigned long qporq;

//
//     FIND MINIMUM OF P AND Q
//
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
S10:
    pp = *q;
S20:
//
//     INITIALIZATION STEP
//
    strtx = cdflib_stvaln(&pp);
    xcur = strtx;
//
//     NEWTON INTERATIONS
//
    for ( i = 1; i <= maxit; i++ )
    {
        cdflib_cumnor(&xcur,&cum,&ccum);
        dx = (cum-pp)/dennor(xcur);
        xcur -= dx;
        if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
//
//     IF WE GET HERE, NEWTON HAS FAILED
//
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
S40:
//
//     IF WE GET HERE, NEWTON HAS SUCCEDED
//
    dinvnr = xcur;
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
# undef maxit
# undef eps
# undef r2pi
# undef nhalf
# undef dennor
}

void
cdflib_dinvr(
    int* status,
    double* x,
    double* fx,
    unsigned long* qleft,
    unsigned long* qhi)
{
  cdflib_E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}

double
cdflib_dlanor(
    double* x)
{
# define dlsqpi 0.91893853320467274177e0

  static const double coef[12] = {
    -1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
    -34459425.0e0,654729075.0e0,-13749310575.e0,316234143225.0e0
  };
  static int K1 = 12;
  static double dlanor,approx,correc,xx,xx2,T2;

  xx = fabs(*x);
  if ( xx < 5.0e0 )
  {
    cdflib_ftnstop(" Argument too small in DLANOR");
  }
  approx = -dlsqpi-0.5e0*xx*xx-log(xx);
  xx2 = xx*xx;
  T2 = 1.0e0/xx2;
  correc = cdflib_eval_pol ( coef, &K1, &T2 ) / xx2;
  correc = cdflib_alnrel ( &correc );
  dlanor = approx+correc;
  return dlanor;
# undef dlsqpi
}

double
cdflib_dpmpar(
    int* i)
{
  static int K1 = 4;
  static int K2 = 8;
  static int K3 = 9;
  static int K4 = 10;
  static double value,b,binv,bm1,one,w,z;
  static int emax,emin,ibeta,m;

    if(*i > 1) goto S10;
    b = cdflib_ipmpar(&K1);
    m = cdflib_ipmpar(&K2);
    value = pow(b,(double)(1-m));
    return value;
S10:
    if(*i > 2) goto S20;
    b = cdflib_ipmpar(&K1);
    emin = cdflib_ipmpar(&K3);
    one = 1.0;
    binv = one/b;
    w = pow(b,(double)(emin+2));
    value = w*binv*binv*binv;
    return value;
S20:
    ibeta = cdflib_ipmpar(&K1);
    m = cdflib_ipmpar(&K2);
    emax = cdflib_ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1.0;
    z = pow(b,(double)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(double)(emax-2));
    value = w*z*b*b;
    return value;
}

void
cdflib_dstinv(
    double* zsmall,
    double* zbig,
    double* zabsst,
    double* zrelst,
    double* zstpmu,
    double* zabsto,
    double* zrelto)
{
  cdflib_E0000(1,NULL,NULL,NULL,NULL,NULL,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
    zstpmu);
}

double
cdflib_dstrem(
    double* z)
{
# define hln2pi 0.91893853320467274178e0
# define ncoef 10

  static const double coef[ncoef] = {
    0.0e0,0.0833333333333333333333333333333e0,
    -0.00277777777777777777777777777778e0,0.000793650793650793650793650793651e0,
    -0.000595238095238095238095238095238e0,
    0.000841750841750841750841750841751e0,-0.00191752691752691752691752691753e0,
    0.00641025641025641025641025641026e0,-0.0295506535947712418300653594771e0,
    0.179644372368830573164938490016e0
  };
  static int K1 = 10;
  static double dstrem,sterl,T2;
//
//    For information, here are the next 11 coefficients of the
//    remainder term in Sterling's formula
//            -1.39243221690590111642743221691
//            13.4028640441683919944789510007
//            -156.848284626002017306365132452
//            2193.10333333333333333333333333
//            -36108.7712537249893571732652192
//            691472.268851313067108395250776
//            -0.152382215394074161922833649589D8
//            0.382900751391414141414141414141D9
//            -0.108822660357843910890151491655D11
//            0.347320283765002252252252252252D12
//            -0.123696021422692744542517103493D14
//
    if(*z <= 0.0e0)
    {
      cdflib_ftnstop ( "Zero or negative argument in DSTREM" );
    }
    if(!(*z > 6.0e0)) goto S10;
    T2 = 1.0e0/pow(*z,2.0);
    dstrem = cdflib_eval_pol ( coef, &K1, &T2 )**z;
    goto S20;
S10:
    sterl = hln2pi+(*z-0.5e0)*log(*z)-*z;
    dstrem = cdflib_gamma_log ( z ) - sterl;
S20:
    return dstrem;
# undef hln2pi
# undef ncoef
}

void
cdflib_dstzr(
    double* zxlo,
    double* zxhi,
    double* zabstl,
    double* zreltl)
{
  cdflib_E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
}

double
cdflib_dt1(
    double* p,
    double* q,
    double* df)
{
  static const double coef[4][5] = {
    {1.0e0,1.0e0,0.0e0,0.0e0,0.0e0},{3.0e0,16.0e0,5.0e0,0.0e0,0.0e0},{-15.0e0,17.0e0,
    19.0e0,3.0e0,0.0e0},{-945.0e0,-1920.0e0,1482.0e0,776.0e0,79.0e0}
  };
  static const double denom[4] = {
    4.0e0,96.0e0,384.0e0,92160.0e0
  };
  static const int ideg[4] = {
    2,3,4,5
  };
  static double dt1,denpow,sum,term,x,xp,xx;
  static int i;

    x = fabs(cdflib_dinvnr(p,q));
    xx = x*x;
    sum = x;
    denpow = 1.0e0;
    for ( i = 0; i < 4; i++ )
    {
        term = cdflib_eval_pol ( &coef[i][0], &ideg[i], &xx ) * x;
        denpow *= *df;
        sum += (term/(denpow*denom[i]));
    }
    if(!(*p >= 0.5e0)) goto S20;
    xp = sum;
    goto S30;
S20:
    xp = -sum;
S30:
    dt1 = xp;
    return dt1;
}

void
cdflib_dzror(
    int* status,
    double* x,
    double* fx,
    double* xlo,
    double* xhi,
    unsigned long* qleft,
    unsigned long* qhi)
{
  cdflib_E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
}

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
    double* zstpmu)
{
# define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))

  static double absstp;
  static double abstol;
  static double big,fbig,fsmall,relstp,reltol,small,step,stpmul,xhi,
    xlb,xlo,xsave,xub,yy;
  static int i99999;
  static unsigned long qbdd,qcond,qdum1,qdum2,qincr,qlim,qup;
    switch(IENTRY){case 0: goto DINVR; case 1: goto DSTINV;}
DINVR:
    if(*status > 0) goto S310;
    qcond = !qxmon(small,*x,big);
    if(qcond)
    {
      cdflib_ftnstop(" SMALL, X, BIG not monotone in INVR");
    }
    xsave = *x;
//
//     See that SMALL and BIG bound the zero and set QINCR
//
    *x = small;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 1;
    goto S300;
S10:
    fsmall = *fx;
    *x = big;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 2;
    goto S300;
S20:
    fbig = *fx;
    qincr = fbig > fsmall;
    if(!qincr) goto S50;
    if(fsmall <= 0.0e0) goto S30;
    *status = -1;
    *qleft = *qhi = 1;
    return;
S30:
    if(fbig >= 0.0e0) goto S40;
    *status = -1;
    *qleft = *qhi = 0;
    return;
S40:
    goto S80;
S50:
    if(fsmall >= 0.0e0) goto S60;
    *status = -1;
    *qleft = 1;
    *qhi = 0;
    return;
S60:
    if(fbig <= 0.0e0) goto S70;
    *status = -1;
    *qleft = 0;
    *qhi = 1;
    return;
S80:
S70:
    *x = xsave;
    step = cdflib_fifdmax1(absstp,relstp*fabs(*x));
//
//      YY = F(X) - Y
//     GET-FUNCTION-VALUE
//
    i99999 = 3;
    goto S300;
S90:
    yy = *fx;
    if(!(yy == 0.0e0)) goto S100;
    *status = 0;
    return;
S100:
    qup = (qincr && yy < 0.0e0) || (!qincr && yy > 0.0e0);
//
//     HANDLE CASE IN WHICH WE MUST STEP HIGHER
//
    if(!qup) goto S170;
    xlb = xsave;
    xub = cdflib_fifdmin1(xlb+step,big);
    goto S120;
S110:
    if(qcond) goto S150;
S120:
//
//      YY = F(XUB) - Y
//
    *x = xub;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 4;
    goto S300;
S130:
    yy = *fx;
    qbdd = (qincr && yy >= 0.0e0) || (!qincr && yy <= 0.0e0);
    qlim = xub >= big;
    qcond = qbdd || qlim;
    if(qcond) goto S140;
    step = stpmul*step;
    xlb = xub;
    xub = cdflib_fifdmin1(xlb+step,big);
S140:
    goto S110;
S150:
    if(!(qlim && !qbdd)) goto S160;
    *status = -1;
    *qleft = 0;
    *qhi = !qincr;
    *x = big;
    return;
S160:
    goto S240;
S170:
//
//     HANDLE CASE IN WHICH WE MUST STEP LOWER
//
    xub = xsave;
    xlb = cdflib_fifdmax1(xub-step,small);
    goto S190;
S180:
    if(qcond) goto S220;
S190:
//
//      YY = F(XLB) - Y
//
    *x = xlb;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 5;
    goto S300;
S200:
    yy = *fx;
    qbdd = (qincr && yy <= 0.0e0) || (!qincr && yy >= 0.0e0);
    qlim = xlb <= small;
    qcond = qbdd || qlim;
    if(qcond) goto S210;
    step = stpmul*step;
    xub = xlb;
    xlb = cdflib_fifdmax1(xub-step,small);
S210:
    goto S180;
S220:
    if(!(qlim && !qbdd)) goto S230;
    *status = -1;
    *qleft = 1;
    *qhi = qincr;
    *x = small;
    return;
S240:
S230:
    cdflib_dstzr(&xlb,&xub,&abstol,&reltol);
//
//  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
//
    *status = 0;
    goto S260;
S250:
    if(!(*status == 1)) goto S290;
S260:
    cdflib_dzror ( status, x, fx, &xlo, &xhi, &qdum1, &qdum2 );
    if(!(*status == 1)) goto S280;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 6;
    goto S300;
S280:
S270:
    goto S250;
S290:
    *x = xlo;
    *status = 0;
    return;
DSTINV:
    small = *zsmall;
    big = *zbig;
    absstp = *zabsst;
    relstp = *zrelst;
    stpmul = *zstpmu;
    abstol = *zabsto;
    reltol = *zrelto;
    return;
S300:
//
//     TO GET-FUNCTION-VALUE
//
    *status = 1;
    return;
S310:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S90;case
      4: goto S130;case 5: goto S200;case 6: goto S270;default: break;}
# undef qxmon
}

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
    double* zxlo)
{
# define ftol(zx) (0.5e0*cdflib_fifdmax1(abstol,reltol*fabs((zx))))

  static double a,abstol,b,c,d,fa,fb,fc,fd,fda;
  static double fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
  static int ext,i99999;
  static unsigned long first,qrzero;
    switch(IENTRY){case 0: goto DZROR; case 1: goto DSTZR;}
DZROR:
    if(*status > 0) goto S280;
    *xlo = xxlo;
    *xhi = xxhi;
    b = *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 1;
    goto S270;
S10:
    fb = *fx;
    *xlo = *xhi;
    a = *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 2;
    goto S270;
S20:
//
//     Check that F(ZXLO) < 0 < F(ZXHI)  or
//                F(ZXLO) > 0 > F(ZXHI)
//
    if(!(fb < 0.0e0)) goto S40;
    if(!(*fx < 0.0e0)) goto S30;
    *status = -1;
    *qleft = *fx < fb;
    *qhi = 0;
    return;
S40:
S30:
    if(!(fb > 0.0e0)) goto S60;
    if(!(*fx > 0.0e0)) goto S50;
    *status = -1;
    *qleft = *fx > fb;
    *qhi = 1;
    return;
S60:
S50:
    fa = *fx;
    first = 1;
S70:
    c = a;
    fc = fa;
    ext = 0;
S80:
    if(!(fabs(fc) < fabs(fb))) goto S100;
    if(!(c != a)) goto S90;
    d = a;
    fd = fa;
S90:
    a = b;
    fa = fb;
    *xlo = c;
    b = *xlo;
    fb = fc;
    c = a;
    fc = fa;
S100:
    tol = ftol(*xlo);
    m = (c+b)*.5e0;
    mb = m-b;
    if(!(fabs(mb) > tol)) goto S240;
    if(!(ext > 3)) goto S110;
    w = mb;
    goto S190;
S110:
    tol = cdflib_fifdsign(tol,mb);
    p = (b-a)*fb;
    if(!first) goto S120;
    q = fa-fb;
    first = 0;
    goto S130;
S120:
    fdb = (fd-fb)/(d-b);
    fda = (fd-fa)/(d-a);
    p = fda*p;
    q = fdb*fa-fda*fb;
S130:
    if(!(p < 0.0e0)) goto S140;
    p = -p;
    q = -q;
S140:
    if(ext == 3) p *= 2.0e0;
    if(!(p*1.0e0 == 0.0e0 || p <= q*tol)) goto S150;
    w = tol;
    goto S180;
S150:
    if(!(p < mb*q)) goto S160;
    w = p/q;
    goto S170;
S160:
    w = mb;
S190:
S180:
S170:
    d = a;
    fd = fa;
    a = b;
    fa = fb;
    b += w;
    *xlo = b;
    *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 3;
    goto S270;
S200:
    fb = *fx;
    if(!(fc*fb >= 0.0e0)) goto S210;
    goto S70;
S210:
    if(!(w == mb)) goto S220;
    ext = 0;
    goto S230;
S220:
    ext += 1;
S230:
    goto S80;
S240:
    *xhi = c;
    qrzero = (fc >= 0.0e0 && fb <= 0.0e0) || (fc < 0.0e0 && fb >= 0.0e0);
    if(!qrzero) goto S250;
    *status = 0;
    goto S260;
S250:
    *status = -1;
S260:
    return;
DSTZR:
    xxlo = *zxlo;
    xxhi = *zxhi;
    abstol = *zabstl;
    reltol = *zreltl;
    return;
S270:
//
//     TO GET-FUNCTION-VALUE
//
    *status = 1;
    return;
S280:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
      default: break;}
# undef ftol
}

void
cdflib_erf_values(
    int* n_data,
    double* x,
    double* fx)
{
  enum { N_MAX = 21 };

  static const double fx_vec[N_MAX] = {
    0.0000000000E+00, 0.1124629160E+00, 0.2227025892E+00, 0.3286267595E+00,
    0.4283923550E+00, 0.5204998778E+00, 0.6038560908E+00, 0.6778011938E+00,
    0.7421009647E+00, 0.7969082124E+00, 0.8427007929E+00, 0.8802050696E+00,
    0.9103139782E+00, 0.9340079449E+00, 0.9522851198E+00, 0.9661051465E+00,
    0.9763483833E+00, 0.9837904586E+00, 0.9890905016E+00, 0.9927904292E+00,
    0.9953222650E+00 };
  static const double x_vec[N_MAX] = {
    0.0E+00, 0.1E+00, 0.2E+00, 0.3E+00,
    0.4E+00, 0.5E+00, 0.6E+00, 0.7E+00,
    0.8E+00, 0.9E+00, 1.0E+00, 1.1E+00,
    1.2E+00, 1.3E+00, 1.4E+00, 1.5E+00,
    1.6E+00, 1.7E+00, 1.8E+00, 1.9E+00,
    2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_error_f(
    double* x)
{
  static const double c = .564189583547756e0;
  static const double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  static const double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  static const double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  static const double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  static const double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  static const double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  static double erf1,ax,bot,t,top,x2;

    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erf1 = *x*(top/bot);
    return erf1;
S10:
    if(ax > 4.0e0) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
      7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
      7];
    erf1 = 0.5e0+(0.5e0-exp(-(*x**x))*top/bot);
    if(*x < 0.0e0) erf1 = -erf1;
    return erf1;
S20:
    if(ax >= 5.8e0) goto S30;
    x2 = *x**x;
    t = 1.0e0/x2;
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
    erf1 = (c-top/(x2*bot))/ax;
    erf1 = 0.5e0+(0.5e0-exp(-x2)*erf1);
    if(*x < 0.0e0) erf1 = -erf1;
    return erf1;
S30:
    erf1 = cdflib_fifdsign(1.0e0,*x);
    return erf1;
}

double
cdflib_error_fc(
    int* ind,
    double* x)
{
  static const double c = .564189583547756e0;
  static const double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  static const double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  static const double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  static const double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  static const double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  static const double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  static int K1 = 1;
  static double erfc1,ax,bot,e,t,top,w;

//
//                     ABS(X) .LE. 0.5
//
    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erfc1 = 0.5e0+(0.5e0-*x*(top/bot));
    if(*ind != 0) erfc1 = exp(t)*erfc1;
    return erfc1;
S10:
//
//                  0.5 .LT. ABS(X) .LE. 4
//
    if(ax > 4.0e0) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
      7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
      7];
    erfc1 = top/bot;
    goto S40;
S20:
//
//                      ABS(X) .GT. 4
//
    if(*x <= -5.6e0) goto S60;
    if(*ind != 0) goto S30;
    if(*x > 100.0e0) goto S70;
    if(*x**x > -cdflib_exparg(&K1)) goto S70;
S30:
    t = pow(1.0e0/ *x,2.0);
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
    erfc1 = (c-t*top/bot)/ax;
S40:
//
//                      FINAL ASSEMBLY
//
    if(*ind == 0) goto S50;
    if(*x < 0.0e0) erfc1 = 2.0e0*exp(*x**x)-erfc1;
    return erfc1;
S50:
    w = *x**x;
    t = w;
    e = w-t;
    erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
    if(*x < 0.0e0) erfc1 = 2.0e0-erfc1;
    return erfc1;
S60:
//
//             LIMIT VALUE FOR LARGE NEGATIVE X
//
    erfc1 = 2.0e0;
    if(*ind != 0) erfc1 = 2.0e0*exp(*x**x);
    return erfc1;
S70:
//
//             LIMIT VALUE FOR LARGE POSITIVE X
//                       WHEN IND = 0
//
    erfc1 = 0.0e0;
    return erfc1;
}

double
cdflib_esum(
    int* mu,
    double* x)
{
  static double esum,w;

    if(*x > 0.0e0) goto S10;
    if(*mu < 0) goto S20;
    w = (double)*mu+*x;
    if(w > 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S10:
    if(*mu > 0) goto S20;
    w = (double)*mu+*x;
    if(w < 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S20:
    w = *mu;
    esum = exp(w)*exp(*x);
    return esum;
}

double
cdflib_eval_pol(
    const double a[],
    const int* n,
    double* x)
{
  static double devlpl,term;
  static int i;

  term = a[*n-1];
  for ( i = *n-1-1; i >= 0; i-- )
  {
    term = a[i]+term**x;
  }

  devlpl = term;
  return devlpl;
}

double
cdflib_exparg(
    int* l)
{
  static int K1 = 4;
  static int K2 = 9;
  static int K3 = 10;
  static double exparg,lnb;
  static int b,m;

    b = cdflib_ipmpar(&K1);
    if(b != 2) goto S10;
    lnb = .69314718055995e0;
    goto S40;
S10:
    if(b != 8) goto S20;
    lnb = 2.0794415416798e0;
    goto S40;
S20:
    if(b != 16) goto S30;
    lnb = 2.7725887222398e0;
    goto S40;
S30:
    lnb = log((double)b);
S40:
    if(*l == 0) goto S50;
    m = cdflib_ipmpar(&K2)-1;
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
S50:
    m = cdflib_ipmpar(&K3);
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
}

void
cdflib_f_cdf_values(
    int* n_data,
    int* a,
    int* b,
    double* x,
    double* fx)
{
  enum { N_MAX = 20 };

  static const int a_vec[N_MAX] = {
    1, 1, 5, 1,
    2, 4, 1, 6,
    8, 1, 3, 6,
    1, 1, 1, 1,
    2, 3, 4, 5 };
  static const int b_vec[N_MAX] = {
     1,  5,  1,  5,
    10, 20,  5,  6,
    16,  5, 10, 12,
     5,  5,  5,  5,
     5,  5,  5,  5 };
  static const double fx_vec[N_MAX] = {
    0.500000E+00, 0.499971E+00, 0.499603E+00, 0.749699E+00,
    0.750466E+00, 0.751416E+00, 0.899987E+00, 0.899713E+00,
    0.900285E+00, 0.950025E+00, 0.950057E+00, 0.950193E+00,
    0.975013E+00, 0.990002E+00, 0.994998E+00, 0.999000E+00,
    0.568799E+00, 0.535145E+00, 0.514343E+00, 0.500000E+00 };
  static const double x_vec[N_MAX] = {
    1.00E+00,  0.528E+00, 1.89E+00,  1.69E+00,
    1.60E+00,  1.47E+00,  4.06E+00,  3.05E+00,
    2.09E+00,  6.61E+00,  3.71E+00,  3.00E+00,
   10.01E+00, 16.26E+00, 22.78E+00, 47.18E+00,
    1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

void
cdflib_f_noncentral_cdf_values(
    int* n_data,
    int* a,
    int* b,
    double* lambda,
    double* x,
    double* fx)
{
  enum { N_MAX = 22 };

  static const int a_vec[N_MAX] = {
     1,  1,  1,  1,
     1,  1,  1,  1,
     1,  1,  2,  2,
     3,  3,  4,  4,
     5,  5,  6,  6,
     8, 16 };
  static const int b_vec[N_MAX] = {
     1,  5,  5,  5,
     5,  5,  5,  5,
     5,  5,  5, 10,
     5,  5,  5,  5,
     1,  5,  6, 12,
    16,  8 };
  static const double fx_vec[N_MAX] = {
    0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
    0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
    0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
    0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
    0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
    0.844589E+00, 0.816368E+00 };
  static const double lambda_vec[N_MAX] = {
    0.00E+00,  0.000E+00, 0.25E+00,  1.00E+00,
    1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    0.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00 };
  static const double x_vec[N_MAX] = {
    1.00E+00,  1.00E+00, 1.00E+00,  0.50E+00,
    1.00E+00,  2.00E+00, 3.00E+00,  4.00E+00,
    5.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  2.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    2.00E+00,  2.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *lambda = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
}

double
cdflib_fifdint(
    double a)
{
  return (double) ((int) a);
}

double
cdflib_fifdmax1(
    double a,
    double b)
{
  if ( a < b )
  {
    return b;
  }
  else
  {
    return a;
  }
}

double
cdflib_fifdmin1(
    double a,
    double b)
{
  if (a < b) return a;
  else return b;
}

double
cdflib_fifdsign(
    double mag,
    double sign)
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;

}

long
cdflib_fifidint(
    double a)
{
  if ( a < 1.0 )
  {
    return (long) 0;
  }
  else
  {
    return ( long ) a;
  }
}

long
cdflib_fifmod(
    long a,
    long b)
{
  return ( a % b );
}

double
cdflib_fpser(
    double* a,
    double* b,
    double* x,
    double* eps)
{
  static int K1 = 1;
  static double fpser,an,c,s,t,tol;

    fpser = 1.0e0;
    if(*a <= 1.e-3**eps) goto S10;
    fpser = 0.0e0;
    t = *a*log(*x);
    if(t < cdflib_exparg(&K1)) return fpser;
    fpser = exp(t);
S10:
//
//                NOTE THAT 1/B(A,B) = B
//
    fpser = *b/ *a*fpser;
    tol = *eps/ *a;
    an = *a+1.0e0;
    t = *x;
    s = t/an;
S20:
    an += 1.0e0;
    t = *x*t;
    c = t/an;
    s += c;
    if(fabs(c) > tol) goto S20;
    fpser *= (1.0e0+*a*s);
    return fpser;
}

void
cdflib_ftnstop(
    char* msg)
{
  printf ( "%s\n", msg );

  exit ( 0 );
}

double
cdflib_gam1(
    double* a)
{
  static const double s1 = .273076135303957e+00;
  static const double s2 = .559398236957378e-01;
  static const double p[7] = {
    .577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
    .597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
    .589597428611429e-03
  };
  static const double q[5] = {
    .100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
    .261132021441447e-01,.423244297896961e-02
  };
  static const double r[9] = {
    -.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
    .118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
    .223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
  };
  static double gam1,bot,d,t,top,w,T1;

    t = *a;
    d = *a-0.5e0;
    if(d > 0.0e0) t = d-0.5e0;
    T1 = t;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S10;
    else  goto S20;
S10:
    gam1 = 0.0e0;
    return gam1;
S20:
    top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
    bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0e0;
    w = top/bot;
    if(d > 0.0e0) goto S30;
    gam1 = *a*w;
    return gam1;
S30:
    gam1 = t/ *a*(w-0.5e0-0.5e0);
    return gam1;
S40:
    top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
      r[0];
    bot = (s2*t+s1)*t+1.0e0;
    w = top/bot;
    if(d > 0.0e0) goto S50;
    gam1 = *a*(w+0.5e0+0.5e0);
    return gam1;
S50:
    gam1 = t*w/ *a;
    return gam1;
}

void
cdflib_gamma_inc(
    double* a,
    double* x,
    double* ans,
    double* qans,
    int* ind)
{
  static const double alog10 = 2.30258509299405e0;
  static const double d10 = -.185185185185185e-02;
  static const double d20 = .413359788359788e-02;
  static const double d30 = .649434156378601e-03;
  static const double d40 = -.861888290916712e-03;
  static const double d50 = -.336798553366358e-03;
  static const double d60 = .531307936463992e-03;
  static const double d70 = .344367606892378e-03;
  static const double rt2pin = .398942280401433e0;
  static const double rtpi = 1.77245385090552e0;
  static const double third = .333333333333333e0;
  static const double acc0[3] = {
    5.e-15,5.e-7,5.e-4
  };
  static const double big[3] = {
    20.0e0,14.0e0,10.0e0
  };
  static const double d0[13] = {
    .833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
    .352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
    -.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
    -.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
    -.438203601845335e-08
  };
  static const double d1[12] = {
    -.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
    .205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
    .764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
    .137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
  };
  static const double d2[10] = {
    -.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
    -.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
    .342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
    .142806142060642e-06
  };
  static const double d3[8] = {
    .229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
    -.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
    -.567495282699160e-05,.142309007324359e-05
  };
  static const double d4[6] = {
    .784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
    .664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
  };
  static const double d5[4] = {
    -.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
    .679778047793721e-04
  };
  static const double d6[2] = {
    -.592166437353694e-03,.270878209671804e-03
  };
  static const double e00[3] = {
    .25e-3,.25e-1,.14e0
  };
  static const double x00[3] = {
    31.0e0,17.0e0,9.7e0
  };
  static int K1 = 1;
  static int K2 = 0;
  static double a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,
    cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
  static int i,iop,m,max,n;
  static double wk[20],T3;
  static int T4,T5;
  static double T6,T7;

//
//  E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
//  NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
//
    e = cdflib_dpmpar(&K1);
    if(*a < 0.0e0 || *x < 0.0e0) goto S430;
    if(*a == 0.0e0 && *x == 0.0e0) goto S430;
    if(*a**x == 0.0e0) goto S420;
    iop = *ind+1;
    if(iop != 1 && iop != 2) iop = 3;
    acc = cdflib_fifdmax1(acc0[iop-1],e);
    e0 = e00[iop-1];
    x0 = x00[iop-1];
//
//  SELECT THE APPROPRIATE ALGORITHM
//
    if(*a >= 1.0e0) goto S10;
    if(*a == 0.5e0) goto S390;
    if(*x < 1.1e0) goto S160;
    t1 = *a*log(*x)-*x;
    u = *a*exp(t1);
    if(u == 0.0e0) goto S380;
    r = u*(1.0e0+cdflib_gam1(a));
    goto S250;
S10:
    if(*a >= big[iop-1]) goto S30;
    if(*a > *x || *x >= x0) goto S20;
    twoa = *a+*a;
    m = cdflib_fifidint(twoa);
    if(twoa != (double)m) goto S20;
    i = m/2;
    if(*a == (double)i) goto S210;
    goto S220;
S20:
    t1 = *a*log(*x)-*x;
    r = exp(t1)/ cdflib_gamma_x(a);
    goto S40;
S30:
    l = *x/ *a;
    if(l == 0.0e0) goto S370;
    s = 0.5e0+(0.5e0-l);
    z = cdflib_rlog(&l);
    if(z >= 700.0e0/ *a) goto S410;
    y = *a*z;
    rta = sqrt(*a);
    if(fabs(s) <= e0/rta) goto S330;
    if(fabs(s) <= 0.4e0) goto S270;
    t = pow(1.0e0/ *a,2.0);
    t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
    t1 -= y;
    r = rt2pin*rta*exp(t1);
S40:
    if(r == 0.0e0) goto S420;
    if(*x <= cdflib_fifdmax1(*a,alog10)) goto S50;
    if(*x < x0) goto S250;
    goto S100;
S50:
//
//  TAYLOR SERIES FOR P/R
//
    apn = *a+1.0e0;
    t = *x/apn;
    wk[0] = t;
    for ( n = 2; n <= 20; n++ )
    {
        apn += 1.0e0;
        t *= (*x/apn);
        if(t <= 1.e-3) goto S70;
        wk[n-1] = t;
    }
    n = 20;
S70:
    sum = t;
    tol = 0.5e0*acc;
S80:
    apn += 1.0e0;
    t *= (*x/apn);
    sum += t;
    if(t > tol) goto S80;
    max = n-1;
    for ( m = 1; m <= max; m++ )
    {
        n -= 1;
        sum += wk[n-1];
    }
    *ans = r/ *a*(1.0e0+sum);
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S100:
//
//  ASYMPTOTIC EXPANSION
//
    amn = *a-1.0e0;
    t = amn/ *x;
    wk[0] = t;
    for ( n = 2; n <= 20; n++ )
    {
        amn -= 1.0e0;
        t *= (amn/ *x);
        if(fabs(t) <= 1.e-3) goto S120;
        wk[n-1] = t;
    }
    n = 20;
S120:
    sum = t;
S130:
    if(fabs(t) <= acc) goto S140;
    amn -= 1.0e0;
    t *= (amn/ *x);
    sum += t;
    goto S130;
S140:
    max = n-1;
    for ( m = 1; m <= max; m++ )
    {
        n -= 1;
        sum += wk[n-1];
    }
    *qans = r/ *x*(1.0e0+sum);
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S160:
//
//  TAYLOR SERIES FOR P(A,X)/X**A
//
    an = 3.0e0;
    c = *x;
    sum = *x/(*a+3.0e0);
    tol = 3.0e0*acc/(*a+1.0e0);
S170:
    an += 1.0e0;
    c = -(c*(*x/an));
    t = c/(*a+an);
    sum += t;
    if(fabs(t) > tol) goto S170;
    j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
    z = *a*log(*x);
    h = cdflib_gam1(a);
    g = 1.0e0+h;
    if(*x < 0.25e0) goto S180;
    if(*a < *x/2.59e0) goto S200;
    goto S190;
S180:
    if(z > -.13394e0) goto S200;
S190:
    w = exp(z);
    *ans = w*g*(0.5e0+(0.5e0-j));
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S200:
    l = cdflib_rexp(&z);
    w = 0.5e0+(0.5e0+l);
    *qans = (w*j-l)*g-h;
    if(*qans < 0.0e0) goto S380;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S210:
//
//  FINITE SUMS FOR Q WHEN A .GE. 1 AND 2*A IS AN INTEGER
//
    sum = exp(-*x);
    t = sum;
    n = 1;
    c = 0.0e0;
    goto S230;
S220:
    rtx = sqrt(*x);
    sum = cdflib_error_fc ( &K2, &rtx );
    t = exp(-*x)/(rtpi*rtx);
    n = 0;
    c = -0.5e0;
S230:
    if(n == i) goto S240;
    n += 1;
    c += 1.0e0;
    t = *x*t/c;
    sum += t;
    goto S230;
S240:
    *qans = sum;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S250:
//
//  CONTINUED FRACTION EXPANSION
//
    tol = cdflib_fifdmax1(5.0e0*e,acc);
    a2nm1 = a2n = 1.0e0;
    b2nm1 = *x;
    b2n = *x+(1.0e0-*a);
    c = 1.0e0;
S260:
    a2nm1 = *x*a2n+c*a2nm1;
    b2nm1 = *x*b2n+c*b2nm1;
    am0 = a2nm1/b2nm1;
    c += 1.0e0;
    cma = c-*a;
    a2n = a2nm1+cma*a2n;
    b2n = b2nm1+cma*b2n;
    an0 = a2n/b2n;
    if(fabs(an0-am0) >= tol*an0) goto S260;
    *qans = r*an0;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S270:
//
//  GENERAL TEMME EXPANSION
//
    if(fabs(s) <= 2.0e0*e && *a*e*e > 3.28e-3) goto S430;
    c = exp(-y);
    T3 = sqrt(y);
    w = 0.5e0 * cdflib_error_fc ( &K1, &T3 );
    u = 1.0e0/ *a;
    z = sqrt(z+z);
    if(l < 1.0e0) z = -z;
    T4 = iop-2;
    if(T4 < 0) goto S280;
    else if(T4 == 0) goto S290;
    else  goto S300;
S280:
    if(fabs(s) <= 1.e-3) goto S340;
    c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
      6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
    c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
      )*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
      d2[2])*z+d2[1])*z+d2[0])*z+d20;
    c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
      d3[0])*z+d30;
    c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
    c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
    c6 = (d6[1]*z+d6[0])*z+d60;
    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
    goto S310;
S290:
    c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
    c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = d2[0]*z+d20;
    t = (c2*u+c1)*u+c0;
    goto S310;
S300:
    t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
S310:
    if(l < 1.0e0) goto S320;
    *qans = c*(w+rt2pin*t/rta);
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S320:
    *ans = c*(w-rt2pin*t/rta);
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S330:
//
//  TEMME EXPANSION FOR L = 1
//
    if(*a*e*e > 3.28e-3) goto S430;
    c = 0.5e0+(0.5e0-y);
    w = (0.5e0-sqrt(y)*(0.5e0+(0.5e0-y/3.0e0))/rtpi)/c;
    u = 1.0e0/ *a;
    z = sqrt(z+z);
    if(l < 1.0e0) z = -z;
    T5 = iop-2;
    if(T5 < 0) goto S340;
    else if(T5 == 0) goto S350;
    else  goto S360;
S340:
    c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
      third;
    c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
    c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
    c4 = (d4[1]*z+d4[0])*z+d40;
    c5 = (d5[1]*z+d5[0])*z+d50;
    c6 = d6[0]*z+d60;
    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
    goto S310;
S350:
    c0 = (d0[1]*z+d0[0])*z-third;
    c1 = d1[0]*z+d10;
    t = (d20*u+c1)*u+c0;
    goto S310;
S360:
    t = d0[0]*z-third;
    goto S310;
S370:
//
//  SPECIAL CASES
//
    *ans = 0.0e0;
    *qans = 1.0e0;
    return;
S380:
    *ans = 1.0e0;
    *qans = 0.0e0;
    return;
S390:
    if(*x >= 0.25e0) goto S400;
    T6 = sqrt(*x);
    *ans = cdflib_error_f ( &T6 );
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S400:
    T7 = sqrt(*x);
    *qans = cdflib_error_fc ( &K2, &T7 );
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S410:
    if(fabs(s) <= 2.0e0*e) goto S430;
S420:
    if(*x <= *a) goto S370;
    goto S380;
S430:
//
//  ERROR RETURN
//
    *ans = 2.0e0;
    return;
}

void
cdflib_gamma_inc_inv(
    double* a,
    double* x,
    double* x0,
    double* p,
    double* q,
    int* ierr)
{
  static const double a0 = 3.31125922108741e0;
  static const double a1 = 11.6616720288968e0;
  static const double a2 = 4.28342155967104e0;
  static const double a3 = .213623493715853e0;
  static const double b1 = 6.61053765625462e0;
  static const double b2 = 6.40691597760039e0;
  static const double b3 = 1.27364489782223e0;
  static const double b4 = .036117081018842e0;
  static const double c = .577215664901533e0;
  static const double ln10 = 2.302585e0;
  static const double tol = 1.e-5;
  static const double amin[2] = {
    500.0e0,100.0e0
  };
  static const double bmin[2] = {
    1.e-28,1.e-13
  };
  static const double dmin[2] = {
    1.e-06,1.e-04
  };
  static const double emin[2] = {
    2.e-03,6.e-03
  };
  static const double eps0[2] = {
    1.e-10,1.e-08
  };
  static int K1 = 1;
  static int K2 = 2;
  static int K3 = 3;
  static int K8 = 0;
  static double am1,amax,ap1,ap2,ap3,apn,b,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,
    r,rta,s,s2,sum,t,u,w,xmax,xmin,xn,y,z;
  static int iop;
  static double T4,T5,T6,T7,T9;

//
//  E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
//            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
//            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
//            LARGEST POSITIVE NUMBER.
//
    e = cdflib_dpmpar(&K1);
    xmin = cdflib_dpmpar(&K2);
    xmax = cdflib_dpmpar(&K3);
    *x = 0.0e0;
    if(*a <= 0.0e0) goto S300;
    t = *p+*q-1.e0;
    if(fabs(t) > e) goto S320;
    *ierr = 0;
    if(*p == 0.0e0) return;
    if(*q == 0.0e0) goto S270;
    if(*a == 1.0e0) goto S280;
    e2 = 2.0e0*e;
    amax = 0.4e-10/(e*e);
    iop = 1;
    if(e > 1.e-10) iop = 2;
    eps = eps0[iop-1];
    xn = *x0;
    if(*x0 > 0.0e0) goto S160;
//
//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
//                       WHEN A .LT. 1
//
    if(*a > 1.0e0) goto S80;
    T4 = *a+1.0e0;
    g = cdflib_gamma_x(&T4);
    qg = *q*g;
    if(qg == 0.0e0) goto S360;
    b = qg/ *a;
    if(qg > 0.6e0**a) goto S40;
    if(*a >= 0.30e0 || b < 0.35e0) goto S10;
    t = exp(-(b+c));
    u = t*exp(t);
    xn = t*exp(u);
    goto S160;
S10:
    if(b >= 0.45e0) goto S40;
    if(b == 0.0e0) goto S360;
    y = -log(b);
    s = 0.5e0+(0.5e0-*a);
    z = log(y);
    t = y-s*z;
    if(b < 0.15e0) goto S20;
    xn = y-s*log(t)-log(1.0e0+s/(t+1.0e0));
    goto S220;
S20:
    if(b <= 0.01e0) goto S30;
    u = ((t+2.0e0*(3.0e0-*a))*t+(2.0e0-*a)*(3.0e0-*a))/((t+(5.0e0-*a))*t+2.0e0);
    xn = y-s*log(t)-log(u);
    goto S220;
S30:
    c1 = -(s*z);
    c2 = -(s*(1.0e0+c1));
    c3 = s*((0.5e0*c1+(2.0e0-*a))*c1+(2.5e0-1.5e0**a));
    c4 = -(s*(((c1/3.0e0+(2.5e0-1.5e0**a))*c1+((*a-6.0e0)**a+7.0e0))*c1+(
      (11.0e0**a-46.0)**a+47.0e0)/6.0e0));
    c5 = -(s*((((-(c1/4.0e0)+(11.0e0**a-17.0e0)/6.0e0)*c1+((-(3.0e0**a)+13.0e0)*
      *a-13.0e0))*c1+0.5e0*(((2.0e0**a-25.0e0)**a+72.0e0)**a-61.0e0))*c1+((
      (25.0e0**a-195.0e0)**a+477.0e0)**a-379.0e0)/12.0e0));
    xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
    if(*a > 1.0e0) goto S220;
    if(b > bmin[iop-1]) goto S220;
    *x = xn;
    return;
S40:
    if(b**q > 1.e-8) goto S50;
    xn = exp(-(*q/ *a+c));
    goto S70;
S50:
    if(*p <= 0.9e0) goto S60;
    T5 = -*q;
    xn = exp((cdflib_alnrel(&T5)+ cdflib_gamma_ln1 ( a ) ) / *a );
    goto S70;
S60:
    xn = exp(log(*p*g)/ *a);
S70:
    if(xn == 0.0e0) goto S310;
    t = 0.5e0+(0.5e0-xn/(*a+1.0e0));
    xn /= t;
    goto S160;
S80:
//
//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
//                       WHEN A .GT. 1
//
    if(*q <= 0.5e0) goto S90;
    w = log(*p);
    goto S100;
S90:
    w = log(*q);
S100:
    t = sqrt(-(2.0e0*w));
    s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.0e0);
    if(*q > 0.5e0) s = -s;
    rta = sqrt(*a);
    s2 = s*s;
    xn = *a+s*rta+(s2-1.0e0)/3.0e0+s*(s2-7.0e0)/(36.0e0*rta)-((3.0e0*s2+7.0e0)*
      s2-16.0e0)/(810.0e0**a)+s*((9.0e0*s2+256.0e0)*s2-433.0e0)/(38880.0e0**a*
      rta);
    xn = cdflib_fifdmax1(xn,0.0e0);
    if(*a < amin[iop-1]) goto S110;
    *x = xn;
    d = 0.5e0+(0.5e0-*x/ *a);
    if(fabs(d) <= dmin[iop-1]) return;
S110:
    if(*p <= 0.5e0) goto S130;
    if(xn < 3.0e0**a) goto S220;
    y = -(w+ cdflib_gamma_log ( a ) );
    d = cdflib_fifdmax1(2.0e0,*a*(*a-1.0e0));
    if(y < ln10*d) goto S120;
    s = 1.0e0-*a;
    z = log(y);
    goto S30;
S120:
    t = *a-1.0e0;
    T6 = -(t/(xn+1.0e0));
    xn = y+t*log(xn)-cdflib_alnrel(&T6);
    T7 = -(t/(xn+1.0e0));
    xn = y+t*log(xn)-cdflib_alnrel(&T7);
    goto S220;
S130:
    ap1 = *a+1.0e0;
    if(xn > 0.70e0*ap1) goto S170;
    w += cdflib_gamma_log ( &ap1 );
    if(xn > 0.15e0*ap1) goto S140;
    ap2 = *a+2.0e0;
    ap3 = *a+3.0e0;
    *x = exp((w+*x)/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2*(1.0e0+*x/ap3))))/ *a);
    xn = *x;
    if(xn > 1.e-2*ap1) goto S140;
    if(xn <= emin[iop-1]*ap1) return;
    goto S170;
S140:
    apn = ap1;
    t = xn/apn;
    sum = 1.0e0+t;
S150:
    apn += 1.0e0;
    t *= (xn/apn);
    sum += t;
    if(t > 1.e-4) goto S150;
    t = w-log(sum);
    xn = exp((xn+t)/ *a);
    xn *= (1.0e0-(*a*log(xn)-xn-t)/(*a-xn));
    goto S170;
S160:
//
//                 SCHRODER ITERATION USING P
//
    if(*p > 0.5e0) goto S220;
S170:
    if(*p <= 1.e10*xmin) goto S350;
    am1 = *a-0.5e0-0.5e0;
S180:
    if(*a <= amax) goto S190;
    d = 0.5e0+(0.5e0-xn/ *a);
    if(fabs(d) <= e2) goto S350;
S190:
    if(*ierr >= 20) goto S330;
    *ierr += 1;
    cdflib_gamma_inc ( a, &xn, &pn, &qn, &K8 );
    if(pn == 0.0e0 || qn == 0.0e0) goto S350;
    r = cdflib_rcomp(a,&xn);
    if(r == 0.0e0) goto S350;
    t = (pn-*p)/r;
    w = 0.5e0*(am1-xn);
    if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S200;
    *x = xn*(1.0e0-t);
    if(*x <= 0.0e0) goto S340;
    d = fabs(t);
    goto S210;
S200:
    h = t*(1.0e0+w*t);
    *x = xn*(1.0e0-h);
    if(*x <= 0.0e0) goto S340;
    if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
    d = fabs(h);
S210:
    xn = *x;
    if(d > tol) goto S180;
    if(d <= eps) return;
    if(fabs(*p-pn) <= tol**p) return;
    goto S180;
S220:
//
//                 SCHRODER ITERATION USING Q
//
    if(*q <= 1.e10*xmin) goto S350;
    am1 = *a-0.5e0-0.5e0;
S230:
    if(*a <= amax) goto S240;
    d = 0.5e0+(0.5e0-xn/ *a);
    if(fabs(d) <= e2) goto S350;
S240:
    if(*ierr >= 20) goto S330;
    *ierr += 1;
    cdflib_gamma_inc ( a, &xn, &pn, &qn, &K8 );
    if(pn == 0.0e0 || qn == 0.0e0) goto S350;
    r = cdflib_rcomp(a,&xn);
    if(r == 0.0e0) goto S350;
    t = (*q-qn)/r;
    w = 0.5e0*(am1-xn);
    if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S250;
    *x = xn*(1.0e0-t);
    if(*x <= 0.0e0) goto S340;
    d = fabs(t);
    goto S260;
S250:
    h = t*(1.0e0+w*t);
    *x = xn*(1.0e0-h);
    if(*x <= 0.0e0) goto S340;
    if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
    d = fabs(h);
S260:
    xn = *x;
    if(d > tol) goto S230;
    if(d <= eps) return;
    if(fabs(*q-qn) <= tol**q) return;
    goto S230;
S270:
//
//                       SPECIAL CASES
//
    *x = xmax;
    return;
S280:
    if(*q < 0.9e0) goto S290;
    T9 = -*p;
    *x = -cdflib_alnrel(&T9);
    return;
S290:
    *x = -log(*q);
    return;
S300:
//
//                       ERROR RETURN
//
    *ierr = -2;
    return;
S310:
    *ierr = -3;
    return;
S320:
    *ierr = -4;
    return;
S330:
    *ierr = -6;
    return;
S340:
    *ierr = -7;
    return;
S350:
    *x = xn;
    *ierr = -8;
    return;
S360:
    *x = xmax;
    *ierr = -8;
    return;
}

void
cdflib_gamma_inc_values(
    int* n_data,
    double* a,
    double* x,
    double* fx)
{
  enum { N_MAX = 20 };

  static const double a_vec[N_MAX] = {
    0.1E+00,  0.1E+00,  0.1E+00,  0.5E+00,
    0.5E+00,  0.5E+00,  1.0E+00,  1.0E+00,
    1.0E+00,  1.1E+00,  1.1E+00,  1.1E+00,
    2.0E+00,  2.0E+00,  2.0E+00,  6.0E+00,
    6.0E+00, 11.0E+00, 26.0E+00, 41.0E+00 };
  static const double fx_vec[N_MAX] = {
    0.7420263E+00, 0.9119753E+00, 0.9898955E+00, 0.2931279E+00,
    0.7656418E+00, 0.9921661E+00, 0.0951626E+00, 0.6321206E+00,
    0.9932621E+00, 0.0757471E+00, 0.6076457E+00, 0.9933425E+00,
    0.0091054E+00, 0.4130643E+00, 0.9931450E+00, 0.0387318E+00,
    0.9825937E+00, 0.9404267E+00, 0.4863866E+00, 0.7359709E+00 };
  static const double x_vec[N_MAX] = {
    3.1622777E-02, 3.1622777E-01, 1.5811388E+00, 7.0710678E-02,
    7.0710678E-01, 3.5355339E+00, 0.1000000E+00, 1.0000000E+00,
    5.0000000E+00, 1.0488088E-01, 1.0488088E+00, 5.2440442E+00,
    1.4142136E-01, 1.4142136E+00, 7.0710678E+00, 2.4494897E+00,
    1.2247449E+01, 1.6583124E+01, 2.5495098E+01, 4.4821870E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_gamma_ln1(
    double* a)
{
  static const double p0 = .577215664901533e+00;
  static const double p1 = .844203922187225e+00;
  static const double p2 = -.168860593646662e+00;
  static const double p3 = -.780427615533591e+00;
  static const double p4 = -.402055799310489e+00;
  static const double p5 = -.673562214325671e-01;
  static const double p6 = -.271935708322958e-02;
  static const double q1 = .288743195473681e+01;
  static const double q2 = .312755088914843e+01;
  static const double q3 = .156875193295039e+01;
  static const double q4 = .361951990101499e+00;
  static const double q5 = .325038868253937e-01;
  static const double q6 = .667465618796164e-03;
  static const double r0 = .422784335098467e+00;
  static const double r1 = .848044614534529e+00;
  static const double r2 = .565221050691933e+00;
  static const double r3 = .156513060486551e+00;
  static const double r4 = .170502484022650e-01;
  static const double r5 = .497958207639485e-03;
  static const double s1 = .124313399877507e+01;
  static const double s2 = .548042109832463e+00;
  static const double s3 = .101552187439830e+00;
  static const double s4 = .713309612391000e-02;
  static const double s5 = .116165475989616e-03;
  static double gamln1,w,x;

    if(*a >= 0.6e0) goto S10;
    w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/((((((q6**a+q5)**a+
      q4)**a+q3)**a+q2)**a+q1)**a+1.0e0);
    gamln1 = -(*a*w);
    return gamln1;
S10:
    x = *a-0.5e0-0.5e0;
    w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
      +1.0e0);
    gamln1 = x*w;
    return gamln1;
}

double
cdflib_gamma_log(
    double* a)
{
  static const double c0 = .833333333333333e-01;
  static const double c1 = -.277777777760991e-02;
  static const double c2 = .793650666825390e-03;
  static const double c3 = -.595202931351870e-03;
  static const double c4 = .837308034031215e-03;
  static const double c5 = -.165322962780713e-02;
  static const double d = .418938533204673e0;
  static double gamln,t,w;
  static int i,n;
  static double T1;

    if(*a > 0.8e0) goto S10;
    gamln = cdflib_gamma_ln1 ( a ) - log ( *a );
    return gamln;
S10:
    if(*a > 2.25e0) goto S20;
    t = *a-0.5e0-0.5e0;
    gamln = cdflib_gamma_ln1 ( &t );
    return gamln;
S20:
    if(*a >= 10.0e0) goto S40;
    n = ( int ) ( *a - 1.25e0 );
    t = *a;
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        t -= 1.0e0;
        w = t*w;
    }
    T1 = t-1.0e0;
    gamln = cdflib_gamma_ln1 ( &T1 ) + log ( w );
    return gamln;
S40:
    t = pow(1.0e0/ *a,2.0);
    w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
    gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
    return gamln;
}

void
cdflib_gamma_rat1(
    double* a,
    double* x,
    double* r,
    double* p,
    double* q,
    double* eps)
{
  static int K2 = 0;
  static double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;

    if(*a**x == 0.0e0) goto S120;
    if(*a == 0.5e0) goto S100;
    if(*x < 1.1e0) goto S10;
    goto S60;
S10:
//
//             TAYLOR SERIES FOR P(A,X)/X**A
//
    an = 3.0e0;
    c = *x;
    sum = *x/(*a+3.0e0);
    tol = 0.1e0**eps/(*a+1.0e0);
S20:
    an += 1.0e0;
    c = -(c*(*x/an));
    t = c/(*a+an);
    sum += t;
    if(fabs(t) > tol) goto S20;
    j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
    z = *a*log(*x);
    h = cdflib_gam1(a);
    g = 1.0e0+h;
    if(*x < 0.25e0) goto S30;
    if(*a < *x/2.59e0) goto S50;
    goto S40;
S30:
    if(z > -.13394e0) goto S50;
S40:
    w = exp(z);
    *p = w*g*(0.5e0+(0.5e0-j));
    *q = 0.5e0+(0.5e0-*p);
    return;
S50:
    l = cdflib_rexp(&z);
    w = 0.5e0+(0.5e0+l);
    *q = (w*j-l)*g-h;
    if(*q < 0.0e0) goto S90;
    *p = 0.5e0+(0.5e0-*q);
    return;
S60:
//
//              CONTINUED FRACTION EXPANSION
//
    a2nm1 = a2n = 1.0e0;
    b2nm1 = *x;
    b2n = *x+(1.0e0-*a);
    c = 1.0e0;
S70:
    a2nm1 = *x*a2n+c*a2nm1;
    b2nm1 = *x*b2n+c*b2nm1;
    am0 = a2nm1/b2nm1;
    c += 1.0e0;
    cma = c-*a;
    a2n = a2nm1+cma*a2n;
    b2n = b2nm1+cma*b2n;
    an0 = a2n/b2n;
    if(fabs(an0-am0) >= *eps*an0) goto S70;
    *q = *r*an0;
    *p = 0.5e0+(0.5e0-*q);
    return;
S80:
//
//                SPECIAL CASES
//
    *p = 0.0e0;
    *q = 1.0e0;
    return;
S90:
    *p = 1.0e0;
    *q = 0.0e0;
    return;
S100:
    if(*x >= 0.25e0) goto S110;
    T1 = sqrt(*x);
    *p = cdflib_error_f ( &T1 );
    *q = 0.5e0+(0.5e0-*p);
    return;
S110:
    T3 = sqrt(*x);
    *q = cdflib_error_fc ( &K2, &T3 );
    *p = 0.5e0+(0.5e0-*q);
    return;
S120:
    if(*x <= *a) goto S80;
    goto S90;
}

void
cdflib_gamma_values(
    int* n_data,
    double* x,
    double* fx)
{
  enum { N_MAX = 18 };

  static const double fx_vec[N_MAX] = {
    4.590845E+00,     2.218160E+00,     1.489192E+00,     1.164230E+00,
    1.0000000000E+00, 0.9513507699E+00, 0.9181687424E+00, 0.8974706963E+00,
    0.8872638175E+00, 0.8862269255E+00, 0.8935153493E+00, 0.9086387329E+00,
    0.9313837710E+00, 0.9617658319E+00, 1.0000000000E+00, 3.6288000E+05,
    1.2164510E+17,    8.8417620E+30 };
  static const double x_vec[N_MAX] = {
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,
    1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00,
    1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00,
    1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00,
   20.0E+00, 30.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_gamma_x(
    double* a)
{
  static const double d = .41893853320467274178e0;
  static const double pi = 3.1415926535898e0;
  static const double r1 = .820756370353826e-03;
  static const double r2 = -.595156336428591e-03;
  static const double r3 = .793650663183693e-03;
  static const double r4 = -.277777777770481e-02;
  static const double r5 = .833333333333333e-01;
  static const double p[7] = {
    .539637273585445e-03,.261939260042690e-02,.204493667594920e-01,
    .730981088720487e-01,.279648642639792e+00,.553413866010467e+00,1.0e0
  };
  static const double q[7] = {
    -.832979206704073e-03,.470059485860584e-02,.225211131035340e-01,
    -.170458969313360e+00,-.567902761974940e-01,.113062953091122e+01,1.0e0
  };
  static int K2 = 3;
  static int K3 = 0;
  static double Xgamm,bot,g,lnx,s,t,top,w,x,z;
  static int i,j,m,n,T1;

    Xgamm = 0.0e0;
    x = *a;
    if(fabs(*a) >= 15.0e0) goto S110;
//
//            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
//
    t = 1.0e0;
    m = cdflib_fifidint(*a)-1;
//
//     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
//
    T1 = m;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S30;
    else  goto S10;
S10:
    for ( j = 1; j <= m; j++ )
    {
        x -= 1.0e0;
        t = x*t;
    }
S30:
    x -= 1.0e0;
    goto S80;
S40:
//
//     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
//
    t = *a;
    if(*a > 0.0e0) goto S70;
    m = -m-1;
    if(m == 0) goto S60;
    for ( j = 1; j <= m; j++ )
    {
        x += 1.0e0;
        t = x*t;
    }
S60:
    x += (0.5e0+0.5e0);
    t = x*t;
    if(t == 0.0e0) return Xgamm;
S70:
//
//     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
//     CODE MAY BE OMITTED IF DESIRED.
//
    if(fabs(t) >= 1.e-30) goto S80;
    if(fabs(t)*cdflib_dpmpar(&K2) <= 1.0001e0) return Xgamm;
    Xgamm = 1.0e0/t;
    return Xgamm;
S80:
//
//     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
//
    top = p[0];
    bot = q[0];
    for ( i = 1; i < 7; i++ )
    {
        top = p[i]+x*top;
        bot = q[i]+x*bot;
    }
    Xgamm = top/bot;
//
//     TERMINATION
//
    if(*a < 1.0e0) goto S100;
    Xgamm *= t;
    return Xgamm;
S100:
    Xgamm /= t;
    return Xgamm;
S110:
//
//  EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
//
    if(fabs(*a) >= 1.e3) return Xgamm;
    if(*a > 0.0e0) goto S120;
    x = -*a;
    n = ( int ) x;
    t = x-(double)n;
    if(t > 0.9e0) t = 1.0e0-t;
    s = sin(pi*t)/pi;
    if(cdflib_fifmod(n,2) == 0) s = -s;
    if(s == 0.0e0) return Xgamm;
S120:
//
//     COMPUTE THE MODIFIED ASYMPTOTIC SUM
//
    t = 1.0e0/(x*x);
    g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x;
//
//     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
//     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
//
    lnx = log(x);
//
//  FINAL ASSEMBLY
//
    z = x;
    g = d+g+(z-0.5e0)*(lnx-1.e0);
    w = g;
    t = g-w;
    if(w > 0.99999e0*cdflib_exparg(&K3)) return Xgamm;
    Xgamm = exp(w)*(1.0e0+t);
    if(*a < 0.0e0) Xgamm = 1.0e0/(Xgamm*s)/x;
    return Xgamm;
}

double
cdflib_gsumln(
    double* a,
    double* b)
{
  static double gsumln,x,T1,T2;

    x = *a+*b-2.e0;
    if(x > 0.25e0) goto S10;
    T1 = 1.0e0+x;
    gsumln = cdflib_gamma_ln1 ( &T1 );
    return gsumln;
S10:
    if(x > 1.25e0) goto S20;
    gsumln = cdflib_gamma_ln1 ( &x ) + cdflib_alnrel ( &x );
    return gsumln;
S20:
    T2 = x-1.0e0;
    gsumln = cdflib_gamma_ln1 ( &T2 ) + log ( x * ( 1.0e0 + x ) );
    return gsumln;
}

int
cdflib_ipmpar(
    int* i)
{
   // MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
   //   3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
   //   PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
   static const int imach[11] = {
       /* imach[ 0] = */ 0,
       /* imach[ 1] = */ 2,
       /* imach[ 2] = */ 31,
       /* imach[ 3] = */ 2147483647,
       /* imach[ 4] = */ 2,
       /* imach[ 5] = */ 24,
       /* imach[ 6] = */-125,
       /* imach[ 7] = */ 128,
       /* imach[ 8] = */ 53,
       /* imach[ 9] = */-1021,
       /* imach[10] = */ 1024,
   };

   return imach[*i];
}

void
cdflib_negative_binomial_cdf_values(
    int* n_data,
    int* f,
    int* s,
    double* p,
    double* cdf)
{
  enum { N_MAX = 27 };

  static const double cdf_vec[N_MAX] = {
    0.6367, 0.3633, 0.1445,
    0.5000, 0.2266, 0.0625,
    0.3438, 0.1094, 0.0156,
    0.1792, 0.0410, 0.0041,
    0.0705, 0.0109, 0.0007,
    0.9862, 0.9150, 0.7472,
    0.8499, 0.5497, 0.2662,
    0.6513, 0.2639, 0.0702,
    1.0000, 0.0199, 0.0001 };
  static const int f_vec[N_MAX] = {
     4,  3,  2,
     3,  2,  1,
     2,  1,  0,
     2,  1,  0,
     2,  1,  0,
    11, 10,  9,
    17, 16, 15,
     9,  8,  7,
     2,  1,  0 };
  static const double p_vec[N_MAX] = {
    0.50, 0.50, 0.50,
    0.50, 0.50, 0.50,
    0.50, 0.50, 0.50,
    0.40, 0.40, 0.40,
    0.30, 0.30, 0.30,
    0.30, 0.30, 0.30,
    0.10, 0.10, 0.10,
    0.10, 0.10, 0.10,
    0.01, 0.01, 0.01 };
  static const int s_vec[N_MAX] = {
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    0, 1, 2 };

  if ( n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *f = 0;
    *s = 0;
    *p = 0.0E+00;
    *cdf = 0.0E+00;
  }
  else
  {
    *f = f_vec[*n_data-1];
    *s = s_vec[*n_data-1];
    *p = p_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
}

void
cdflib_normal_cdf_values(
    int* n_data,
    double* x,
    double* fx)
{
  enum { N_MAX = 13 };

  static const double fx_vec[N_MAX] = {
    0.500000000000000E+00, 0.539827837277029E+00, 0.579259709439103E+00,
    0.617911422188953E+00, 0.655421741610324E+00, 0.691462461274013E+00,
    0.725746882249927E+00, 0.758036347776927E+00, 0.788144601416604E+00,
    0.815939874653241E+00, 0.841344746068543E+00, 0.933192798731142E+00,
    0.977249868051821E+00 };
  static const double x_vec[N_MAX] = {
    0.00E+00, 0.10E+00, 0.20E+00,
    0.30E+00, 0.40E+00, 0.50E+00,
    0.60E+00, 0.70E+00, 0.80E+00,
    0.90E+00, 1.00E+00, 1.50E+00,
    2.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
}

void
cdflib_poisson_cdf_values(
    int* n_data,
    double* a,
    int* x,
    double* fx)
{
  enum { N_MAX = 21 };

  static const double a_vec[N_MAX] = {
    0.02E+00, 0.10E+00, 0.10E+00, 0.50E+00,
    0.50E+00, 0.50E+00, 1.00E+00, 1.00E+00,
    1.00E+00, 1.00E+00, 2.00E+00, 2.00E+00,
    2.00E+00, 2.00E+00, 5.00E+00, 5.00E+00,
    5.00E+00, 5.00E+00, 5.00E+00, 5.00E+00,
    5.00E+00 };
  static const double fx_vec[N_MAX] = {
    0.980E+00, 0.905E+00, 0.995E+00, 0.607E+00,
    0.910E+00, 0.986E+00, 0.368E+00, 0.736E+00,
    0.920E+00, 0.981E+00, 0.135E+00, 0.406E+00,
    0.677E+00, 0.857E+00, 0.007E+00, 0.040E+00,
    0.125E+00, 0.265E+00, 0.441E+00, 0.616E+00,
    0.762E+00 };
  static const int x_vec[N_MAX] = {
     0, 0, 1, 0,
     1, 2, 0, 1,
     2, 3, 0, 1,
     2, 3, 0, 1,
     2, 3, 4, 5,
     6 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *x = 0;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_psi(
    double* xx)
{
  static const double dx0 = 1.461632144968362341262659542325721325e0;
  static const double piov4 = .785398163397448e0;
  static const double p1[7] = {
    .895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
    .118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
    .130560269827897e+04
  };
  static const double p2[4] = {
    -.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
    -.648157123766197e+00
  };
  static const double q1[6] = {
    .448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
    .364127349079381e+04,.190831076596300e+04,.691091682714533e-05
  };
  static const double q2[4] = {
    .322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
    .777788548522962e+01
  };
  static int K1 = 3;
  static int K2 = 1;
  static double psi,aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z;
  static int i,m,n,nq;
//
//     MACHINE DEPENDENT CONSTANTS ...
//        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
//                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
//                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
//                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
//                 PSI MAY BE REPRESENTED AS ALOG(X).
//        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
//                 MAY BE REPRESENTED BY 1/X.
//
    xmax1 = cdflib_ipmpar(&K1);
    xmax1 = cdflib_fifdmin1(xmax1,1.0e0/cdflib_dpmpar(&K2));
    xsmall = 1.e-9;
    x = *xx;
    aug = 0.0e0;
    if(x >= 0.5e0) goto S50;
//
//     X .LT. 0.5,  USE REFLECTION FORMULA
//     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
//
    if(fabs(x) > xsmall) goto S10;
    if(x == 0.0e0) goto S100;
//
//     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
//     FOR  PI*COTAN(PI*X)
//
    aug = -(1.0e0/x);
    goto S40;
S10:
//
//     REDUCTION OF ARGUMENT FOR COTAN
//
    w = -x;
    sgn = piov4;
    if(w > 0.0e0) goto S20;
    w = -w;
    sgn = -sgn;
S20:
//
//     MAKE AN ERROR EXIT IF X .LE. -XMAX1
//
    if(w >= xmax1) goto S100;
    nq = cdflib_fifidint(w);
    w -= (double)nq;
    nq = cdflib_fifidint(w*4.0e0);
    w = 4.0e0*(w-(double)nq*.25e0);
//
//     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
//     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
//     QUADRANT AND DETERMINE SIGN
//
    n = nq/2;
    if(n+n != nq) w = 1.0e0-w;
    z = piov4*w;
    m = n/2;
    if(m+m != n) sgn = -sgn;
//
//     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
//
    n = (nq+1)/2;
    m = n/2;
    m += m;
    if(m != n) goto S30;
//
//     CHECK FOR SINGULARITY
//
    if(z == 0.0e0) goto S100;
//
//     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
//     SIN/COS AS A SUBSTITUTE FOR TAN
//
    aug = sgn*(cos(z)/sin(z)*4.0e0);
    goto S40;
S30:
    aug = sgn*(sin(z)/cos(z)*4.0e0);
S40:
    x = 1.0e0-x;
S50:
    if(x > 3.0e0) goto S70;
//
//     0.5 .LE. X .LE. 3.0
//
    den = x;
    upper = p1[0]*x;
    for ( i = 1; i <= 5; i++ )
    {
        den = (den+q1[i-1])*x;
        upper = (upper+p1[i+1-1])*x;
    }
    den = (upper+p1[6])/(den+q1[5]);
    xmx0 = x-dx0;
    psi = den*xmx0+aug;
    return psi;
S70:
//
//     IF X .GE. XMAX1, PSI = LN(X)
//
    if(x >= xmax1) goto S90;
//
//     3.0 .LT. X .LT. XMAX1
//
    w = 1.0e0/(x*x);
    den = w;
    upper = p2[0]*w;
    for ( i = 1; i <= 3; i++ )
    {
        den = (den+q2[i-1])*w;
        upper = (upper+p2[i+1-1])*w;
    }
    aug = upper/(den+q2[3])-0.5e0/x+aug;
S90:
    psi = aug+log(x);
    return psi;
S100:
//
//     ERROR RETURN
//
    psi = 0.0e0;
    return psi;
}

void
cdflib_psi_values(
    int* n_data,
    double* x,
    double* fx)
{
  enum { N_MAX = 11 };

  static const double fx_vec[N_MAX] = {
    -0.5772156649E+00, -0.4237549404E+00, -0.2890398966E+00,
    -0.1691908889E+00, -0.0613845446E+00, -0.0364899740E+00,
     0.1260474528E+00,  0.2085478749E+00,  0.2849914333E+00,
     0.3561841612E+00,  0.4227843351E+00 };
  static const double x_vec[N_MAX] = {
    1.0E+00,  1.1E+00,  1.2E+00,
    1.3E+00,  1.4E+00,  1.5E+00,
    1.6E+00,  1.7E+00,  1.8E+00,
    1.9E+00,  2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
}

double
cdflib_rcomp(
    double* a,
    double* x)
{
  static const double rt2pin = .398942280401433e0;
  static double rcomp,t,t1,u;
    rcomp = 0.0e0;
    if(*a >= 20.0e0) goto S20;
    t = *a*log(*x)-*x;
    if(*a >= 1.0e0) goto S10;
    rcomp = *a*exp(t)*(1.0e0+cdflib_gam1(a));
    return rcomp;
S10:
    rcomp = exp(t)/ cdflib_gamma_x(a);
    return rcomp;
S20:
    u = *x/ *a;
    if(u == 0.0e0) return rcomp;
    t = pow(1.0e0/ *a,2.0);
    t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
    t1 -= (*a*cdflib_rlog(&u));
    rcomp = rt2pin*sqrt(*a)*exp(t1);
    return rcomp;
}

double
cdflib_rexp(
    double* x)
{
  static const double p1 = .914041914819518e-09;
  static const double p2 = .238082361044469e-01;
  static const double q1 = -.499999999085958e+00;
  static const double q2 = .107141568980644e+00;
  static const double q3 = -.119041179760821e-01;
  static const double q4 = .595130811860248e-03;
  static double rexp,w;

    if(fabs(*x) > 0.15e0) goto S10;
    rexp = *x*(((p2**x+p1)**x+1.0e0)/((((q4**x+q3)**x+q2)**x+q1)**x+1.0e0));
    return rexp;
S10:
    w = exp(*x);
    if(*x > 0.0e0) goto S20;
    rexp = w-0.5e0-0.5e0;
    return rexp;
S20:
    rexp = w*(0.5e0+(0.5e0-1.0e0/w));
    return rexp;
}

double
cdflib_rlog(
    double* x)
{
  static const double a = .566749439387324e-01;
  static const double b = .456512608815524e-01;
  static const double p0 = .333333333333333e+00;
  static const double p1 = -.224696413112536e+00;
  static const double p2 = .620886815375787e-02;
  static const double q1 = -.127408923933623e+01;
  static const double q2 = .354508718369557e+00;
  static double rlog,r,t,u,w,w1;

    if(*x < 0.61e0 || *x > 1.57e0) goto S40;
    if(*x < 0.82e0) goto S10;
    if(*x > 1.18e0) goto S20;
//
//              ARGUMENT REDUCTION
//
    u = *x-0.5e0-0.5e0;
    w1 = 0.0e0;
    goto S30;
S10:
    u = *x-0.7e0;
    u /= 0.7e0;
    w1 = a-u*0.3e0;
    goto S30;
S20:
    u = 0.75e0**x-1.e0;
    w1 = b+u/3.0e0;
S30:
//
//               SERIES EXPANSION
//
    r = u/(u+2.0e0);
    t = r*r;
    w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
    rlog = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
    return rlog;
S40:
    r = *x-0.5e0-0.5e0;
    rlog = r-log(*x);
    return rlog;
}

double
cdflib_rlog1(
    double* x)
{
  static const double a = .566749439387324e-01;
  static const double b = .456512608815524e-01;
  static const double p0 = .333333333333333e+00;
  static const double p1 = -.224696413112536e+00;
  static const double p2 = .620886815375787e-02;
  static const double q1 = -.127408923933623e+01;
  static const double q2 = .354508718369557e+00;
  static double rlog1,h,r,t,w,w1;

    if(*x < -0.39e0 || *x > 0.57e0) goto S40;
    if(*x < -0.18e0) goto S10;
    if(*x > 0.18e0) goto S20;
//
//              ARGUMENT REDUCTION
//
    h = *x;
    w1 = 0.0e0;
    goto S30;
S10:
    h = *x+0.3e0;
    h /= 0.7e0;
    w1 = a-h*0.3e0;
    goto S30;
S20:
    h = 0.75e0**x-0.25e0;
    w1 = b+h/3.0e0;
S30:
//
//               SERIES EXPANSION
//
    r = h/(h+2.0e0);
    t = r*r;
    w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
    rlog1 = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
    return rlog1;
S40:
    w = *x+0.5e0+0.5e0;
    rlog1 = *x-log(w);
    return rlog1;
}

void
cdflib_student_cdf_values(
    int* n_data,
    int* a,
    double* x,
    double* fx)
{
  enum { N_MAX = 13 };

  static const int a_vec[N_MAX] = {
    1, 2, 3, 4,
    5, 2, 5, 2,
    5, 2, 3, 4,
    5 };
  static const double fx_vec[N_MAX] = {
    0.60E+00, 0.60E+00, 0.60E+00, 0.60E+00,
    0.60E+00, 0.75E+00, 0.75E+00, 0.95E+00,
    0.95E+00, 0.99E+00, 0.99E+00, 0.99E+00,
    0.99E+00 };
  static const double x_vec[N_MAX] = {
    0.325E+00, 0.289E+00, 0.277E+00, 0.271E+00,
    0.267E+00, 0.816E+00, 0.727E+00, 2.920E+00,
    2.015E+00, 6.965E+00, 4.541E+00, 3.747E+00,
    3.365E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
}

double
cdflib_stvaln(
    double* p)
{
  static const double xden[5] = {
    0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
    0.38560700634e-2
  };
  static const double xnum[5] = {
    -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
    -0.453642210148e-4
  };
  static int K1 = 5;
  static double stvaln,sign,y,z;

    if(!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
S10:
    sign = 1.0e0;
    z = 1.0e0-*p;
S20:
    y = sqrt(-(2.0e0*log(z)));
    stvaln = y+ cdflib_eval_pol ( xnum, &K1, &y ) / cdflib_eval_pol ( xden, &K1, &y );
    stvaln = sign*stvaln;
    return stvaln;
}
