/** \file
 * \copydoc as3.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "as3.h"

#include <assert.h>
#include <math.h>

// Algorithm AS 3 from Applied Statistics (1968) Volume 17 Page 189
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
