#include <math.h>

int gaussg_(int, double *, double *);
int ordleg_(double *,double *, int *);

int gaussg_(int nzero, double *f, double *wt)  {

    /* System generated locals */
    int i__1;

    /* Local variables */
    static double xlim, piov2, a, b, g;
    static int i;
    static double ftemp, gtemp, fi, dn, fn, gm, gp, pi;
    static int ir;
    static double gt;
    extern /* Subroutine */ int ordleg_();
    static double fi1, dn1;
    static int irm, irp;



/* 189 "blkbox.for" */
    /* Parameter adjustments */
    --wt;
    --f;

    /* Function Body */
    xlim = 1e-12;
    ir = nzero + nzero;
    pi = atan(1.) * 4.;
    fi = (double) ir;
    fi1 = fi + 1.;
    piov2 = pi * .5;
    fn = piov2 / nzero;
    i__1 = nzero;
    for (i = 1; i <= i__1; ++i) {
/* L10: */
	wt[i] = i - .5;
    }
    i__1 = nzero;
    for (i = 1; i <= i__1; ++i) {
/* L20: */
	f[i] = sin(wt[i] * fn + piov2);
    }
    dn = fi / sqrt(fi * 4. * fi - 1.);
    dn1 = fi1 / sqrt(fi1 * 4. * fi1 - 1.);
    a = dn1 * fi;
    b = dn * fi1;
    irp = ir + 1;
    irm = ir - 1;
    i__1 = nzero;
    for (i = 1; i <= i__1; ++i) {
L5:
	ordleg_(&g, &f[i], &ir);
	ordleg_(&gm, &f[i], &irm);
	ordleg_(&gp, &f[i], &irp);
	gt = (f[i] * f[i] - 1.) / (a * gp - b * gm);
	ftemp = f[i] - g * gt;
	gtemp = f[i] - ftemp;
	f[i] = ftemp;
	if (fabs(gtemp) > xlim) {
	    goto L5;
	}
/* L2: */
    }
    i__1 = nzero;
    for (i = 1; i <= i__1; ++i) {
	a = (1. - f[i] * f[i]) * 2.;
	ordleg_(&b, &f[i], &irm);
	b = b * b * fi * fi;
	wt[i] = a * (fi - .5) / b;
/* L6: */
    }
    return 0;
} /* gaussg_ */

int ordleg_(double *sx,double *coa, int *ir)  {

    /* System generated locals */
    int i__1;

    /* Local variables */
    static int irpp;
    static double fn2sq, a, b;
    static int k, n;
    static double delta, theta, c1;
    static int irppm;
    static double c4;
    static int n1;
    static double s1, fk, fn;
    static int kk;
    static double fn2, ang, sqr2;

/* 359 "blkbox.for" */
    irpp = *ir + 1;
    irppm = irpp - 1;
    delta = acos(*coa);
    sqr2 = sqrt(2.);
/* 364 "blkbox.for" */
    theta = delta;
    c1 = sqr2;
    i__1 = irppm;
    for (n = 1; n <= i__1; ++n) {
	fn = (double) n;
	fn2 = fn + fn;
	fn2sq = fn2 * fn2;
	c1 *= sqrt(1. - 1. / fn2sq);
/* L20: */
    }
/* 373 "blkbox.for" */
    n = irppm;
    ang = fn * theta;
    s1 = 0.;
    c4 = 1.;
    a = -1.;
    b = 0.;
    n1 = n + 1;
    i__1 = n1;
    for (kk = 1; kk <= i__1; kk += 2) {
	k = kk - 1;
	if (k == n) {
	    c4 *= .5;
	}
	s1 += c4 * cos(ang);
	a += 2.;
	b += 1.;
	fk = (double) k;
	ang = theta * (fn - fk - 2.);
	c4 = a * (fn - b + 1.) / (b * (fn2 - a)) * c4;
/* L27: */
    }
    *sx = s1 * c1;
/* 392 "blkbox.for" */
    return 0;
} /* ordleg_ */

