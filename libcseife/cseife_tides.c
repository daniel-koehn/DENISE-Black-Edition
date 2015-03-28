/*! \file cseife_tides.c
 * \brief remove tides (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cseife_tides.c 3373 2010-11-14 13:29:20Z tforb $
 * \date 28/06/2005
 * 
 * remove tides (implementation)
 * 
 * Copyright 1984 by Erhard Wielandt
 * This algorithm was part of seife.f. A current version of seife.f can be
 * obtained from http://www.software-for-seismometry.de/
 * 
 * REVISIONS and CHANGES 
 *  - 28/06/2005   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_CSEIFE_TIDES_C_VERSION \
  "TF_CSEIFE_TIDES_C   V1.0   "
#define TF_CSEIFE_TIDES_C_CVSID \
  "$Id: cseife_tides.c 3373 2010-11-14 13:29:20Z tforb $"

#include <cseife.h>
#include <math.h>

/* subs/seife_tides.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

   the code was derived through f2c, but modified thereafter
*/

/* detide with synthetic tides interpolated over ni samples */
void seife_tid(double* x, int n, double dt, int nstep)
{
    /* Initialized data */

    static double omega[6] = { 1.93227,.92954,2.,1.00274,1.89567,.89293 };
    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    double tdif, omeg[6];
    int ndim;
    double step, tint, omeg0, xgez1, xgez2, xgez3, a[169]	/* 
	    was [13][13] */, c__[13], d__[13], e[13], f[13];
    int i__, j, k;
    double t;
    int nfreq, k2, k3;
    double tstep2;
    int jj;
    double rs[13], sx, dth;
    int nco;
    double cor[6], dur;

/*  remove tides. number of frequencies is automatically chosen according */
/*  to the total length of the record. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ndim = 13;
    dth = dt / two;

    if (nstep == 0) {
	nstep = (float)300. / dt;
    }
    nstep = nstep > 1 ? nstep : 1;
    step = (double) nstep;
    tstep2 = (step - one) * dth;
    tint = tstep2 + dth;
/*  determine the number of frequencies required for a good fit */
    dur = n * dt / 3600.;
    nfreq = 6;
    if (dur < 35.) {
	nfreq = 5;
    }
    if (dur < 18.) {
	nfreq = 4;
    }
    if (dur < 14.) {
	nfreq = 3;
    }
    if (dur < 5.) {
	nfreq = 2;
    }
    nco = (nfreq << 1) + 1;
    omeg0 = atan(one) * 8. / 86400.;
    i__1 = nfreq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L101: */
	omeg[i__ - 1] = omeg0 * omega[i__ - 1];
    }
/*  determine partial amplitudes by least-squares fit */
    i__1 = nco;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rs[i__ - 1] = 0.;
	i__2 = nco;
	for (k = 1; k <= i__2; ++k) {
/* L1: */
	    a[i__ + k * 13 - 14] = 0.;
	}
    }
/*  correction for averaging over nstep samples */
    i__2 = nfreq;
    for (j = 1; j <= i__2; ++j) {
/* L102: */
	cor[j - 1] = step * sin(omeg[j - 1] * dth) / sin(step * omeg[j - 1] * 
		dth);
    }
/*  set up system of linear equations */
    c__[0] = one;
    i__2 = n - nstep + 1;
    i__1 = nstep;
    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	sx = zero;
	i__3 = j + nstep - 1;
	for (jj = j; jj <= i__3; ++jj) {
/* L12: */
	    sx += x[jj];
	}
	sx /= step;
	t = (j - 1) * dt + tstep2;
	i__3 = nfreq;
	for (k = 1; k <= i__3; ++k) {
	    c__[(k << 1) - 1] = cos(omeg[k - 1] * t);
/* L103: */
	    c__[k * 2] = sin(omeg[k - 1] * t);
	}
	i__3 = nco;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    rs[i__ - 1] += sx * c__[i__ - 1];
	    i__4 = nco;
	    for (k = 1; k <= i__4; ++k) {
/* L2: */
		a[i__ + k * 13 - 14] += c__[i__ - 1] * c__[k - 1];
	    }
	}
    }
/*  solve for partial amplitudes */
    seife_gauss(a, nco, ndim, rs, f);
    i__4 = n - nstep + 1;
    i__3 = nstep;
    for (j = 1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
	t = (j - 1) * dt + tstep2;
/*  remove average and tides */
	i__1 = nfreq;
	for (k = 1; k <= i__1; ++k) {
	    k2 = k << 1;
	    k3 = k2 + 1;
	    c__[k2 - 1] = cos(omeg[k - 1] * t);
	    c__[k3 - 1] = sin(omeg[k - 1] * t);
	    d__[k2 - 1] = -omeg[k - 1] * c__[k3 - 1];
	    d__[k3 - 1] = omeg[k - 1] * c__[k2 - 1];
	    e[k2 - 1] = -omeg[k - 1] * d__[k3 - 1];
/* L104: */
	    e[k3 - 1] = omeg[k - 1] * d__[k2 - 1];
	}
	xgez1 = f[0];
	xgez2 = zero;
	xgez3 = zero;
	i__1 = nfreq;
	for (k = 1; k <= i__1; ++k) {
	    k2 = k << 1;
	    k3 = k2 + 1;
	    xgez1 += cor[k - 1] * (f[k2 - 1] * c__[k2 - 1] + f[k3 - 1] * c__[
		    k3 - 1]);
	    xgez2 += cor[k - 1] * (f[k2 - 1] * d__[k2 - 1] + f[k3 - 1] * d__[
		    k3 - 1]);
/* L105: */
	    xgez3 += cor[k - 1] * (f[k2 - 1] * e[k2 - 1] + f[k3 - 1] * e[k3 - 
		    1]) / two;
	}
	if (j > n - (nstep << 1) + 1) {
	    nstep = n + 1 - j;
	}
	i__1 = j + nstep - 1;
	for (jj = j; jj <= i__1; ++jj) {
	    tdif = (jj - j) * dt - tstep2;
/* L3: */
	    x[jj] = x[jj] - xgez1 - xgez2 * tdif - xgez3 * tdif * tdif;
	}
    }
} /* seife_tides__ */


/* ----- END OF cseife_tides.c ----- */
