/*! \file cseife_deriv.c
 * \brief calculate time derivative (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cseife_deriv.c 3373 2010-11-14 13:29:20Z tforb $
 * \date 28/06/2005
 * 
 * calculate time derivative (implementation)
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
#define TF_CSEIFE_DERIV_C_VERSION \
  "TF_CSEIFE_DERIV_C   V1.0   "
#define TF_CSEIFE_DERIV_C_CVSID \
  "$Id: cseife_deriv.c 3373 2010-11-14 13:29:20Z tforb $"

#include <cseife.h>

/* subs/seife_deriv.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

   the code was derived through f2c, but modified thereafter
*/

/* derivative (time constant t0) */
void seife_dif(double* x, int n, double dt, double tau)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j;
    static double twodt;

/*  derivative (non-recursive, non-causal, symmetric-difference) */
    /* Parameter adjustments */
    --x;

    if (tau < 1.e-23) {
	tau = 1.;
    }
    twodt = dt * 2. / tau;
    i__1 = n - 2;
    for (j = 1; j <= i__1; ++j) {
/* L1: */
	x[j] = (x[j + 2] - x[j]) / twodt;
    }
    for (j = n - 1; j >= 2; --j) {
/* L2: */
	x[j] = x[j - 1];
    }
    x[n] = x[n - 1];
} /* seife_dif */


/* ----- END OF cseife_deriv.c ----- */
