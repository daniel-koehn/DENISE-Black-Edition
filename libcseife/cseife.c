/*! \file cseife.c
 * \brief make seife functions available within C code (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cseife.c 3381 2010-11-15 12:33:00Z tforb $
 * \author Thomas Forbriger
 * \date 14/01/2005
 * 
 * make seife functions available within C code (implementation)
 * 
 * Copyright (c) 2005 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 14/01/2005   V1.0   Thomas Forbriger
 *  - 11/07/2005   V1.1   support debug mode
 *  - 15/11/2010   V1.2   do not use tfmacros.h
 * 
 * ============================================================================
 */
#define TF_CSEIFE_C_VERSION \
  "TF_CSEIFE_C   V1.2"
#define TF_CSEIFE_C_CVSID \
  "$Id: cseife.c 3381 2010-11-15 12:33:00Z tforb $"

#include <cseife.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* initialize global flags */
struct seife_flags seife_global_flags={ 0 };

/* enter debug mode */
void seife_debug_mode_on()
{
  seife_global_flags.debug=1;
}

/* Butterworth lowpass (period t0, order o) */
void seife_lpb(double* x, int n, double dt, double t0, int o)
{
  seife_butterworth(x, n, dt, IIRlpb, t0, o);
}

/* Butterworth highpass (period t0, order o) */
void seife_hpb(double* x, int n, double dt, double t0, int o)
{
  seife_butterworth(x, n, dt, IIRhpb, t0, o);
}

/* 2nd order lowpass (period t0, damping h) */
void seife_lp2(double* x, int n, double dt, double t0, double h)
{
  seife_rekfl(x, x, n, seife_rfk(IIRlp2, t0/dt, h, 0., 0.));
}

/* 2nd order highpass (period t0, damping h) */
void seife_hp2(double* x, int n, double dt, double t0, double h)
{
  seife_first(x, n);
  seife_rekfl(x, x, n, seife_rfk(IIRhp2, t0/dt, h, 0., 0.));
}

/* 2nd order bandpass (period t0, damping h) */
void seife_bp2(double* x, int n, double dt, double t0, double h)
{
  seife_rekfl(x, x, n, seife_rfk(IIRbp2, t0/dt, h, 0., 0.));
}

/* 1st order lowpass (period t0) */
void seife_lp1(double* x, int n, double dt, double t0)
{
  seife_rekfl(x, x, n, seife_rfk(IIRlp1, t0/dt, 0., 0., 0.));
}

/* 1st order highpass (period t0) */
void seife_hp1(double* x, int n, double dt, double t0)
{
  seife_first(x, n);
  seife_rekfl(x, x, n, seife_rfk(IIRhp1, t0/dt, 0., 0., 0.));
}

/* integration (time constant t0) */
void seife_int(double* x, int n, double dt, double t0)
{
  double t00=t0;
  if (t00 < 1.e-22) { t00=1.; }
  seife_rekfl(x, x, n, seife_rfk(IIRint, t00/dt, 0., 0., 0.));
}

/* 1st order highpass equalizer (former period t0s, new period t0) */
void seife_he1(double* x, int n, double dt, double t0s, double t0)
{
  seife_rekfl(x, x, n, seife_rfk(IIRhe1, t0/dt, 0., t0s/dt, 0.));
}

/* 1st order lowpass equalizer (former period t0s, new period t0) */
void seife_le1(double* x, int n, double dt, double t0s, double t0)
{
  struct seife_IIRcoef coef=seife_rfk(IIRle1, t0/dt, 0., t0s/dt, 0.);
  double fac=t0s/t0;
  coef.f0 *= fac;
  coef.f1 *= fac;
  seife_rekfl(x, x, n, coef);
}

/* 2nd order highpass equalizer (former: period t0s and damping hs,
 *                               new: period t0 and damping h) 
 */
void seife_he2(double* x, int n, double dt, 
               double t0s, double hs, double t0, double h)
{
  if (seife_global_flags.debug)
  {
    fprintf(stderr, 
            "DEBUG (seife_he2): n=%d, dt=%f, t0s=%f, hs=%f, t0=%f, h=%f\n",
            n, dt, t0s, hs, t0, h);
  }
  seife_rekfl(x, x, n, seife_rfk(IIRhe2, t0/dt, h, t0s/dt, hs));
}

/* 2nd order lowpass equalizer (former: period t0s and damping hs,
 *                              new: period t0 and damping h) 
 */
void seife_le2(double* x, int n, double dt, 
               double t0s, double hs, double t0, double h)
{
  struct seife_IIRcoef coef=seife_rfk(IIRle2, t0/dt, h, t0s/dt, hs);
  double fac=t0s/t0;
  fac *= fac;
  coef.f0 *= fac;
  coef.f1 *= fac;
  coef.f2 *= fac;
  seife_rekfl(x, x, n, coef);
}

/* set baseline to first value */
void seife_first(double* x, int n)
{
  double x0=x[0];
  int i;
  for (i=0; i<n; ++i) { x[i] -= x0; }
}

/* general Butterworth filter framework */
void seife_butterworth(double* x, int n, double dt,
                       enum seife_filter_type IIRtype, 
                       double t0, int m)
{
  enum seife_filter_type t1st, t2nd;
  if (IIRtype == IIRhpb)
  {
    seife_first(x, n);
    t1st=IIRhp1;
    t2nd=IIRhp2;
  }
  else if (IIRtype == IIRlpb)
  {
    t1st=IIRlp1;
    t2nd=IIRlp2;
  }
  else
  {
    SEIFE_CHECKERROR( 1, "seife_butterworth", "illegal IIRtype" )
  }
  int mm=m/2;
  if (m > (2*mm))
  { seife_rekfl(x, x, n, seife_rfk(t1st, t0/dt, 0., 0., 0.)); }
  double pih=2.*atan(1.);
  int j;
  for (j=1; j<=mm; ++j)
  {
    double h=sin(pih*(2*j-1)/m);
    seife_rekfl(x, x, n, seife_rfk(t2nd, t0/dt, h, 0., 0.)); 
  }
}


/* ----- END OF cseife.c ----- */
