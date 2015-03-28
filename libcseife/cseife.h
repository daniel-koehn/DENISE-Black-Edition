/*! \file cseife.h
 * \brief make seife functions available within C code (prototypes)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cseife.h 3381 2010-11-15 12:33:00Z tforb $
 * \author Thomas Forbriger
 * \date 14/01/2005
 * 
 * make seife functions available within C code (prototypes)
 * 
 * Copyright (c) 2005 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 14/01/2005   V1.0   Thomas Forbriger
 *  - 11/07/2005   V1.1   introduce global flags for debug mode
 *  - 15/11/2010   V1.2   provide CHECKERROR macro
 * 
 * ============================================================================
 */

/* include guard */
#ifndef TF_CSEIFE_H_VERSION

#define TF_CSEIFE_H_VERSION \
  "TF_CSEIFE_H   V1.2"
#define TF_CSEIFE_H_CVSID \
  "$Id: cseife.h 3381 2010-11-15 12:33:00Z tforb $"

#define SEIFE_EXIT_FAILURE 1

#define SEIFE_CHECKERROR( EXPR , SUB, STR )\
    if ( EXPR ) { fprintf(stderr, "ERROR (%s):\n   %s\n", SUB, STR );\
        exit(SEIFE_EXIT_FAILURE); }


/* seife filter types */
enum seife_filter_type {
  IIRlpb, IIRhpb, IIRlp2, IIRhp2, IIRint, IIRhe1, IIRle1, IIRhe2, IIRle2,
  IIRlp1, IIRhp1, IIRbp2
}; /* enum seife_filter_type */

/* global flags (to be initialized in cseife.c) */
extern struct seife_flags { int debug; } seife_global_flags;

/* structure to hold IIR filter coefficients */
struct seife_IIRcoef { double f0, f1, f2, g1, g2; };

/* enter debug mode */
void seife_debug_mode_on();
/* Butterworth lowpass (period t0, order o) */
void seife_lpb(double* x, int n, double dt, double t0, int o);
/* Butterworth highpass (period t0, order o) */
void seife_hpb(double* x, int n, double dt, double t0, int o);
/* 2nd order lowpass (period t0, damping h) */
void seife_lp2(double* x, int n, double dt, double t0, double h);
/* 2nd order highpass (period t0, damping h) */
void seife_hp2(double* x, int n, double dt, double t0, double h);
/* 2nd order bandpass (period t0, damping h) */
void seife_bp2(double* x, int n, double dt, double t0, double h);
/* 1st order lowpass (period t0) */
void seife_lp1(double* x, int n, double dt, double t0);
/* 1st order highpass (period t0) */
void seife_hp1(double* x, int n, double dt, double t0);
/* integration (time constant t0) */
void seife_int(double* x, int n, double dt, double t0);
/* 1st order highpass equalizer (former period t0s, new period t0) */
void seife_he1(double* x, int n, double dt, double t0s, double t0);
/* 1st order lowpass equalizer (former period t0s, new period t0) */
void seife_le1(double* x, int n, double dt, double t0s, double t0);
/* 2nd order highpass equalizer (former: period t0s and damping hs,
 *                               new: period t0 and damping h) 
 */
void seife_he2(double* x, int n, double dt, 
               double t0s, double hs, double t0, double h);
/* 2nd order lowpass equalizer (former: period t0s and damping hs,
 *                              new: period t0 and damping h) 
 */
void seife_le2(double* x, int n, double dt, 
               double t0s, double hs, double t0, double h);
/* detide with synthetic tides interpolated over ni samples */
void seife_tid(double* x, int n, double dt, int ni);
/* derivative (time constant t0) */
void seife_dif(double* x, int n, double dt, double t0);
/* set baseline to first value */
void seife_first(double* x, int n);
/* solve system of linear equations */
void seife_gauss(double *aik, int m, int n, double* rs, double* f);
/* apply IIR filter coefficients */
void seife_rekfl(double* x, double* y, int n, struct seife_IIRcoef c);
/* determine IIR filter coefficients */
struct seife_IIRcoef seife_rfk(enum seife_filter_type IIRtype, 
                               double t0, double h, double t0s, double hs);
/* general Butterworth filter framework */
void seife_butterworth(double* x, int n, double dt,
                       enum seife_filter_type IIRtype, 
                       double t0, int o);

#endif /* TF_CSEIFE_H_VERSION (includeguard)  */

/* ----- END OF cseife.h ----- */
