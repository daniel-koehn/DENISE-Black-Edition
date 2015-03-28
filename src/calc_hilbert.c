/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
 *
 * This file is part of DENISE.
 * 
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Calculate the Hilbert transform of a signal                                  
 *   last update 27/1/10, L. Rehor
 *   update to FFTW3 01/07/2013 Martin Schaefer
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include <fftw3.h>

void calc_hilbert(float ** datatrace, float ** envelope, int ns, int ntr){

	/* declaration of variables */
	int i,j, nfreq, npad, k;
	float xr, yr, x, y, dump, a, *h;
	double npadd;
		
	/* declaration of variables for FFTW3*/
	fftw_complex *in, *out;
	fftw_plan p;

	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	/* calculation of the Hilbert transform */
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

	/* calculation of the FFT */
	npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
        npadd = (double)npad;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	
	/* calculation of the vector h */
	     h = vector(1,npad);
	     if ((2*(npad/2))==npad) {
		/* npad is even */
		h[1]=1.0;
		h[npad/2+1]=1.0;
		for (i=2;i<=(npad/2);i++) {
	   	     h[i] = 2.0; }
	     } else {
		/* npad is odd */
		h[1]=1.0;
		for (i=2;i<=((npad+1)/2);i++) {
	   	     h[i]=2.0;}
	     }
		
	for(k=1;k<=ntr;k++){
						
	     for (j=0;j<ns;j++){
	     in[j][0]=datatrace[k][j+1];
	     in[j][1]=0.0;
		     }
	     for (j=ns;j<npad;j++){
	     in[j][0]=0.0;
	     in[j][1]=0.0;
		     }
	     /* FFT */
	     p = fftw_plan_dft_1d(npad, in, out, 1, FFTW_ESTIMATE);
	     fftw_execute(p);
	     fftw_destroy_plan(p);
	     
	     /* elementwise multiplication of FFT and h */
	     for (i=0;i<npad;i++){
		out[i][0] *= h[i+1];
		out[i][1] *= h[i+1];
	     }
	     
	     p = fftw_plan_dft_1d(npad, out, in, -1, FFTW_ESTIMATE);
	     fftw_execute(p);
	     fftw_destroy_plan(p);		    

	     /* %%%%%%%%%%%%%%%%%%%%%% */
       	     /* taking complex part    */
             /* %%%%%%%%%%%%%%%%%%%%%% */
             for (j=0;j<=ns;j++){
             envelope[k][j+1]=(1.0/npadd)*in[j][1]*(-1);
             }

	}
	fftw_free(in);
	fftw_free(out);
	free_vector(h,1,npad);
}
