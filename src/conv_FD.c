/*------------------------------------------------------------------------
 *   Function for convolving           
 *   19/01/02, T. Bohlen
 *   last update 19.09.11 M. Schäfer => FFTW3 check http://fftw.org/
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"
#include <fftw3.h>

void  conv_FD(float * temp_TS, float * temp_TS1, float * temp_conv, int ns){

	/* declaration of local variables */
	int i,j, h, nfreq, npad;
	float xr, yr, x, y, dump, a;
	double npadd;
	
	/* declaration of variables for FFTW3*/
	fftw_complex *in1, *in2, *out1, *out2;
	fftw_plan p1, p2, p3;
	
    
	npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
        npadd = (double)npad;
	
	in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	
	in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	
	
	for (j=0;j<ns;j++){
	in1[j][0]=temp_TS[j];
	in1[j][1]=0.0;
		}
	for (j=ns;j<npad;j++){
	in1[j][0]=0.0;
	in1[j][1]=0.0;
		}
	
	for (j=0;j<ns;j++){
	in2[j][0]=temp_TS1[j];
	in2[j][1]=0.0;
		}
	for (j=ns;j<npad;j++){
	in2[j][0]=0.0;
	in2[j][1]=0.0;
		}
	
	
	p1 = fftw_plan_dft_1d(npad, in1, out1, 1, FFTW_ESTIMATE);
	fftw_execute(p1);
	
	p2 = fftw_plan_dft_1d(npad, in2, out2, 1, FFTW_ESTIMATE);
	fftw_execute(p2);
		
	/* multiplication of complex vectors */
	for(i=0;i<npad;i++){
	
	out1[i][0]    = (out1[i][0] * out2[i][0]) - (out1[i][1] * out2[i][1]);
	out1[i][1] = (out1[i][0] * out2[i][0]) + (out1[i][1] * out2[i][1]);
		   
			}
	
	p3 = fftw_plan_dft_1d(npad, out1, in1, -1, FFTW_ESTIMATE);
	fftw_execute(p3);
		
	 /* write output into data matrix */
        for(j=0;j<ns;j++){
	       temp_conv[j] = in1[j][0];
	    }
	
	
	fftw_free(in1);
	fftw_free(out1);
	
	fftw_free(in2);
	fftw_free(out2);
	
	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	
}
