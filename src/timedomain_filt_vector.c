/*------------------------------------------------------------------------
 *   Filter seismograms in time domain with a Butterworth filter
 *   Lowpass or highpass filtering can be applied                                  
 *   last update 14/06/11, L. Rehor
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"
#include "cseife.h"

void  timedomain_filt_vector(float * data, float fc, int order, int ntr, int ns, int method){

	/* 
	data	: 	2-dimensional array containing seismograms (
	fc	:	corner frequency in Hz
	order	:	order of filter
	ntr	:	number of traces
	ns	:	number of samples
	method	:	definition of filter
			1: lowpass filter
			2: highpass filter
	*/
			

	/* declaration of extern variables */
	extern float DT;
	
	/* declaration of local variables */
	int itr, j;
	double *seismogram, T0;
	
	
	
	seismogram=dvector(1,ns);
	
	T0=1.0/fc;
	
	if ((order%2)!=0){
		err("Order of timedomain filter must be an even number!");}
		
			
	
		for (j=1;j<=ns;j++){
			seismogram[j]=(double)data[j];}
		
		if (method==1){		/*lowpass filter*/
			seife_lpb(seismogram,ns,DT,T0,order);}
		
		
		if (method==2){		/*highpass filter*/
			seife_hpb(seismogram,ns,DT,T0,order);}
		
		for (j=1;j<=ns;j++){
			data[j]=(float)seismogram[j];}
	
	
	free_dvector(seismogram,1,ns);
}
		
