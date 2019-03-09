/*------------------------------------------------------------------------
 * Smooth gradient
 * 
 * Daniel Koehn
 * Kiel, 12.08.2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void smooth_grad(float ** waveconv, float ** vel_mod){

	extern int SPATFILTER;
		
	/* Smooth gradients */
	/* ---------------- */

	/* apply wavenumber damping */
	if(SPATFILTER==1){
	  wavenumber(waveconv);
	}

	if(SPATFILTER==2){
	  smooth2(waveconv);
	}

	if(SPATFILTER==3){
	  gauss_filt(waveconv);
	}

	if(SPATFILTER==4){
	  gauss_filt_var(waveconv,vel_mod);
	}


}
