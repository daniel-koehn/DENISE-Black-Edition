/*
 * Forward, FWI, RTM and RTMOD modules (PSV TTI problem)  
 *
 * Daniel Koehn
 * Kiel, 04/02/2017
 */

#include "fd.h"

void physics_TTI(){

	/* global variables */
	extern int MODE;

	/* 2D PSV TTI Forward Problem */
	if(MODE==0){
	   FD_TTI();
	}

	/* 2D PSV Full Waveform Inversion */
	/*if(MODE==1){
	   FWI_PSV();
	}*/

        /* 2D PSV TTI Reverse Time Migration */
	if(MODE==2){
	   RTM_TTI();
	}

}



