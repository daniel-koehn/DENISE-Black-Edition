/*
 * Forward, FWI, RTM and RTMOD modules (VTI problem)  
 *
 * Daniel Koehn
 * Kiel, 01/02/2017
 */

#include "fd.h"

void physics_VTI(){

	/* global variables */
	extern int MODE;

	/* 2D VTI Forward Problem */
	if(MODE==0){
	   FD_VTI();
	}

	/* 2D PSV Full Waveform Inversion */
	/*if(MODE==1){
	   FWI_PSV();
	}*/

        /* 2D PSV Reverse Time Migration */
	if(MODE==2){
	   RTM_VTI();
	}

}



