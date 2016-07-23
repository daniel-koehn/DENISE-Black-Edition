/*
 * Forward, FWI, RTM and RTMOD modules (PSV problem)  
 *
 * Daniel Koehn
 * Kiel, 26/04/2016
 */

#include "fd.h"

void physics_PSV(){

	/* global variables */
	extern int MODE;

	/* 2D PSV Forward Problem */
	if(MODE==0){
	   FD_PSV();
	}

	/* 2D PSV Full Waveform Inversion */
	if(MODE==1){
	   FWI_PSV();
	}

        /* 2D PSV Reverse Time Migration */
	if(MODE==2){
	   RTM_PSV();
	}

}



