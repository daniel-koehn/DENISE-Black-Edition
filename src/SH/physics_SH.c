/*
 * Forward, FWI, RTM and RTMOD modules (isotropic SH problem)  
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

void physics_SH(){

	/* global variables */
	extern int MODE;

	/* 2D SH Forward Problem */
	if(MODE==0){
	   FD_SH();
	}

	/* 2D SH Full Waveform Inversion */
	if(MODE==1){
	   FWI_SH();
	}

        /* 2D SH Reverse Time Migration */
	/*if(MODE==2){
	   RTM_AC();
	}*/

}



