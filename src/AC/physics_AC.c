/*
 * Forward, FWI, RTM and RTMOD modules (acoustic problem)  
 *
 * Daniel Koehn
 * Kiel, 10/06/2017
 */

#include "fd.h"

void physics_AC(){

	/* global variables */
	extern int MODE;

	/* 2D AC Forward Problem */
	if(MODE==0){
	   FD_AC();
	}

	/* 2D AC Full Waveform Inversion */
	if(MODE==1){
	   FWI_AC();
	}

        /* 2D AC Reverse Time Migration */
	if(MODE==2){
	   RTM_AC();
	}

}



