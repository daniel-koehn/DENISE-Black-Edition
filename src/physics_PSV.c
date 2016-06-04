/*
 * Forward, FWI, RTM and RTMOD modules (PSV problem)  
 *
 * Daniel Koehn
 * Kiel, 26/04/2016
 */

#include "fd.h"

void physics_PSV(){

	/* global variables */
	extern int INVMAT;

	/* 2D PSV Forward Problem */
	if(INVMAT==10){
	   FD_PSV();
	}

	/* 2D PSV Full Waveform Inversion */
	if(INVMAT==0){
	   FWI_PSV();
	}

        /* 2D PSV Reverse Time Migration */
	/* if(INVMAT==2){
	   FWI_PSV(fileinp1,req_send,req_rec,send_statuses,rec_statuses);
	}*/

}



