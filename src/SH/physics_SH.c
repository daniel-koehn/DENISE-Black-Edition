/*
 * Forward, FWI, RTM and FD-based FWI gradient modules (isotropic SH problem)  
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

extern int L;

void physics_SH(){

	/* global variables */
	extern int MODE;

	/* 2D SH Forward Problem */
	if(MODE==0){
	   FD_SH();
	}

	/* 2D SH Full Waveform Inversion */
	if(MODE==1){
	   if(L==0){FWI_SH();}   /* elastic SH-FWI */
	   if(L){FWI_SH_visc();} /* visco-elastic SH-FWI */
	}

        /* 2D SH Reverse Time Migration */
	/*if(MODE==2){
	   RTM_AC();
	}*/

        /* 2D SH FD-based FWI gradient computation */
	if(MODE==3){
	   FD_grad_SH();
	}


}



