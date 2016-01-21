/*
 * Allocate memory for PSV problem 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void alloc_PSV(struct wavePSV_el wavePSV_el, struct wavePSV_visc wavePSV_visc, struct wavePSV_PML wavePSV_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* memory allocation for elastic wavefield variables */
        wavePSV_el.psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.uy   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.uxy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.uyx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	wavePSV_el.utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

        /* memory allocation for visco-elastic wavefield variables */
	if (L > 0) {
	    wavePSV_visc.pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    wavePSV_visc.pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    wavePSV_visc.pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	}

        /* memory allocation for PML variables */
        if(FW>0){

	  wavePSV_PML.d_x = vector(1,2*FW);
	  wavePSV_PML.K_x = vector(1,2*FW);
	  wavePSV_PML.alpha_prime_x = vector(1,2*FW);
	  wavePSV_PML.a_x = vector(1,2*FW);
	  wavePSV_PML.b_x = vector(1,2*FW);
	  
	  wavePSV_PML.d_x_half = vector(1,2*FW);
	  wavePSV_PML.K_x_half = vector(1,2*FW);
	  wavePSV_PML.alpha_prime_x_half = vector(1,2*FW);
	  wavePSV_PML.a_x_half = vector(1,2*FW);
	  wavePSV_PML.b_x_half = vector(1,2*FW);

	  wavePSV_PML.d_y = vector(1,2*FW);
	  wavePSV_PML.K_y = vector(1,2*FW);
	  wavePSV_PML.alpha_prime_y = vector(1,2*FW);
	  wavePSV_PML.a_y = vector(1,2*FW);
	  wavePSV_PML.b_y = vector(1,2*FW);
	  
	  wavePSV_PML.d_y_half = vector(1,2*FW);
	  wavePSV_PML.K_y_half = vector(1,2*FW);
	  wavePSV_PML.alpha_prime_y_half = vector(1,2*FW);
	  wavePSV_PML.a_y_half = vector(1,2*FW);
	  wavePSV_PML.b_y_half = vector(1,2*FW);

	  wavePSV_PML.psi_sxx_x =  matrix(1,NY,1,2*FW); 
	  wavePSV_PML.psi_syy_y =  matrix(1,2*FW,1,NX);
	  wavePSV_PML.psi_sxy_y =  matrix(1,2*FW,1,NX);
	  wavePSV_PML.psi_sxy_x =  matrix(1,NY,1,2*FW);
	  wavePSV_PML.psi_vxx   =  matrix(1,NY,1,2*FW);
	  wavePSV_PML.psi_vxxs  =  matrix(1,NY,1,2*FW); 
	  wavePSV_PML.psi_vyy   =  matrix(1,2*FW,1,NX);
	  wavePSV_PML.psi_vxy   =  matrix(1,2*FW,1,NX);
	  wavePSV_PML.psi_vyx   =  matrix(1,NY,1,2*FW);
   
        }
	
}



