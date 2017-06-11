/*
 * Allocate memory for AC problem 
 *
 * Daniel Koehn
 * Kiel, 10/06/2017
 */

#include "fd.h"

void alloc_AC(struct waveAC *waveAC, struct waveAC_PML *waveAC_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* memory allocation for elastic wavefield variables */
        (*waveAC).p =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveAC).utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

        /* memory allocation for visco-elastic wavefield variables */
	if (L > 0) {
	    (*waveAC).pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    (*waveAC).pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    (*waveAC).pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	}

        /* memory allocation for PML variables */
        if(FW>0){

	  (*waveAC_PML).d_x = vector(1,2*FW);
	  (*waveAC_PML).K_x = vector(1,2*FW);
	  (*waveAC_PML).alpha_prime_x = vector(1,2*FW);
	  (*waveAC_PML).a_x = vector(1,2*FW);
	  (*waveAC_PML).b_x = vector(1,2*FW);
	  
	  (*waveAC_PML).d_x_half = vector(1,2*FW);
	  (*waveAC_PML).K_x_half = vector(1,2*FW);
	  (*waveAC_PML).alpha_prime_x_half = vector(1,2*FW);
	  (*waveAC_PML).a_x_half = vector(1,2*FW);
	  (*waveAC_PML).b_x_half = vector(1,2*FW);

	  (*waveAC_PML).d_y = vector(1,2*FW);
	  (*waveAC_PML).K_y = vector(1,2*FW);
	  (*waveAC_PML).alpha_prime_y = vector(1,2*FW);
	  (*waveAC_PML).a_y = vector(1,2*FW);
	  (*waveAC_PML).b_y = vector(1,2*FW);
	  
	  (*waveAC_PML).d_y_half = vector(1,2*FW);
	  (*waveAC_PML).K_y_half = vector(1,2*FW);
	  (*waveAC_PML).alpha_prime_y_half = vector(1,2*FW);
	  (*waveAC_PML).a_y_half = vector(1,2*FW);
	  (*waveAC_PML).b_y_half = vector(1,2*FW);

	  (*waveAC_PML).psi_p_x =  matrix(1,NY,1,2*FW); 
	  (*waveAC_PML).psi_p_y =  matrix(1,2*FW,1,NX);
	  (*waveAC_PML).psi_vxx   =  matrix(1,NY,1,2*FW);
	  (*waveAC_PML).psi_vxxs  =  matrix(1,NY,1,2*FW); 
	  (*waveAC_PML).psi_vyy   =  matrix(1,2*FW,1,NX);
   
        }
	
}



