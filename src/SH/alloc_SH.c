/*
 * Allocate memory for SH problem 
 *
 * Daniel Koehn
 * Kiel, 03/12/2017
 */

#include "fd.h"

void alloc_SH(struct waveSH *waveSH, struct waveSH_PML *waveSH_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* memory allocation for elastic wavefield variables */
        (*waveSH).psxz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).psyz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).pvz  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).pvzp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).pvzm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).uz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).uzx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*waveSH).uttz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

        /* memory allocation for visco-elastic wavefield variables */
	if (L > 0) {
	    (*waveSH).pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    (*waveSH).pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	    (*waveSH).pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	}

        /* memory allocation for PML variables */
        if(FW>0){

	  (*waveSH_PML).d_x = vector(1,2*FW);
	  (*waveSH_PML).K_x = vector(1,2*FW);
	  (*waveSH_PML).alpha_prime_x = vector(1,2*FW);
	  (*waveSH_PML).a_x = vector(1,2*FW);
	  (*waveSH_PML).b_x = vector(1,2*FW);
	  
	  (*waveSH_PML).d_x_half = vector(1,2*FW);
	  (*waveSH_PML).K_x_half = vector(1,2*FW);
	  (*waveSH_PML).alpha_prime_x_half = vector(1,2*FW);
	  (*waveSH_PML).a_x_half = vector(1,2*FW);
	  (*waveSH_PML).b_x_half = vector(1,2*FW);

	  (*waveSH_PML).d_y = vector(1,2*FW);
	  (*waveSH_PML).K_y = vector(1,2*FW);
	  (*waveSH_PML).alpha_prime_y = vector(1,2*FW);
	  (*waveSH_PML).a_y = vector(1,2*FW);
	  (*waveSH_PML).b_y = vector(1,2*FW);
	  
	  (*waveSH_PML).d_y_half = vector(1,2*FW);
	  (*waveSH_PML).K_y_half = vector(1,2*FW);
	  (*waveSH_PML).alpha_prime_y_half = vector(1,2*FW);
	  (*waveSH_PML).a_y_half = vector(1,2*FW);
	  (*waveSH_PML).b_y_half = vector(1,2*FW);

	  (*waveSH_PML).psi_syz_y =  matrix(1,2*FW,1,NX);
	  (*waveSH_PML).psi_sxz_x =  matrix(1,NY,1,2*FW);
	  (*waveSH_PML).psi_vzy   =  matrix(1,2*FW,1,NX);
	  (*waveSH_PML).psi_vzx   =  matrix(1,NY,1,2*FW);
   
        }
	
}



