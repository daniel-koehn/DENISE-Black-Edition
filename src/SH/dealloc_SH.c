/*
 * Deallocate memory for SH problem 
 *
 * Daniel Koehn
 * Kiel, 22/12/2017
 */

#include "fd.h"

void dealloc_SH(struct waveSH *waveSH, struct waveSH_PML *waveSH_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* deallocate memory */

        /* PML variables */
               if(FW>0){

		free_vector((*waveSH_PML).d_x,1,2*FW);
		free_vector((*waveSH_PML).K_x,1,2*FW);
		free_vector((*waveSH_PML).alpha_prime_x,1,2*FW);
		free_vector((*waveSH_PML).a_x,1,2*FW);
		free_vector((*waveSH_PML).b_x,1,2*FW);

		free_vector((*waveSH_PML).d_x_half,1,2*FW);
		free_vector((*waveSH_PML).K_x_half,1,2*FW);
		free_vector((*waveSH_PML).alpha_prime_x_half,1,2*FW);
		free_vector((*waveSH_PML).a_x_half,1,2*FW);
		free_vector((*waveSH_PML).b_x_half,1,2*FW);

		free_vector((*waveSH_PML).d_y,1,2*FW);
		free_vector((*waveSH_PML).K_y,1,2*FW);
		free_vector((*waveSH_PML).alpha_prime_y,1,2*FW);
		free_vector((*waveSH_PML).a_y,1,2*FW);
		free_vector((*waveSH_PML).b_y,1,2*FW);

		free_vector((*waveSH_PML).d_y_half,1,2*FW);
		free_vector((*waveSH_PML).K_y_half,1,2*FW);
		free_vector((*waveSH_PML).alpha_prime_y_half,1,2*FW);
		free_vector((*waveSH_PML).a_y_half,1,2*FW);
		free_vector((*waveSH_PML).b_y_half,1,2*FW);

		free_matrix((*waveSH_PML).psi_sxz_x,1,NY,1,2*FW);
		free_matrix((*waveSH_PML).psi_syz_y,1,2*FW,1,NX);
		free_matrix((*waveSH_PML).psi_vzy,1,2*FW,1,NX);
		free_matrix((*waveSH_PML).psi_vzx,1,NY,1,2*FW);

		// free_matrix((*waveSH_PML).absorb_coeff,1,NY,1,NX);

       }

	/* elastic wavefield variables */
	free_matrix((*waveSH).psxz,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).psyz,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).pvz,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).pvzp1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).pvzm1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).uz,-nd+1,NY+nd,-nd+1,NX+nd);
	// free_matrix((*waveSH).uxz,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).uzx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveSH).uttz,-nd+1,NY+nd,-nd+1,NX+nd);

        /*visco-elastic wavefield variables */
        if(L > 0) {
		free_f3tensor((*waveSH).pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor((*waveSH).pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor((*waveSH).pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	} 

}



