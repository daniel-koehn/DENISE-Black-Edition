/*
 * Deallocate memory for acoustic problem 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void dealloc_AC(struct waveAC *waveAC, struct waveAC_PML *waveAC_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* deallocate memory */

        /* PML variables */
        if(FW>0){

		free_vector((*waveAC_PML).d_x,1,2*FW);
		free_vector((*waveAC_PML).K_x,1,2*FW);
		free_vector((*waveAC_PML).alpha_prime_x,1,2*FW);
		free_vector((*waveAC_PML).a_x,1,2*FW);
		free_vector((*waveAC_PML).b_x,1,2*FW);

		free_vector((*waveAC_PML).d_x_half,1,2*FW);
		free_vector((*waveAC_PML).K_x_half,1,2*FW);
		free_vector((*waveAC_PML).alpha_prime_x_half,1,2*FW);
		free_vector((*waveAC_PML).a_x_half,1,2*FW);
		free_vector((*waveAC_PML).b_x_half,1,2*FW);

		free_vector((*waveAC_PML).d_y,1,2*FW);
		free_vector((*waveAC_PML).K_y,1,2*FW);
		free_vector((*waveAC_PML).alpha_prime_y,1,2*FW);
		free_vector((*waveAC_PML).a_y,1,2*FW);
		free_vector((*waveAC_PML).b_y,1,2*FW);

		free_vector((*waveAC_PML).d_y_half,1,2*FW);
		free_vector((*waveAC_PML).K_y_half,1,2*FW);
		free_vector((*waveAC_PML).alpha_prime_y_half,1,2*FW);
		free_vector((*waveAC_PML).a_y_half,1,2*FW);
		free_vector((*waveAC_PML).b_y_half,1,2*FW);

		free_matrix((*waveAC_PML).psi_p_x,1,NY,1,2*FW);
		free_matrix((*waveAC_PML).psi_p_y,1,2*FW,1,NX);
		free_matrix((*waveAC_PML).psi_vxx,1,NY,1,2*FW);
		free_matrix((*waveAC_PML).psi_vxxs,1,NY,1,2*FW);
		free_matrix((*waveAC_PML).psi_vyy,1,2*FW,1,NX);

       }

	/* elastic wavefield variables */
	free_matrix((*waveAC).p,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).ux,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).uttx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix((*waveAC).utty,-nd+1,NY+nd,-nd+1,NX+nd);

        /*visco-elastic wavefield variables */
        if(L > 0) {
		free_f3tensor((*waveAC).pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor((*waveAC).pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor((*waveAC).pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        }	

}



