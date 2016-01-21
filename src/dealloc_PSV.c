/*
 * Deallocate memory for PSV problem 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void dealloc_PSV(struct wavePSV_el wavePSV_el, struct wavePSV_visc wavePSV_visc, struct wavePSV_PML wavePSV_PML){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

        /* deallocate memory */

        /* PML variables */
        if(FW>0){

		free_vector(wavePSV_PML.d_x,1,2*FW);
		free_vector(wavePSV_PML.K_x,1,2*FW);
		free_vector(wavePSV_PML.alpha_prime_x,1,2*FW);
		free_vector(wavePSV_PML.a_x,1,2*FW);
		free_vector(wavePSV_PML.b_x,1,2*FW);

		free_vector(wavePSV_PML.d_x_half,1,2*FW);
		free_vector(wavePSV_PML.K_x_half,1,2*FW);
		free_vector(wavePSV_PML.alpha_prime_x_half,1,2*FW);
		free_vector(wavePSV_PML.a_x_half,1,2*FW);
		free_vector(wavePSV_PML.b_x_half,1,2*FW);

		free_vector(wavePSV_PML.d_y,1,2*FW);
		free_vector(wavePSV_PML.K_y,1,2*FW);
		free_vector(wavePSV_PML.alpha_prime_y,1,2*FW);
		free_vector(wavePSV_PML.a_y,1,2*FW);
		free_vector(wavePSV_PML.b_y,1,2*FW);

		free_vector(wavePSV_PML.d_y_half,1,2*FW);
		free_vector(wavePSV_PML.K_y_half,1,2*FW);
		free_vector(wavePSV_PML.alpha_prime_y_half,1,2*FW);
		free_vector(wavePSV_PML.a_y_half,1,2*FW);
		free_vector(wavePSV_PML.b_y_half,1,2*FW);

		free_matrix(wavePSV_PML.psi_sxx_x,1,NY,1,2*FW);
		free_matrix(wavePSV_PML.psi_syy_y,1,2*FW,1,NX);
		free_matrix(wavePSV_PML.psi_sxy_x,1,NY,1,2*FW);
		free_matrix(wavePSV_PML.psi_sxy_y,1,2*FW,1,NX);
		free_matrix(wavePSV_PML.psi_vxx,1,NY,1,2*FW);
		free_matrix(wavePSV_PML.psi_vxxs,1,NY,1,2*FW);
		free_matrix(wavePSV_PML.psi_vyy,1,2*FW,1,NX);
		free_matrix(wavePSV_PML.psi_vxy,1,2*FW,1,NX);
		free_matrix(wavePSV_PML.psi_vyx,1,NY,1,2*FW);
		/*free_matrix(wavePSV_PML.absorb_coeff,1,NY,1,NX);*/

       }

	/* elastic wavefield variables */
	free_matrix(wavePSV_el.psxx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.psxy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.psyy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.ux,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.uy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.uxy,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.uyx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.uttx,-nd+1,NY+nd,-nd+1,NX+nd);
	free_matrix(wavePSV_el.utty,-nd+1,NY+nd,-nd+1,NX+nd);

        /*visco-elastic wavefield variables */
        if(L > 0) {
		free_f3tensor(wavePSV_visc.pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(wavePSV_visc.pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(wavePSV_visc.pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        }	

}



