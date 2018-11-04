/*
 * Assemble inverse Pseudo-Hessians (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 27/10/2018
 */

#include "fd.h"

void apply_inv_hessSH(struct fwiSH *fwiSH, struct matSH *matSH, int nshots){

        /* global variables */
	extern int NX, NY, IDX, IDY, NTDTINV, MYID;
	extern float DT, C_vs, C_rho, C_taus, SCALERHO, SCALEQS;

	/* local variables */
	int i, j;
	float max_hess_vs2, max_hess_rho2, max_hess_ts2;
	float hess_max, hg_vs, hg_rho, hg_ts, vs2, rhovs;
	float **A, **Ainv;
	float EPS_HESS;
	
	A =  matrix(1,3,1,3);
	Ainv =  matrix(1,3,1,3);
	
	EPS_HESS = 1e-4;

	max_hess_vs2 = 0.0;
	max_hess_rho2 = 0.0;
	max_hess_ts2 = 0.0;
	
	/* Assemble Pseudo-Hessians for Vs-density-ts parametrization */
	/* ---------------------------------------------------------- */
	for (i=1;i<=NX;i=i+IDX){
   		for (j=1;j<=NY;j=j+IDY){
		 
		 	vs2 = (*matSH).pu[j][i] * (*matSH).pu[j][i];
		 
		        /* Calculate Pseudo-Hessian off-diagonal elements */
			/* ---------------------------------------------- */
			 (*fwiSH).hess_muts[j][i] = DT * (*fwiSH).hess_muts[j][i] / (nshots * NTDTINV);			      					       						
			(*fwiSH).hess_murho[j][i] = DT * (*fwiSH).hess_murho[j][i] / (nshots * NTDTINV);
			(*fwiSH).hess_tsrho[j][i] = DT * (*fwiSH).hess_tsrho[j][i] / (nshots * NTDTINV);
		 
      			/* Calculate Pseudo-Hessian vs-vs */ 
			/* ------------------------------ */
      			(*fwiSH).hess_mu2[j][i] = DT * (*fwiSH).hess_mu2[j][i] / (nshots * NTDTINV);
			(*fwiSH).hess_vs2[j][i] = 4.0 * (*matSH).prho[j][i] * (*matSH).prho[j][i] * vs2 * (*fwiSH).hess_mu2[j][i];
			
			/* Calculate Pseudo-Hessian rho-rho */ 
			/* -------------------------------- */
       			(*fwiSH).hess_rho2[j][i] = DT * (*fwiSH).hess_rho2[j][i] / (nshots * NTDTINV);
			(*fwiSH).hess_rho2p[j][i] = (vs2 * vs2 * (*fwiSH).hess_vs2[j][i]) + (2.0 * vs2 * (*fwiSH).hess_murho[j][i]) + (*fwiSH).hess_rho2[j][i];

			/* Calculate Pseudo-Hessian taus-taus */ 
			/* ---------------------------------- */
       			(*fwiSH).hess_ts2[j][i] = DT * (*fwiSH).hess_ts2[j][i] / (nshots * NTDTINV);						
			
			/* Estimate maximum Pseudo-Hessian values */
			/* -------------------------------------- */
			if(fabs((*fwiSH).hess_vs2[j][i])>max_hess_vs2){max_hess_vs2 = fabs((*fwiSH).hess_vs2[j][i]);}
			if(fabs((*fwiSH).hess_rho2p[j][i])>max_hess_rho2){max_hess_rho2 = fabs((*fwiSH).hess_rho2p[j][i]);}
			if(fabs((*fwiSH).hess_ts2[j][i])>max_hess_ts2){max_hess_ts2 = fabs((*fwiSH).hess_ts2[j][i]);}
		                                                                       
   		}
	}
	
	hess_max = 0.0;
        MPI_Allreduce(&max_hess_vs2,&hess_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	max_hess_vs2 = hess_max;
	
	hess_max = 0.0;
        MPI_Allreduce(&max_hess_rho2,&hess_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	max_hess_rho2 = hess_max;
	
	hess_max = 0.0;
        MPI_Allreduce(&max_hess_ts2,&hess_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	max_hess_ts2 = hess_max;
	
	/* Apply inverse Hessian to gradient */
	/* --------------------------------- */	
	for (i=1;i<=NX;i=i+IDX){
   		for (j=1;j<=NY;j=j+IDY){
		 
		/* Assemble Hessian (main diagonal) */		
		A[1][1] = (*fwiSH).hess_vs2[j][i] + EPS_HESS * max_hess_vs2;
		A[2][2] = SCALERHO * SCALERHO * ((*fwiSH).hess_rho2p[j][i] + EPS_HESS * max_hess_rho2);
		A[3][3] =  SCALEQS * SCALEQS  * ((*fwiSH).hess_ts2[j][i] + EPS_HESS * max_hess_ts2);				
		
		/* Assemble Hessian (off-diagonal blocks) */		
		vs2 = (*matSH).pu[j][i] * (*matSH).pu[j][i];
		rhovs = 2.0 * (*matSH).prho[j][i] * (*matSH).pu[j][i];
		
		A[1][2] = rhovs * ((*fwiSH).hess_murho[j][i]  + vs2 * (*fwiSH).hess_mu2[j][i]);
		A[1][3] = rhovs * (*fwiSH).hess_muts[j][i]; 
		A[2][3] = (*fwiSH).hess_tsrho[j][i] + vs2 * (*fwiSH).hess_muts[j][i];
		
	        A[1][2] = 0.0;
		A[1][3] = 0.0;
		A[2][3] = 0.0;
		
		A[2][1] = A[1][2];		
		A[3][1] = A[1][3];
		A[3][2] = A[2][3];		
		
		/* invert 3x3 Hessian matrix */
		mat_inv_3x3(A, Ainv);
		
		/*if((MYID==0)&&(i==1)&&(j==1)){		
		   
		   printf("Matrix A\n");
		   printf("%e \t %e \t %e \n",A[1][1],A[1][2],A[1][3]);
		   printf("%e \t %e \t %e \n",A[2][1],A[2][2],A[2][3]);		   
		   printf("%e \t %e \t %e \n\n",A[3][1],A[3][2],A[3][3]);
		   
		   printf("Matrix Ainv\n");
		   printf("%e \t %e \t %e \n",Ainv[1][1],Ainv[1][2],Ainv[1][3]);
		   printf("%e \t %e \t %e \n",Ainv[2][1],Ainv[2][2],Ainv[2][3]);		   
		   printf("%e \t %e \t %e \n",Ainv[3][1],Ainv[3][2],Ainv[3][3]);		
		   
		}*/
			
		 hg_vs = Ainv[1][1] * (*fwiSH).waveconv_u[j][i] + Ainv[1][2] * (*fwiSH).waveconv_rho[j][i] + Ainv[1][3] * (*fwiSH).waveconv_ts[j][i];
		hg_rho = Ainv[2][1] * (*fwiSH).waveconv_u[j][i] + Ainv[2][2] * (*fwiSH).waveconv_rho[j][i] + Ainv[2][3] * (*fwiSH).waveconv_ts[j][i];
                 hg_ts = Ainv[3][1] * (*fwiSH).waveconv_u[j][i] + Ainv[3][2] * (*fwiSH).waveconv_rho[j][i] + Ainv[3][3] * (*fwiSH).waveconv_ts[j][i];
		 
		  (*fwiSH).waveconv_u[j][i] = hg_vs;
		(*fwiSH).waveconv_rho[j][i] = hg_rho;
		 (*fwiSH).waveconv_ts[j][i] = hg_ts;		 
			
   		}
	}
	
	free_matrix(A,1,3,1,3);
	free_matrix(Ainv,1,3,1,3);
	
}
