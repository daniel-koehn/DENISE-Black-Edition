/*------------------------------------------------------------------------
 *   zero wavefield
 *  
 *  
 *   last update 06/10/11, L. Groos
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_fdveps_visc(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sxx, 
                 float ** syy, float ** sxy, float ** vxm1, float ** vym1, float ** vxp1, float ** vyp1,
                 float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy,
                 float ** psi_vxxs, float ***pr, float ***pp, float ***pq){



	register int i, j, k, l;
	extern int FW, NX, NY, L;

	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
				vx[j][i]=0.0;
				vy[j][i]=0.0;
				sxx[j][i]=0.0;
				syy[j][i]=0.0;
                                sxy[j][i]=0.0;
                                vxm1[j][i]=0.0;
				vym1[j][i]=0.0;
				vxp1[j][i]=0.0;
				vyp1[j][i]=0.0;
				
			}
		}
		
		for (j=1;j<=NY;j++){
		         for (i=1;i<=2*FW;i++){
		 
		                psi_sxx_x[j][i] = 0.0;
		                psi_sxy_x[j][i] = 0.0;
		                psi_vxx[j][i] = 0.0;
		                psi_vxxs[j][i] = 0.0;
		                psi_vyx[j][i] = 0.0;   
		 
		         }
		}
		
		for (j=1;j<=2*FW;j++){
		         for (i=1;i<=NX;i++){
		                
		                
		                psi_syy_y[j][i] = 0.0;
		                psi_sxy_y[j][i] = 0.0;
		                psi_vyy[j][i] = 0.0;
		                psi_vxy[j][i] = 0.0;
		                
		         }
		}
		
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (l=1;l<=L;l++){
					pr[j][i][l] = 0.0;
					pp[j][i][l] = 0.0;
					pq[j][i][l] = 0.0;
				}
			}
		}
					            
}
