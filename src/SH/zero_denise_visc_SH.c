/*-----------------------------------------------------------------------------
 *   Initialize wavefield and PML variables for the visco-elastic SH problem
 *  
 *  
 *   D. Koehn
 *   Kiel, 04.12.2017
 *
 *  --------------------------------------------------------------------------- */

#include "fd.h"

void zero_denise_visc_SH(int ny1, int ny2, int nx1, int nx2, float ** vz, float ** sxz, float ** syz, float ** vzm1, 
			 float ** vzp1, float ** psi_sxz_x, float ** psi_syz_y, float ** psi_vzx,  float ** psi_vzy, 
                         float ***pr, float ***pp, float ***pq, float ***Rxz, float ***Ryz){



	register int i, j, k, l;
	extern int FW, NX, NY, L, MODE;

	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
				vz[j][i]=0.0;
				sxz[j][i]=0.0;
				syz[j][i]=0.0;
                                vzm1[j][i]=0.0;
				vzp1[j][i]=0.0;
				
			}
		}
		
		for (j=1;j<=NY;j++){
		         for (i=1;i<=2*FW;i++){
		 
		                psi_sxz_x[j][i] = 0.0;
		                psi_vzx[j][i] = 0.0;   
		 
		         }
		}
		
		for (j=1;j<=2*FW;j++){
		         for (i=1;i<=NX;i++){
		                
		                
		                psi_syz_y[j][i] = 0.0;
		                psi_vzy[j][i] = 0.0;
		                
		         }
		}
		
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (l=1;l<=L;l++){
					pr[j][i][l] = 0.0;
					pp[j][i][l] = 0.0;
					pq[j][i][l] = 0.0;

					if(MODE==1 || MODE==2){
					  Rxz[j][i][l] = 0.0;
					  Ryz[j][i][l] = 0.0;
					}
				}
			}
		}		
					            
}
