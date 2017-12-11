/*-----------------------------------------------------------------------------
 *   Initialize wavefield and PML variables for the elastic SH problem
 *  
 *  
 *   D. Koehn
 *   Kiel, 04.12.2017
 *
 *  --------------------------------------------------------------------------- */

#include "fd.h"

void zero_denise_elast_SH(int ny1, int ny2, int nx1, int nx2, float ** vz, float ** sxz, float ** syz, float ** vzm1, 
			 float ** vzp1, float ** psi_sxz_x, float ** psi_syz_y, float ** psi_vzx,  float ** psi_vzy){



	register int i, j, k, l;
	extern int FW, NX, NY, L;

	
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
					            
}
