/*------------------------------------------------------------------------
 *   Initialize wavefield and PML variables for the acoustic problem
 *  
 *  
 *   D. Koehn
 *   Kiel, 10.06.2017
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_denise_acoustic_AC(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** p, 
                 float ** vxm1, float ** vxp1, float ** vyp1,
                 float ** psi_p_x, float ** psi_vxx, float ** psi_p_y, float ** psi_vyy, float ** psi_vxxs){


	register int i, j, k;
	extern int FW, NX, NY;

	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
				vx[j][i]=0.0;
				vy[j][i]=0.0;
				p[j][i]=0.0;
                                vxm1[j][i]=0.0;
				vxp1[j][i]=0.0;
				vyp1[j][i]=0.0;
				
			}
		}
		
		for (j=1;j<=NY;j++){
		         for (i=1;i<=2*FW;i++){
		 
		                psi_p_x[j][i] = 0.0;
		                psi_vxx[j][i] = 0.0;
		                psi_vxxs[j][i] = 0.0;  
		 
		         }
		}
		
		for (j=1;j<=2*FW;j++){
		         for (i=1;i<=NX;i++){
		                
		                
		                psi_p_y[j][i] = 0.0;
		                psi_vyy[j][i] = 0.0;
		                
		         }
		}            
	
}
