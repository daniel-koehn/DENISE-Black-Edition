/* Extend density model */
/*------------------------------------------------------------------------
 * extend density model for gravity modelling to reduce boundary effects 
 *
 * Daniel Koehn
 * Kiel, the 8th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void extend_mod(float  **rho_grav, float  **rho_grav_ext, int nxgrav, int nygrav){

	int i, j;
	extern int NXG, NYG, NGRAVB;
        
	/* copy original model into larger model */
        for(i=1;i<=NXG;i++){       
		for(j=1;j<=NYG;j++){            
                
                    rho_grav_ext[j][i+NGRAVB] = rho_grav[j][i];     
		
        	}
	}
        
        
        /* add left boundary */
	 for(i=1;i<=NGRAVB;i++){       
		for(j=1;j<=nygrav;j++){            
                 
                   rho_grav_ext[j][i] = rho_grav_ext[j][NGRAVB+1];
                        
        	}
	}

        /* add right boundary */
	 for(i=NXG+NGRAVB+1;i<=nxgrav;i++){       
		for(j=1;j<=nygrav;j++){            
                 
                   rho_grav_ext[j][i] = rho_grav_ext[j][NXG+NGRAVB];
                        
        	}
	}
		
        /* add bottom boundary */
	 for(i=1;i<=nxgrav;i++){       
		for(j=NYG+1;j<=nygrav;j++){            
                 
                   rho_grav_ext[j][i] = rho_grav_ext[NYG][i];
                        
        	}
	}

}
