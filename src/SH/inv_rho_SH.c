
#include "fd.h"

void inv_rho_SH(float ** rho, float ** rhoi){

	extern int NX, NY;
	int i, j;
	
        for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		   rhoi[j][i] = 1.0 / rho[j][i];	       	   

		   if(rho[j][i]<1e-4){

		       rhoi[j][i] = 0.0;

		   }      	
	
		}
	}	
		
}
