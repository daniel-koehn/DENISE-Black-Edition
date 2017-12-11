
#include "fd.h"

void av_mu_SH(float ** u, float ** uip, float ** ujp, float ** rho){

	extern int NX, NY, INVMAT1;
	int i, j;
	float u1, u2, u3, u4;
	
   if(INVMAT1==3){
	
        for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
		       u1 = u[j][i];
		       u2 = u[j][i+1];
		       u3 = u[j+1][i];
		        
		       uip[j][i]=2.0/((1.0/u1)+(1.0/u2));
		       ujp[j][i]=2.0/((1.0/u1)+(1.0/u3));   
		      	
	
		}
	}
	
   }
	
   if(INVMAT1==1){
	
	for (j=1;j<=NY;j++){
	        for (i=1;i<=NX;i++){
	       
		       u1 = rho[j][i] * u[j][i] * u[j][i];
		       u2 = rho[j][i+1] * u[j][i+1] * u[j][i+1];
		       u3 = rho[j+1][i] * u[j+1][i] * u[j+1][i];
		        
		       uip[j][i]=2.0/((1.0/u1)+(1.0/u2));
		       ujp[j][i]=2.0/((1.0/u1)+(1.0/u3));
		       		
		       if((u1==0.0)||(u2==0.0)){
		           uip[j][i]=0.0;
		       }
		       
		       if((u1==0.0)||(u3==0.0)){
                           ujp[j][i]=0.0;		           
		       }
		      	
	
		}
	}
	
   }
	
}
