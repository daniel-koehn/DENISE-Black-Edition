#include "fd.h"
void eprecond_SH(float ** W, float ** vz){
  
	extern int NX, NY, IDX, IDY;
	int i, j, k, l;
	
		
	        for (i=1;i<=NX;i=i+IDX){
		        for (j=1;j<=NY;j=j+IDY){
					     
			       W[j][i]+=(vz[j][i]*vz[j][i]);
		
			    }
		    }				
}
