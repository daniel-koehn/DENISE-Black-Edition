/*
 * Normalize Gradient
 *
 * Daniel Koehn
 */

#include "fd.h"

void norm(float ** waveconv)
{

   /* extern variables */
   extern int NX, NY;
   extern int MYID;
	
   /* local variables */
   int i, j;
   float max, norm, tmp_max;
	
   max = fabs(waveconv[1][1]);

   /* estimate maximum of the global gradient */
   for (i=1;i<=NX;i++){
      for (j=1;j<=NY;j++){
	            
	 if(fabs(waveconv[j][i])>max){max=fabs(waveconv[j][i]);}
	            
      }
   }
	
   /* Compute maximum value on all CPUs */
   tmp_max = 0.0;
   MPI_Allreduce(&max,&tmp_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
   norm = tmp_max;
	
   for (i=1;i<=NX;i++){
      for (j=1;j<=NY;j++){	    	    
	  waveconv[j][i] = waveconv[j][i]/norm;	    
      }
   }                
	
   /* printf("MYID=%d \t norm_fac=%e \n",MYID,norm); */
}
