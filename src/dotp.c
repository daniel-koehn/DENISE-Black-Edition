/* Calculate dot product of a part of two vectors
 *
 * sw = 0 -> vectors have the same length
 * sw = 1 -> vec1 has length h1:h2, vec2 has length 1:n
 *
 * Daniel Koehn
 * Kiel, the 22nd of July 2013
 */

#include "fd.h"

float dotp(float * vec1, float *vec2, int n1, int n2, int sw)
{

	/* extern variables */
	/*extern int NPROCX, NPROCY, MYID, POS[3];*/
	
	/* local variables */
	int i,h;
	float tmp,sum;
	
	tmp=0.0;
	
	if(sw==0){
	  /* estimate dot product for vectors with same length*/
	  for (i=n1;i<=n2;i++){
	       tmp += vec1[i]*vec2[i];       
          }
	}
	
	if(sw==1){  
	  h=1;
          /* estimate dot product for vector of diffent length*/
          for (i=n1;i<=n2;i++){
               tmp += vec1[i]*vec2[h];
               h++;
          }
        }   
	
   /* Compute maximum value on all CPUs */
   sum = 0.0;
   MPI_Allreduce(&tmp,&sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
   
   return sum;
}



