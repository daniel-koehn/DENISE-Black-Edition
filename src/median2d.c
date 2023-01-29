/*****************************************************************************/
/*         calculate the median of a 2D matrix                               */
/* Median code from https://codingpointer.com/c-tutorial/median              */
/*****************************************************************************/


#include "fd.h"

float median2d(float **mat, int ny, int nx){

	int i, j, k, n;
	float *a, temp, med;

	n = nx*ny;
	a = vector(0,n-1);
	
	k = 0;
	for (j=1;j<=ny;j++){
		for (i=1;i<=nx;i++){
		
			a[k] = mat[j][i];
			k++;
		}
	}			
		
	temp = 0;                                                                        
	/* sorting */                                                                   
	for(i=0;i<n-1;i++){                                                                              
   	   for(j=0; j<n-i-1; j++){                                                                            
    		if(a[j]<a[j+1]){                                                                           
     			
		    temp = a[j];                                                               
     		    a[j] = a[j+1];                                                             
     		    a[j+1] = temp;
		                                                                 
    		}                                                                           
   	   }                                                                            
 	}               
	
	 /* calculate median */                                                       
 	if( n%2 == 0)                                                                  
  		med = (a[(n/2)-1]+a[(n/2)])/2.0;                                           
 	else                                                                           
  		med = a[(n/2)];
	
	free_vector(a,0,n-1);
	
	return med;
}
