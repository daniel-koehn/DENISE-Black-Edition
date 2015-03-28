/*------------------------------------------------------------------------
 * send L2 norm of each PE to PE 0
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float exchange_L2(float L2,int sw, int bcast_l2){

        extern int NP, MYID;
	float l2, *L2s, *L2c;
	int i;
	const int tag=1;
	MPI_Status  status;
	
	L2s=vector(1,2);    /* vector for L2 exchange */
	L2c=vector(0,NP+1); /* vector for merging of L2 norms */
	
L2s[1]=L2;
if(MYID==0){L2c[0]=L2;}

/* send L2 norm values to PE 0 */

  for(i=1;i<=NP-1;i++){
      
       /*    if(MYID==i){
           printf("Process %d sending %e to %d\n",i, L2s[1], 0);}*/
        
       
           if(MYID==i){
	     MPI_Send(&L2s[1],2,MPI_FLOAT,0,tag,MPI_COMM_WORLD);}
	     
	     /*MPI_Barrier(MPI_COMM_WORLD);*/
	  
	   if(MYID==0){  
	     MPI_Recv(&L2s[1],2,MPI_FLOAT,i,tag,MPI_COMM_WORLD,&status);
	     L2c[i]=L2s[1];
	     }
            	         
  
  }
  
  /*MPI_Barrier(MPI_COMM_WORLD);*/
  
  if(MYID==0){
/* Merge L2 norms */  
    if(sw==1){
       l2=0.0;
         for(i=0;i<=NP-1;i++){
           l2=l2+L2c[i];
      /*     printf("i=%d \t l2=%e \n",i,l2);
           printf("i=%d \t L2=%e \n",i,L2c[i]);*/
         }
       L2=l2;
    }
    
/* find maximum mun */  
    if(sw==2){
	for(i=0;i<=NP-1;i++){
		if (i == 0) {
			l2 = L2c[i];
	} else {
		if(L2c[i]>l2){
			l2 = L2c[i];
		}
	}
/*	   
           printf("i=%d \t mun_max=%e \n",i,l2);
           printf("i=%d \t mun=%e \n",i,L2c[i]);
 */        }
       L2=l2;
    } 
    
/* find minimum mun */  
    if(sw==3){
	for(i=0;i<=NP-1;i++){
		if (i == 0) {
			l2 = L2c[i];
		} else {
			if ((fabs(L2c[i])<fabs(l2)) && (fabs(L2c[i])>0.0)) {
				l2 = L2c[i];
			}
		}
/*		printf("i=%d \t mun_min=%e \n",i,l2);
		printf("i=%d \t mun=%e \n",i,L2c[i]);
*/	}
	L2=l2;
    } 
        
/* calculate average value */  
    if(sw==4){
       l2 = 0.0;
         for(i=0;i<=NP-1;i++){
	    
           l2 += L2c[i];}
	   l2 /= NP;
           printf("i=%d \t mean(mun)=%e \n",i,l2);
	   L2=l2;
         }
       
/* calculate median */  
    if(sw==5){
    	l2 = median(L2c,NP);
    	printf("i=%d \t median(mun)=%e \n",i,l2);
        L2=l2;
    }
       
    }
        
  
  	
  free_vector(L2s,1,2);		
  free_vector(L2c,0,NP+1);
  
  if (bcast_l2)
  {
  	MPI_Barrier(MPI_COMM_WORLD);
  	MPI_Bcast(&l2,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  	MPI_Barrier(MPI_COMM_WORLD);
  }
  
  return l2;
}

