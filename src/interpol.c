/*------------------------------------------------------------------------
 * Interpolation of unknown wavefield variables at the cfgt
 *   last update 10/01/05, D. Koehn
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void interpol(int IDXI, int IDYI, float **  intvar, int cfgt_check){


	int i,j;
	extern int NX, NY, POS[3], NPROCX, NPROCY, INDEX[5], BOUNDARY;
	MPI_Status  status;	
        float * int_top_to_bot, * int_lef_to_rig;
        
        int_top_to_bot=vector(1,NX);
        int_lef_to_rig=vector(1,NY);
        
/* exchange gridpoint NY+1 */
  if (POS[2]!=0){  /* no boundary exchange at top of global grid */
     for (i=1;i<=NX;i++){int_top_to_bot[i] = intvar[1][i];}
  }          
                                                                                                                                                                       
  MPI_Bsend(&int_top_to_bot[1],NX,MPI_FLOAT,INDEX[3],3,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Recv(&int_top_to_bot[1],NX,MPI_FLOAT,INDEX[4],3,MPI_COMM_WORLD,&status);
                                                                              
  if (POS[2]!=NPROCY-1){   /* no boundary exchange at bottom of global grid */
     for (i=1;i<=NX;i++){intvar[NY+1][i]=int_top_to_bot[i];}
  }   
                                                                                                                                                        
 /* interpolate in y-direction */	
   for (i=1;i<=NX;i=i+IDXI){ 
	for (j=2;j<=NY;j=j+IDYI){
	
	    intvar[j][i] = (intvar[j-1][i]+intvar[j+1][i])/2.0; 
          
	}
   }

/* exchange gridpoint NX+1 */
  if ((BOUNDARY) || (POS[1]!=0)){  /* no boundary exchange at top of global grid */
       for (j=1;j<=NY;j++){int_lef_to_rig[j] = intvar[j][1];}
  }
  
  MPI_Bsend(&int_lef_to_rig[1],NY,MPI_FLOAT,INDEX[1],3,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Recv(&int_lef_to_rig[1],NY,MPI_FLOAT,INDEX[2],3,MPI_COMM_WORLD,&status);
  
  if ((BOUNDARY) || (POS[1]!=NPROCX-1)){   /* no boundary exchange at bottom of global grid */
      for (j=1;j<=NY;j++){intvar[j][NX+1]=int_lef_to_rig[j];}
  }                
             
/* interpolate in x-direction */
   for (i=2;i<=NX;i=i+IDXI){ 
	for (j=1;j<=NY;j++){
	
	    intvar[j][i] = (intvar[j][i-1]+intvar[j][i+1])/2.0;
	    /*if(i==NX){intvar[j][i]=intvar[j][i-1];}*/
	    
          
	}
   }

free_vector(int_top_to_bot,1,NX);
free_vector(int_lef_to_rig,1,NY);
   
   
}
