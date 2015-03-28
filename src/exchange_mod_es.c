/**/
/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 *   last update 29/06/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_mod_es(float ** matmod, int ncptot, int nparameter){


	/* definition of local variables */
	extern int MYID;
	int ntotal = ncptot * nparameter;
	int i,j,h; 
	float fdum1[ntotal];
	
	printf("MYID = %d \t ntotal = %d \n",MYID,ntotal);
	 
	if (MYID == 0){ 
        
	h=1;
	for(i=1;i<=ncptot;i++){
	    for(j=1;j<=nparameter;j++){
	        fdum1[h] = matmod[i][j];
		h++;
	    }
	}                                                                 
        
	} /** if (MYID == 0) **/
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&fdum1,ntotal,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	h=1;
	for(i=1;i<=ncptot;i++){
	    for(j=1;j<=nparameter;j++){
	        matmod[i][j] = fdum1[h] ;
		h++;
	    }
	}                                                                 
        
	
}
