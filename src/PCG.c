/*------------------------------------------------------------------------
 * Module for the Preconditioned Conjugate Gradient Method (PCG)
 * 
 * Daniel Koehn
 * Kiel, 12.08.2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void PCG(float * PCG_new, float * PCG_old, float * PCG_dir, int PCG_class){

	extern int NX, NY, IDX, IDY, SPATFILTER, GRAD_FILTER;
	extern int MYID, PCG_BETA;
	int i, j, k, h;
	float betaz, betan, tmp, *beta;

	beta  =  vector(1,PCG_class);	

	/* calculate CG-beta */
	/* ----------------- */

	/* loop over parameter classes */
        h=1;
	for (k=1;k<=PCG_class;k++){   

	    /* calculate beta */
	    betaz = 0.0;
	    betan = 0.0;

	    for (i=1;i<=NX;i=i+IDX){
	        for (j=1;j<=NY;j=j+IDY){	  	  

		    /* Fletcher-Reeves */
		    if(PCG_BETA==1){
	               betaz += PCG_new[h] * PCG_new[h];
	               betan += PCG_old[h] * PCG_old[h];
		    }	  

	            /* Polak-Ribiere */
		    if(PCG_BETA==2){
	               betaz += PCG_new[h] * (PCG_new[h] - PCG_old[h]);
	               betan += PCG_old[h] * PCG_old[h];
		    }

	            /* Hestenes-Stiefel */
		    if(PCG_BETA==3){
	               betaz += - PCG_new[h] * (PCG_new[h] - PCG_old[h]);
	               betan += PCG_dir[h] * (PCG_new[h] - PCG_old[h]);
		    }

		    /* Dai-Yuan */
		    if(PCG_BETA==4){
	               betaz += - PCG_new[h] * PCG_new[h];
	               betan += PCG_dir[h] * (PCG_new[h] - PCG_old[h]);
		    }	  

	  
                    h++;
	  	  
	        }
	    }

	    tmp = 0.0;
	    MPI_Allreduce(&betaz,&tmp,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	    betaz = tmp;
	     
	    tmp = 0.0;
	    MPI_Allreduce(&betan,&tmp,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	    betan = tmp;
	     
	    beta[k] = 0.0f;
	    if(betan !=0.0f) beta[k] = betaz/betan;
	     
	    /* direction reset */
	    if((beta[k]<0.0)&&(PCG_BETA==2)){beta[k] = 0.0;}

            if(MYID==0){
	        printf("\n\n parameter class = %d \t beta = %e\n",k,beta[k]);
	    } 

        }
	     
	/* calculate conjugate gradient direction */
	/* -------------------------------------- */

	/* loop over parameter classes */
        h=1;
	for (k=1;k<=PCG_class;k++){   

	    for (i=1;i<=NX;i=i+IDX){
	        for (j=1;j<=NY;j=j+IDY){	  	  
	  
		    PCG_dir[h] = PCG_new[h] + beta[k] * PCG_dir[h]; 
	  
                    h++;
	  	  
	        }
	    }

        }
	          
	/* dealllocate memory */
	free_vector(beta,1,PCG_class);

}
