/* Gradient of the gravity data objective function */
/*------------------------------------------------------------------------
 *
 * Daniel Koehn
 * Kiel, the 10th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void grav_grad(int ngrav, float **gravpos, float **grad_grav, float *gz_res){

	int i, j, k, ii, jj, HESS_GRAV;
        double time1, time2;
        double G = 6.6738480e-11; /* gravitational constant [N m^2 kg ^{-2}]*/
        double r, x, y;
        float ** grad_grav_glob=NULL, HESS_GRAV_EPS, ** hess_grav_glob=NULL, max_hess;
	char jac[225];
        extern float DH;
	extern int MYID, NX, NY, NXG, NYG, POS[3], IDX, IDY;
        extern char JACOBIAN[STRING_SIZE];
	FILE *FP;        

        HESS_GRAV=1;
        HESS_GRAV_EPS=0.0005; 

        grad_grav_glob =  matrix(1,NYG,1,NXG); 

        if(HESS_GRAV==1){
          hess_grav_glob = matrix(1,NYG,1,NXG);
          max_hess = 0.0;
        }

	if (MYID==0){
		time1=MPI_Wtime();
		printf("\n **Message from grav_grad (printed by PE %d):\n",MYID);
		printf(" Calculate gradient of the gravity objective function ...");
	}
        
        /* model gravity field by Newtons law */
        for(i=1;i<=ngrav;i++){   /* loop over all gravimeter positions */
    
                /* loop over all subsurface points */
		for(k=1;k<=NXG;k++){       
			for(j=1;j<=NYG;j++){            
                          
        			/* radius vector */
                        	x = (double) k*DH;
				y = (double) j*DH;
  
        			r = sqrt(pow(x-gravpos[1][i],(double) 2) + pow(y-gravpos[2][i],(double) 2));
                        
                        	/* calculate vertical component of the gravity field */
        			grad_grav_glob[j][k] += G*gz_res[i]*(y-gravpos[2][i])/r;

                                /* calculate Hessian */
                                if(HESS_GRAV==1){
           
                                  hess_grav_glob[j][k] += pow((G*(y-gravpos[2][i])/r),(double) 2); 
                                  if(hess_grav_glob[j][k]>max_hess){max_hess=hess_grav_glob[j][k];}                                  

                                }
    
                	}
        	}

       }

       /* distribute gradient results on the different MPI-processes */
	for(i=1;i<=NXG;i++){       
		for(j=1;j<=NYG;j++){            
                  
                        /* global -> local array */
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				
                             ii=i-POS[1]*NX;
			     jj=j-POS[2]*NY;

                             if(HESS_GRAV==0){grad_grav[jj][ii]=grad_grav_glob[j][i];}
                             if(HESS_GRAV==1){grad_grav[jj][ii]=grad_grav_glob[j][i]/(hess_grav_glob[j][i]+HESS_GRAV_EPS*max_hess);}
                                        
			}

        	}
	}

        /* output of the gravity gradients */
	sprintf(jac,"%s_grav.%i%i",JACOBIAN,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
	       fwrite(&grad_grav[j][i],sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge model file */ 
	sprintf(jac,"%s_grav",JACOBIAN);
	if (MYID==0) mergemod(jac,3);


	if (MYID==0){
		time2=MPI_Wtime();
		printf(" finished gravity gradient calculation (real time: %4.2f s).\n",time2-time1);
	}

	/* free memory */
	free_matrix(grad_grav_glob,1,NYG,1,NXG);
        
        if(HESS_GRAV==1){
          free_matrix(hess_grav_glob,1,NYG,1,NXG);
        }
        

}
