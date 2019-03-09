/* Gradient of the gravity data objective function */
/*------------------------------------------------------------------------
 *
 * Daniel Koehn, Daniel Wehner
 * last update: Kiel, the 24th of May 2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void grav_grad(int ngrav, float **gravpos, float **grad_grav, float *gz_res){

	int i, j, k, ii, jj, HESS_GRAV, p, m, n, l, s, j_sea;
        double time1, time2;
        double G = 6.6738480e-11; /* gravitational constant [N m^2 kg ^{-2}]*/
        double r, x, y, z, r_center, term1, term2, K;
        float ** grad_grav_glob=NULL, HESS_GRAV_EPS, ** hess_grav_glob=NULL, max_hess, * xprism=NULL, * yprism=NULL, * zprism=NULL;
	char jac[225];
        extern float DH;
	extern int MYID, NX, NY, NXG, NYG, POS[3], IDX, IDY, NGRAVB, NZGRAV, GRAV_TYPE;
        extern char JACOBIAN[STRING_SIZE];

	FILE *FP;        

        HESS_GRAV=0;
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
        
	/* grid point below water column --> avoids update within the water column*/
	j_sea = 28;
	
	
	/* invert gravity data */
	
	/* define corner positions of the prisms */
	xprism=vector(1,2);
	yprism=vector(1,2);
	zprism=vector(1,2);
	
	for(p=1;p<=2;p++){
		xprism[p] = 0.0;
		yprism[p] = 0.0;
		zprism[p] = 0.0;
	}
	
	/* loop over all gravimeter stations */
	for(i=1;i<=ngrav;i++){			
		
		/* loop over all subsurface points */
		for(k=1;k<=NXG;k++){
			for(j=j_sea;j<=NYG;j++){
				
				x = (double) (k*DH)+(NGRAVB*DH);
				y = (double) j*DH;
				z = (double) NZGRAV;
				
				/*radius to the center of the prism*/
				r_center = sqrt(pow(gravpos[1][i]-x,(double) 2) + pow(gravpos[2][i]-y,(double) 2));
				
				/* Defining the prism model parameters and its edges */
				xprism[1] = (double) x-(DH/2);
				xprism[2] = (double) x+(DH/2);
				zprism[1] = (double) -z;
				zprism[2] = (double) z;
				yprism[1] = (double) y-(DH/2);
				yprism[2] = (double) y+(DH/2);
				
				/* geometrical kernel K */
				K = 0;
				
				/* radius to the corners of the prism*/
				for(l=1;l<=2;l++){
					for(m=1;m<=2;m++){
						for(n=1;n<=2;n++){
					
							r = sqrt(pow(gravpos[1][i]-xprism[l],(double) 2) + pow(gravpos[2][i]-yprism[n],(double) 2) + pow(zprism[m],(double) 2));
					
							/* modulo operator */
							if(((l+m+n)%2) == 0){
								s = -1;
							}else{
								s = 1;
							}
						
							/* calculate geometrical Kernel K for gravity inversion*/
							if(GRAV_TYPE==1){
						        K += (double) s * ( ((gravpos[1][i]-xprism[l]) * log(0-zprism[m]+r)) + ((0-zprism[m]) * log(gravpos[1][i]-xprism[l]+r)) - ((gravpos[2][i]-yprism[n]) * atan( ((gravpos[1][i]-xprism[l])*(0-zprism[m]))/((gravpos[2][i]-yprism[n])*r))) );
							
							/* calculate Hessian */
							/*if(HESS_GRAV==1){*/
								/* implementation of the calcualtion */
							/*}*/
							}
							
							/* calculate geometrical kernel K for gravity gradient inversion */
							if(GRAV_TYPE==2){
							K += (double) s * ( -atan( ((gravpos[1][i]-xprism[l])*(0-zprism[m])) / ((gravpos[2][i]-yprism[n])*r) ) );
							
							/* calculate Hessian */
							/*if(HESS_GRAV==1){*/
								/* implementation of the calcualtion */
							/*}*/
							}
							
						}
					}
				}
				
				/* Preconditioning of the kernel K by depth weighting function */
				K *= (double) pow((y+0.0005),(double) 0.8);
				
				/* calculate gradient */
				grad_grav_glob[j][k] += (double) (G*K*gz_res[i]);
								
			}
		}
		
		
	}
		
	
	/*----------------------------------------------------------------------------------------------------------------------------*/						
       
       
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
	sprintf(jac,"%s_grav.%i.%i",JACOBIAN,POS[1],POS[2]);
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
