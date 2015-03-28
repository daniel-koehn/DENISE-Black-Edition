/* Model gravity field data */
/*------------------------------------------------------------------------
 * model z-component of the gravity field by Newtons law
 *
 * Daniel Koehn
 * Kiel, the 7th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void grav_mod(float  **rho, int ngrav, float **gravpos, float *gz, int NXGRAV, int NYGRAV, int NZGRAV){

	int i, j, k, l, h, h1, n1, n2;
        double time1, time2;
        double G = 6.6738480e-11; /* gravitational constant [N m^2 kg ^{-2}]*/
        double r, DH3, x, y, z, zpos;
        double gztmp, gzall;
	extern float DH;
	extern int MYID, NP;
        extern char GRAV_DATA_OUT[STRING_SIZE];
	FILE *FP;

        DH3 = DH * DH * DH;        
 
	if (MYID==0){
		time1=MPI_Wtime();
		printf("\n **Message from grav_mod (printed by PE %d):\n",MYID);
		printf(" Calculate vertical component of the gravity field ...");
	}

        /* z-position of the gravity stations */
        zpos = (NZGRAV/2) * DH;

        /* parallelization of sum in z-direction */
        /* calculate indices for summation on each core */
	n1 = (NZGRAV / NP) * MYID;
	if (NZGRAV % NP > MYID){
	   n1 += MYID;
           n2 = n1 + (NZGRAV / NP) + 1;
	}else{
	   n1 += NZGRAV % NP;
	   n2 = n1 + (NZGRAV / NP);
	}
        
        /* model gravity field by Newtons law */
        for(i=1;i<=ngrav;i++){           /* loop over all gravimeter positions */

            gztmp = 0;
    
                /* loop over all subsurface points */
                for(l=n1;l<n2;l++){
    			for(k=1;k<=NXGRAV;k++){       
        			for(j=1;j<=NYGRAV;j++){            
                                  
                			/* radius vector */
                                	x = (double) k*DH;
					y = (double) j*DH;
                                        z = (double) l*DH;

                			r = sqrt(pow(x-gravpos[1][i],(double) 2) + pow(y-gravpos[2][i],(double) 2) + pow(z-zpos,(double) 2));
                                
                                	/* calculate vertical component of the gravity field */
                			gztmp += (double)(DH3*G*rho[j][k]*(y-gravpos[2][i])/pow(r,(double) 3.0));
            
                        	}
                	}

		}

                /* collect results from all cores and calculate the sum*/
                MPI_Allreduce(&gztmp,&gzall,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		/*printf("%d \t %d \t %d \t %e \n",n1,n2,MYID,gzall);*/

		gz[i] = gzall;

       }

       /* output of the gravity modelling results */
       if(MYID==0){

         FP=fopen(GRAV_DATA_OUT,"w");         

         for(i=1;i<=ngrav;i++){
         	fprintf(FP,"%d \t %e \n",i,gz[i]);
         }

         fclose(FP);

       }

	if (MYID==0){
		time2=MPI_Wtime();
		printf(" finished gravity modelling (real time: %4.2f s).\n",time2-time1);
	}
}
