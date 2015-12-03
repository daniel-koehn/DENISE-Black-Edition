/* Model gravity field data */
/*------------------------------------------------------------------------
 * model z-component of the gravity field by Newtons law
 *
 * Daniel Koehn, Daniel Wehner
 * last update: Kiel, the 13th of October 2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void grav_mod(float  **rho, int ngrav, float **gravpos, float *gz_mod, int NXGRAV, int NYGRAV, int NZGRAV){

	int i, j, k, l, h, h1, n1, n2, s, m, n, p;
        double time1, time2;
        double G = 6.6738480e-11; /* gravitational constant [N m^2 kg ^{-2}]*/
        double r, DH3, x, y, z, zpos;
        double gztmp, gzall;
	float * rho_mean=NULL, * xprism=NULL, * yprism=NULL, * zprism=NULL;
	extern float DH;
	extern int MYID, NP, NYG, NGRAVB, GRAV_TYPE, NXG;
        extern char GRAV_DATA_OUT[STRING_SIZE];
	float ** rho_back=NULL, ** rho_back_ext=NULL;
	
	FILE *FP;

        DH3 = DH * DH * DH;
	
	/* reset gz_mod */
	for(i=1;i<=ngrav;i++){
		gz_mod[i] = 0.0;
	}
	 
	       
	if (MYID==0){
		time1=MPI_Wtime();
		printf("\n **Message from grav_mod (printed by PE %d):\n",MYID);
		printf(" Calculate vertical component of the gravity field ... \n");
	}
	
        
	/* load background density model */
	rho_back = matrix(1,NYG,1,NXG);
  	rho_back_ext = matrix(1,NYG,1,NXGRAV);
  	read_back_density(rho_back);
  	extend_mod(rho_back,rho_back_ext,NXGRAV,NYG);
		
	
	/* model gravity field by Newtons law for rectangular prisms */
	
	/* define corner positions of the prisms */
	xprism=vector(1,2);
	yprism=vector(1,2);
	zprism=vector(1,2);
	
	for(p=1;p<=2;p++){
		xprism[p] = 0.0;
		yprism[p] = 0.0;
		zprism[p] = 0.0;
	}
	
	/* loop over all gravimeter positions*/
	for (i=1;i<=ngrav;i++){
	
	gztmp = 0.0;
	
		/* loop over all subsurface points */
		for(k=1;k<=NXGRAV;k++){
			for(j=1;j<=NYG;j++){
			
				x = (double) k*DH;
				y = (double) j*DH;
				z = (double) NZGRAV;
			
				/* Defining the prism model parameters and its edges */
				xprism[1] = (double) x-(DH/2);
				xprism[2] = (double) x+(DH/2);
				zprism[1] = (double) -z;
				zprism[2] = (double) z;
				yprism[1] = (double) y-(DH/2);
				yprism[2] = (double) y+(DH/2);
			
			
				/* radius to the corners of the prism */
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
											
						/* calculate gz component (Nagy, 2000)*/
						if(GRAV_TYPE==1){
						gztmp += (double)(s*G*(rho[j][k]-rho_back_ext[j][k])) * ( ((gravpos[1][i]-xprism[l]) * log(0-zprism[m]+r)) + ((0-zprism[m]) * log(gravpos[1][i]-xprism[l]+r)) - ((gravpos[2][i]-yprism[n]) * atan( ((gravpos[1][i]-xprism[l])*(0-zprism[m]))/((gravpos[2][i]-yprism[n])*r))) );
						}
						
						/* calculate gzz component (by Nagy, 2000)*/
						if(GRAV_TYPE==2){
						gztmp += (double)(s*G*(rho[j][k]-rho_back_ext[j][k])) * ( -atan( ((gravpos[1][i]-xprism[l])*(0-zprism[m])) / ((gravpos[2][i]-yprism[n])*r) ) );
						}
																
						}
					}
				}
			
			
			}
		}
		
		gz_mod [i] = gztmp;
	
	} 
	
        /*------------------------------------------------------------------------------------------------------------------*/

       	/* output of the gravity modelling results */
	if(MYID==0){

       		FP=fopen(GRAV_DATA_OUT,"w");         

         	for(i=1;i<=ngrav;i++){
         		fprintf(FP,"%d \t %e \n",i,gz_mod[i]);
         	}

         	fclose(FP);

       	}

	if (MYID==0){
		time2=MPI_Wtime();
		printf(" finished gravity modelling (real time: %4.2f s).\n",time2-time1);
	}
	
	/* free memory */
	free_matrix(rho_back,1,NYG,1,NXG);
  	free_matrix(rho_back_ext,1,NYGRAV,1,NXGRAV);
}
