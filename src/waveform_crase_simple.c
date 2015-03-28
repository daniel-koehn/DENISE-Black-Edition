/*------------------------------------------------------------------------
 *   Jastram-Test-Model for Adaptive FD-Grid
 *   
 *   Daniel Koehn
 *   last update 23.10.2004
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void model_elastic(float  **  rho, float **  pi, float **  u){


	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern float DH;
	extern int NX, NY, NXG, NYG,  POS[3], MYID;


	/* local variables */

	float Rho, Vp, Vs, Vpnm1, x, y, undf, r;
	float aund, ampund, FW, shiftx;
	int i, j, ii, jj;
	char modfile[STRING_SIZE];
	
        /* parameters for background */
	const float vp2=2000.0, vs2=vp2/sqrt(3.0), rho2=1000.0*0.31*pow(vp2,(1.0/4.0));
	
	/* parameters for sphere 1 and 2 */
	const float vp3=1500.0, vs3=vp3/sqrt(3.0), rho3=1000.0*0.31*pow(vp3,(1.0/4.0));
	
	/* location of the spheres */
	const float X01 = 80.0;
	const float Y01 = 130.0;
	
	/* radii of spheres */
        float A0, A1, A3, A4, lambda0, lambda1, lambda3, lambda4;
	float y0, y1, y2, y3, y4, y5, undy0, undy1, undy3, undy4;

	 
	 lambda0=1600.0;
              A0=50.0;
              y0=1200.0;
 
         lambda1=lambda0/2.0;
              A1=100.0;
              y1=1000.0;
 
              y2=800.0;
 
         lambda3=lambda0/2.0;
              A3=150.0;
              y3=600.0;
 
         lambda4=lambda0/2.0;
              A4=50.0;
              y4=410.0;
 
              y5=100.0;


        /*FP1=fopen("/stripe1/koehn/projects/fullwaveform/full_waveform/test_forward/par/model/vp.bin","rb");*/
	        
	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){
		
                        x=(float)i*DH;
		        y=(float)j*DH;
				
		        if(y<=NXG*DH){  
                          Vp = 2178.0;
                          Vs = Vp/sqrt(3.0);
                         Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
   
     
                        undy0 = A0*sin(2.0*PI*x/lambda0); 
                        if(y<=(y0+undy0)){  
                          Vp  = 1980.0;
                          Vs  = Vp/sqrt(3.0);
                          Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
     
                        undy1 = A1*cos(2.0*PI*x/lambda1); 
                        if(y<=(y1+undy1)){  
                          Vp = 1980.0;
                          Vs = Vp/sqrt(3.0);
                         Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
    
                        if(y<=y2){  
                          Vp = 1980.0;
                          Vs = Vp/sqrt(3.0);
                         Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
     
                        undy3 = A3*sin(2.0*PI*x/lambda3); 
                        if(y<=(y3+undy3)){  
                          Vp = 1800.0;
                          Vs = Vp/sqrt(3.0);
                         Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
     
                        undy4 = A4*sin(2.0*PI*x/lambda4); 
                        if(y<=(y4+undy4)){  
                          Vp = 1800.0;
                          Vs = Vp/sqrt(3.0);
                         Rho = 1000.0*0.31*pow(Vp,(1.0/4.0));}
     
                        if(y<=y5){  
                          Vp = 1500.0;
                          Vs = 1e-6;
                         Rho = 1000.0;}

				
			   			
			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				u[jj][ii]=Vs*Vs*Rho;
				rho[jj][ii]=Rho;
				pi[jj][ii] = Vp*Vp*Rho - 2.0 * u[jj][ii];
			}
		}
	}	

		

	
	/* each PE writes his model to disk */
        sprintf(modfile,"model/waveform_test_model_u.bin");
        writemod(modfile,u,3);
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
	
	
        sprintf(modfile,"model/waveform_test_model_pi.bin");
        writemod(modfile,pi,3);
	

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); 
	
	sprintf(modfile,"model/waveform_test_model_rho.bin");
        writemod(modfile,rho,3);
	

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
	
	/*fclose(FP1);*/
}

