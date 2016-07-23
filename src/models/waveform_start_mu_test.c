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
	const float r1 = 10.0;
	const float r2 = 5.5;
	

        /*FP1=fopen("/stripe1/koehn/projects/fullwaveform/full_waveform/test_forward/par/model/vp.bin","rb");*/
	        
	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){
		
                                x=(float)i*DH;
				y=(float)j*DH;
				
				/*fread(&Vpnm1,sizeof(float),1,FP1);*/
				
				Vp=vp2; Vs=vs2; Rho=rho2;		
		                
				r=sqrt(((x-X01)*(x-X01))+((y-Y01)*(y-Y01)));
				if(r<=r1){
				  Vp=vp3;
				}
			   			
			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				u[jj][ii]=Vs*Vs*Rho;
				rho[jj][ii]=Rho;
				pi[jj][ii]=Vp*Vp*Rho - (4.0/3.0)*u[jj][ii];
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

