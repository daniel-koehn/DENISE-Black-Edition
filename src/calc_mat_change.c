/*------------------------------------------------------------------------
 *   calculate step length for material parameter update
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void calc_mat_change(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float  **  rhonp1, float **  pi, float **  pinp1, float **  u, float **  unp1, int iter, 
                     int epstest, int INVMAT, float eps_scale_vp, float eps_scale_vs){


	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern float DH, DT;
	extern float EPSILON, MUN, EPSILON_u;
	extern int NX, NY, NXG, NYG,  POS[3], MYID;


	/* local variables */

	float Rho, Vp, Vs, Vsnp1, Vpnp1, x, y, undf, r, pi0, K, mu;
	float dpi, pimax, rhomax, umax,gradmax, gradmax_rho, gradmax_u, epsilon1, pimaxr, gradmaxr, umaxr, gradmaxr_u;
	float VPUPPERLIM, VSLOWERLIM;
	int i, j, ii, jj;
	char modfile[STRING_SIZE];
	
        VPUPPERLIM = 2200.0;
        VSLOWERLIM = 800.0;
	 
        /* find maximum of pi and gradient waveconv */
	if((INVMAT==1)||(INVMAT==0)){
	pimax = 0.0;
	gradmax = 0.0;
	
	    for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		Vp = sqrt((pi[j][i] + 2.0 * u[j][i])/rho[j][i]);
		
		if(Vp>pimax){pimax=Vp;}
		
		if((i*j == 1) || (gradmax == 0.0)) {
				gradmax = fabs(waveconv[j][i]);		
		} else { 
		        if(fabs(waveconv[j][i]) > gradmax){
			   gradmax = fabs(waveconv[j][i]);
			}
		                               
		}
	}}}
	
	
	 /* find maximum of u and gradient waveconv_u */
	if((INVMAT==3)||(INVMAT==0)){
	umax = 0.0;
	gradmax_u = 0.0;
	
	    for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		Vs = sqrt(u[j][i]/rho[j][i]);
		
		if(Vs>umax){umax=Vs;}
		
		if((i*j == 1) || (gradmax_u == 0.0)) {
				gradmax_u = fabs(waveconv_u[j][i]);		
		} else { 
		        if(fabs(waveconv_u[j][i]) > gradmax_u){
			   gradmax_u = fabs(waveconv_u[j][i]);
			}
		                               
		}
	}}}
	
        /* calculate test step length epsilon for pi */
	if(INVMAT==0){
	/*eps_scale_vp = 150.0;*/
	
	pimaxr=0.0; 
	gradmaxr=0.0;
	MPI_Allreduce(&pimax,  &pimaxr,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&gradmax,&gradmaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	EPSILON = (pimaxr/(eps_scale_vp*(gradmaxr)));
	epsilon1 = EPSILON;
	
	MPI_Allreduce(&EPSILON,&epsilon1,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	if (MYID==0)  EPSILON=epsilon1;

	exchange_par();
	
	}
	
	if(INVMAT==0){
	
	umaxr=0.0;
	gradmaxr_u=0.0;
	
	MPI_Allreduce(&umax,&umaxr,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&gradmax_u,&gradmaxr_u,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	EPSILON_u = (umaxr/(eps_scale_vp*gradmaxr_u));
	epsilon1 = EPSILON_u;
	
	MPI_Allreduce(&EPSILON_u,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	if (MYID==0)  EPSILON_u=epsilon1;

	exchange_par();
	
	}
	
	if(epstest==1){	
	/* loop over local grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		    /* update material parameters */
		    
		    if(INVMAT==0){
		     
		      Vs = sqrt(u[j][i]/rho[j][i]);
		      Vp = sqrt((pi[j][i] + 2.0 * u[j][i])/rho[j][i]);
 
		      Vp = Vp + EPSILON*waveconv[j][i];
		      Vs = Vs + EPSILON_u*waveconv_u[j][i];
		      /*Vs = Vp/sqrt(3.0);*/
		      
		      if (Vp>VPUPPERLIM)  Vp = VPUPPERLIM;
                      if (Vs<VSLOWERLIM)  Vs = VSLOWERLIM;
		      
		      /*rho[j][i] = 1000.0*0.31*pow(Vp,(1.0/4.0));*/
		       pi[j][i] = Vp*Vp*rho[j][i] - (2.0 * Vs*Vs*rho[j][i]);
                        u[j][i] = Vs*Vs*rho[j][i];
		      
		    }
		    
		    /*if((INVMAT==2)){
		    rho[j][i] = rho[j][i] - EPSILON_rho*waveconv_rho[j][i];}*/
		    
		    /*if(INVMAT==3){
		    mu = u[j][i];
		    mu = mu + EPSILON_u*waveconv_u[j][i];
		    u[j][i] = mu;}*/
		      
		    /* output of Vp and Vs */
		    unp1[j][i] = sqrt(u[j][i]/rho[j][i]);
		    pinp1[j][i] = sqrt((pi[j][i]+2.0*u[j][i])/rho[j][i]);                        
		}
	}
	 
	sprintf(modfile,"model/waveform_test_model_vp_it_%d.bin",iter);
        writemod(modfile,pinp1,3);
	
	MPI_Barrier(MPI_COMM_WORLD);

        if (MYID==0) mergemod(modfile,3);

        sprintf(modfile,"model/waveform_test_model_vs_it_%d.bin",iter);
        writemod(modfile,unp1,3);
	
	MPI_Barrier(MPI_COMM_WORLD);

        if (MYID==0) mergemod(modfile,3);
	

	sprintf(modfile,"model/waveform_test_model_rho_it_%d.bin",iter);
        writemod(modfile,rho,3);
	
	MPI_Barrier(MPI_COMM_WORLD);

        if (MYID==0) mergemod(modfile,3);
	}	
		
}

