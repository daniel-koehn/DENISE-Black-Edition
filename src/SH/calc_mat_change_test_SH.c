/*------------------------------------------------------------------------
 *   calculate test step length for material parameter update
 *   
 *   Daniel Koehn
 *   last update 22.12.2017
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
float calc_mat_change_test_SH(float  **  waveconv_rho, float  **  waveconv_u, float  **  waveconv_ts, float  **  rho, float  **  rhonp1, float **  u, float **  unp1, float **  ts, float **  tsnp1, int iter, int epstest, float eps_scale, int itest){


	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern float DH, DT;
	extern float EPSILON, EPSILON_u, EPSILON_rho, EPSILON_ts, MUN;
	extern float C_vs, C_rho, C_taus, C_vs_min, C_rho_min, C_taus_min;
	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	
	extern int INV_RHO_ITER, INV_VS_ITER, INV_QS_ITER;
	
	extern char INV_MODELFILE[STRING_SIZE];
	
	extern float VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM, QSUPPERLIM, QSLOWERLIM;

	/* local variables */

	float Rho, Vs, Vsnp1, x, y, undf, r, K, mu, Zs, eps_true;
	float rhomax, umax, tausmax, gradmax_rho, gradmax_u, gradmax_ts, epsilon1, gradmaxr_u, umaxr;
	float gradmaxr_rho, rhomaxr, gradmaxr_ts, tausmaxr;
	int i, j, ii, jj, testuplow;
	char modfile[STRING_SIZE];		
	
	/* find maximum of Vs and gradient waveconv_u */
	

	umax = 0.0;
	gradmax_u = 0.0;
	
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		//Zs = (u[j][i]-C_vs_min) / C_vs;     /* normalize Vs */
		Zs = u[j][i];		

		if(Zs>umax){umax=Zs;}
		
		if((i*j == 1) || (gradmax_u == 0.0)) {
		        gradmax_u = fabs(waveconv_u[j][i]);
		} else {
		if(fabs(waveconv_u[j][i]) > gradmax_u)
		{
		gradmax_u = fabs(waveconv_u[j][i]);
					
				}
		                               
		}
	}}	
	
	/* find maximum of rho and gradient waveconv_rho */
	rhomax = 0.0;
	gradmax_rho = 0.0;
	
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		//if((rho[j][i]-C_rho_min)/C_rho>rhomax){rhomax = (rho[j][i]-C_rho_min)/C_rho;}
		if(rho[j][i]>rhomax){rhomax = rho[j][i];}		

		if((i*j == 1) || (gradmax_rho == 0.0)) {
		        gradmax_rho = fabs(waveconv_rho[j][i]);
		} else {
		if(fabs(waveconv_rho[j][i]) > gradmax_rho)
		{
		gradmax_rho = fabs(waveconv_rho[j][i]);
					}
				}
		                               
		}}

	/* find maximum of taus and gradient waveconv_ts */
	/*tausmax = 0.0;
	gradmax_ts = 0.0;
	
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		if((ts[j][i]-C_taus_min)/C_taus>tausmax){tausmax=(ts[j][i]-C_taus_min)/C_taus;}
		
		if((i*j == 1) || (gradmax_ts == 0.0)) {
		        gradmax_ts = fabs(waveconv_ts[j][i]);
		} else {
		if(fabs(waveconv_ts[j][i]) > gradmax_ts)
		{
		gradmax_ts = fabs(waveconv_ts[j][i]);
					}
				}
		                               
		}}*/

		
	/* calculate scaling factor for the gradients */
        /* --------------------------------------------- */	
	
	/* parameter 1 */
	   MPI_Allreduce(&umax,&umaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax_u,&gradmaxr_u,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	   EPSILON_u = eps_scale * (umaxr/gradmaxr_u);
	   if (iter<INV_VS_ITER){EPSILON_u = 0.0;}
           epsilon1=EPSILON_u;
	   MPI_Allreduce(&EPSILON_u,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	   EPSILON_u=epsilon1;		
	
	
	/* parameter 2 */

	   MPI_Allreduce(&rhomax,&rhomaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax_rho,&gradmaxr_rho,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
    	   EPSILON_rho = eps_scale * (rhomaxr/gradmaxr_rho);
	   if (iter<INV_RHO_ITER){EPSILON_rho = 0.0;}
           epsilon1=EPSILON_rho;
	   MPI_Allreduce(&EPSILON_rho,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	   EPSILON_rho=epsilon1;	

	/* parameter 3 */

	   /*MPI_Allreduce(&tausmax,&tausmaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax_ts,&gradmaxr_ts,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
    	   EPSILON_ts = eps_scale * (tausmaxr/gradmaxr_ts);
	   if (iter<INV_QS_ITER){EPSILON_ts = 0.0;}
           epsilon1=EPSILON_ts;
	   MPI_Allreduce(&EPSILON_ts,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	   EPSILON_ts=epsilon1;	*/

	   
        if(MYID==0){
	  printf("MYID = %d \t umaxr = %e \t gradmaxr_u = %e \n",MYID,umaxr,gradmaxr_u);
	  printf("MYID = %d \t rhomaxr = %e \t gradmaxr_rho = %e \n",MYID,rhomaxr,gradmaxr_rho);
	  //printf("MYID = %d \t tsmaxr = %e \t gradmaxr_ts = %e \n",MYID,tausmaxr,gradmaxr_ts);
	}
      
      /* save true step length */
      eps_true = EPSILON_u;
	
	/* loop over local grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){		
		    
		    /* update normalized material parameters Vs, rho, taus */
		    if((INVMAT1==3) || (INVMAT1==1)){
		      
		      testuplow=0;
		                         
		        /*unp1[j][i] = ((u[j][i]-C_vs_min) / C_vs) - EPSILON_u * waveconv_u[j][i];   	
		      rhonp1[j][i] = ((rho[j][i]-C_rho_min) / C_rho) - EPSILON_rho * waveconv_rho[j][i];
		       tsnp1[j][i] = ((ts[j][i]-C_taus_min) / C_taus) - EPSILON_ts * waveconv_ts[j][i];*/
		      
		        unp1[j][i] = u[j][i] - EPSILON_u * waveconv_u[j][i];   	
		      rhonp1[j][i] = rho[j][i] - EPSILON_rho * waveconv_rho[j][i];

			
		      /* apply bound constraints */      
                      if(INVMAT1==1){
		      
		        /* Denormalize updated material parameters */
		        /*unp1[j][i] = C_vs_min + unp1[j][i] * C_vs;
		        rhonp1[j][i] = C_rho_min + rhonp1[j][i] * C_rho;			  
		        tsnp1[j][i] = C_taus_min + tsnp1[j][i] * C_taus;*/	
		      		  		      		      
			/* S-wave velocities */
		        if((unp1[j][i]<VSLOWERLIM)&&(unp1[j][i]>1e-6)){
		        	unp1[j][i] = u[j][i];
			}

			if(unp1[j][i]>VSUPPERLIM){
		        	unp1[j][i] = u[j][i];
			}
		      
		      	/* densities */
			if(rhonp1[j][i]<RHOLOWERLIM){
		        	rhonp1[j][i] = rho[j][i];
			}
			
			if(rhonp1[j][i]>RHOUPPERLIM){
		        	rhonp1[j][i] = rho[j][i];
			}

			/* taus */
			/*if(tsnp1[j][i]>2.0/QSLOWERLIM){
		        	tsnp1[j][i] = ts[j][i];
			}
			
			if(tsnp1[j][i]<2.0/QSUPPERLIM){
		        	tsnp1[j][i] = ts[j][i];
			}*/

		      }
		      
		      /* None of these parameters should be smaller than zero */
		      if(unp1[j][i]<0.0){
		        unp1[j][i] = u[j][i];
		      } 

		      if(rhonp1[j][i]<0.0){
		        rhonp1[j][i] = rho[j][i];
		      }

		      /*if(tsnp1[j][i]<0.0){
		        tsnp1[j][i] = ts[j][i];
		      }	*/	      		  
					  
		      if(itest==0){
			  rho[j][i] = rhonp1[j][i];
                            u[j][i] = unp1[j][i];
                            /*ts[j][i] = tsnp1[j][i];*/
			   
		      } 
		      
		    }
                              
		}
	}
	
	if(itest==0){

	   sprintf(modfile,"%s_vs.bin",INV_MODELFILE);
	
           writemod(modfile,unp1,3);
	
	   MPI_Barrier(MPI_COMM_WORLD);

           if (MYID==0) mergemod(modfile,3);
	
	   sprintf(modfile,"%s_rho.bin",INV_MODELFILE);
	   writemod(modfile,rho,3);
	
	   MPI_Barrier(MPI_COMM_WORLD);

           if (MYID==0) mergemod(modfile,3);

	   /*sprintf(modfile,"%s_ts.bin",INV_MODELFILE);
	   writemod(modfile,ts,3);
	
	   MPI_Barrier(MPI_COMM_WORLD);

           if (MYID==0) mergemod(modfile,3);*/

        }

	return eps_true;
}

