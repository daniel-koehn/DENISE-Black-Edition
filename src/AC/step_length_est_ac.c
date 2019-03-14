/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search for the AC case
 *
 * D. Koehn
 * Kiel, 12.06.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"


float step_length_est_ac(struct waveAC *waveAC, struct waveAC_PML *waveAC_PML, struct matAC *matAC, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int iter, int nsrc, int ns, int ntr, int ntr_glob, float * epst1, 
         double * L2t, int nsrc_glob, int nsrc_loc, int *step1, int *step3, int nxgrav, int nygrav, int ngrav, float **gravpos, float *gz_mod, int NZGRAV, int ntr_loc, 
         float **Ws, float **Wr, int hin, int *DTINV_help, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern int MYID,MIN_ITER,TIME_FILT,STEPMAX, GRAVITY, IDX, IDY, NX, NY, NXG, NYG, POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
        extern float EPS_SCALE, SCALEFAC, LAM_GRAV, GAMMA_GRAV, L2_GRAV_IT1;
        extern float FC;

        /* local variables */
	float opteps_vp, ** rho_grav, ** rho_grav_ext;
	int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;
        int step2, itest, itests, iteste, stepmax, countstep;
        float scalefac, eps_scale, L2_grav, tmp, tmp1;
	double L2sum1;
	float * gz_res;
	char jac_grav[STRING_SIZE];

	FILE *FP_GRAV, *FP;

	scalefac = SCALEFAC;  /* scale factor for the step length */
	stepmax  = STEPMAX;   /* number of maximum misfit calculations/steplength 2/3*/ 

	*step1=0;
	step2=0;

	/* start with first guess for step length alpha */
	eps_scale=EPS_SCALE; /* maximum model change = 1% of the maximum model value */
	countstep=0;         /* count number of forward calculations */

	itests=2;
	iteste=2;

        /* store current AC models */
        copy_mat((*matAC).prho,(*fwiPSV).prho_old);
        copy_mat((*matAC).ppi,(*fwiPSV).ppi_old);

	while((step2!=1)||(*step1!=1)){

	for (itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

        /* update material parameters for test step eps_scale */
	tmp=calc_mat_change_test_AC((*fwiPSV).waveconv,(*fwiPSV).waveconv_rho,(*fwiPSV).prho_old,(*matAC).prho,(*fwiPSV).ppi_old,(*matAC).ppi,iter,1,eps_scale,1);

        (*seisPSVfwi).L2 = obj_ac(waveAC,waveAC_PML,matAC,fwiPSV,mpiPSV,seisPSV,seisPSVfwi,acq,hc,nsrc,nsrc_loc,nsrc_glob,ntr,ntr_glob,ns,itest,iter,Ws,Wr,hin,DTINV_help,eps_scale,req_send,req_rec);
        L2t[itest] = (*seisPSVfwi).L2;

        epst1[itest]=eps_scale;
	epst1[1] = 0.0;

	if(GRAVITY==2){

	  /* save seismic L2-norm of seismic data residuals */
	  L2sum1 = L2t[itest];
	  gz_res = vector(1,ngrav);
	  
	  /* global density model */
	  rho_grav =  matrix(1,NYG,1,NXG);
	  rho_grav_ext =  matrix(1,nygrav,1,nxgrav);

	  /* model gravity data */
	  /* save current density model */
	  sprintf(jac_grav,"%s_tmp.rho.%i.%i",JACOBIAN,POS[1],POS[2]);
	  FP_GRAV=fopen(jac_grav,"wb");

	  for (i=1;i<=NX;i=i+IDX){
	      for (j=1;j<=NY;j=j+IDY){
                  tmp1 = (*matAC).prho[j][i];
		  fwrite(&tmp1,sizeof(float),1,FP_GRAV);
	      }
	  }
	
	  fclose(FP_GRAV);

	  MPI_Barrier(MPI_COMM_WORLD);
		  
	  /* merge model file */ 
	  sprintf(jac_grav,"%s_tmp.rho",JACOBIAN);
	  if (MYID==0) mergemod(jac_grav,3);
	  
	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  read_density_glob(rho_grav,2);
	  extend_mod(rho_grav,rho_grav_ext,nxgrav,nygrav);
	  grav_mod(rho_grav_ext,ngrav,gravpos,gz_mod,nxgrav,nygrav,NZGRAV);

	  /* calculate gravity data residuals */
	  L2_grav=calc_res_grav(ngrav,gz_mod,gz_res);

	  /* calculate lambda_grav */
	  /* LAM_GRAV = GAMMA_GRAV * (L2sum1/L2_GRAV_IT1); */

	  /* add gravity penalty term to the seismic objective function */
	  L2t[itest]+=LAM_GRAV * L2_grav;

	  /* free memory */
	  free_matrix(rho_grav,1,NYG,1,NXG);
	  free_matrix(rho_grav_ext,1,nygrav,1,nxgrav);
	  free_vector(gz_res,1,ngrav);

	}
	     
	} /* end of L2 test */

	/* Did not found a step size which reduces the misfit function */
	if((*step1==0)&&(L2t[1]<=L2t[2])){
	 eps_scale = eps_scale/scalefac; 
	 countstep++;
	}

	/* Found a step size with L2t[2] < L2t[3]*/
	if((*step1==1)&&(L2t[2]<L2t[3])){
	 epst1[3]=eps_scale;
	 step2=1;
	}

	/* Could not found a step size with L2t[2] < L2t[3]*/
	if((*step1==1)&&(L2t[2]>=L2t[3])){
	 epst1[3]=eps_scale;
	 /* increase step length to find  a larger misfit function than L2t[2]*/
	 eps_scale = eps_scale + (eps_scale/scalefac);
	 countstep++;                       
	}         

	/* found a step size which reduces the misfit function */
	if((*step1==0)&&(L2t[1]>L2t[2])){
	 epst1[2]=eps_scale; 
	 *step1=1;
	 iteste=3;
	 itests=3;
	 countstep=0;
	 /* find a second step length with a larger misfit function than L2t[2]*/
	 eps_scale = eps_scale + (eps_scale/scalefac);
	}

	*step3=0;

	if((*step1==0)&&(countstep>stepmax)){
	  if(MYID==0){
	  printf(" Steplength estimation failed!");} 
	  *step3=1;
	  break;
	}

	if((*step1==1)&&(countstep>stepmax)){
	  if(MYID==0){
	  printf("Could not found a proper 3rd step length which brackets the minimum\n");}
	  *step1=1;
	  step2=1;
	}

	if(MYID==0){
        printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",iteste,itests,*step1,step2,eps_scale,countstep,stepmax,scalefac,MYID,L2t[1],L2t[2],L2t[3]);}

	} /* end of while loop */

	if(*step1==1){ /* only find an optimal step length if step1==1 */

		if(MYID==0){
		   printf("================================================= \n");
		   printf("calculate optimal step length epsilon for Vp and Vs \n");
		   printf("================================================= \n");
		}

		opteps_vp=calc_opt_step(L2t,epst1,1);
		eps_scale = opteps_vp;

	}

return eps_scale;
}

