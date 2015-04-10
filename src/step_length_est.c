/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search
 *
 * D. Koehn
 * Kiel, 18.1.2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"


float step_length_est(FILE *fprec, float ** waveconv, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, float ** ppi, float ** ppinp1, int iter, int nfstart,
	int nsrc, float ** puipjp, float ** prip, float ** prjp, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvx, float ** pvy, float ** psxx, float ** psyy, float ** psxy, float ** ux, float ** uy, float ** pvxp1, float ** pvyp1, float ** psi_sxx_x, float ** psi_sxy_x,
	float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ** pvxm1, float ** pvym1, float ** uttx, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uxy, float ** uyx, int ntr, int **recpos_loc, float **sectionvx, float **sectionvy, float **sectionp, float **sectioncurl, 
	float **sectiondiv, float **sectionread, int ntr_glob, float ** sectionvxdata, float ** sectionvxdiff, float ** sectionvxdiffold, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, float ** sectionpdata, float ** sectionpdiff, float ** sectionpdiffold, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, 
        float ** bufferrig_to_lef, float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, int itest,int nsrc_glob, int nsrc_loc, MPI_Request * req_send, MPI_Request * req_rec, 
        float ***pr, float ***pp, float ***pq, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, float **ptaup, float **ptaus, 
        float *etajm, float *peta, float *etaip, float **ptausipjp, int **recpos, int *step1, int *step3, float C_vp, float **gradg, float FC, 
	int nxgrav, int nygrav, int ngrav, float **gravpos, float *gz_mod, int NZGRAV, int * recswitch, FILE *FP, int ntr_loc){

	extern int MYID,MIN_ITER,TIME_FILT,STEPMAX, GRAVITY, IDX, IDY, NX, NY, NXG, NYG, POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
        extern float EPS_SCALE, SCALEFAC, LAM_GRAV, GAMMA_GRAV, L2_GRAV_IT1;

	float opteps_vp, ** rho_grav, ** rho_grav_ext;
	int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;

	/* Variables for step length calculation */
        int step2, itests, iteste, stepmax, countstep;
        float scalefac, eps_scale, L2_grav, L2sum1;
	float * gz_res;
	char jac_grav[STRING_SIZE];
	FILE *FP_GRAV;

scalefac = SCALEFAC;  /* scale factor for the step length */
stepmax  = STEPMAX;   /* number of maximum misfit calculations/steplength 2/3*/ 

*step1=0;
step2=0;

/* start with first guess for step length alpha */
eps_scale=EPS_SCALE; /* maximum model change = 1% of the maximum model value */
countstep=0; /* count number of forward calculations */

itests=2;
iteste=2;

while((step2!=1)||(*step1!=1)){

for (itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,sectionpdata,sectionpdiff,sectionpdiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec,pr,pp,pq,fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,ptaup,ptaus,etajm,peta,etaip,ptausipjp,recpos,FC,recswitch,FP,ntr_loc);

if(GRAVITY==2){

  /* save seismic L2-norm of seismic data residuals */
  L2sum1 = L2t[itest];
  gz_res = vector(1,ngrav);
  
  /* global density model */
  rho_grav =  matrix(1,NYG,1,NXG);
  rho_grav_ext =  matrix(1,nygrav,1,nxgrav);

  /* model gravity data */
  /* save current density model */
  sprintf(jac_grav,"%s_tmp.rho.%i%i",JACOBIAN,POS[1],POS[2]);
  FP_GRAV=fopen(jac_grav,"wb");

  for (i=1;i<=NX;i=i+IDX){
      for (j=1;j<=NY;j=j+IDY){
          fwrite(&prhonp1[j][i],sizeof(float),1,FP_GRAV);
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
  LAM_GRAV = GAMMA_GRAV * (L2sum1/L2_GRAV_IT1); 

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
  if(TIME_FILT==0){
  err(" ");
  } 
  *step3=1;
  break;
}

if((*step1==1)&&(countstep>stepmax)){
  if(MYID==0){
  printf("Could not found a proper 3rd step length which brackets the minimum\n");}
  *step1=1;
  step2=1;
}

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",iteste,itests,*step1,step2,eps_scale,countstep,stepmax,scalefac,MYID,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

if(*step1==1){ /* only find an optimal step length if step1==1 */
/* calculate optimal step length epsilon for Vp and Vs*/
if(MYID==0){
printf("================================================= \n");
printf("calculate optimal step length epsilon for Vp and Vs \n");
printf("================================================= \n");
}
opteps_vp=calc_opt_step(L2t,waveconv,gradg,epst1,1,C_vp);
eps_scale = opteps_vp;
}

return eps_scale;
}

