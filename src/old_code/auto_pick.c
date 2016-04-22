/*------------------------------------------------------------------------
 *  Auto_Pick: STA/LTA Autopicker 
 *
 *  Beta 0.1 release
 *
 *  Copyright (c) D. Koehn 
 *  
 *  July, the 31st 2012
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */
#include "cseife.h"

#include "stfinv/stfinv.h" /* libstfinv - inversion for source time function */

int main(int argc, char **argv){
/* variables in main */
int ns, nseismograms=0, nt, nd, fdo3, j, i, ii, jj, shotid, recid, k, nc, iter, h, infoout, SHOTINC, TIMEWIN, test_eps, lq, iq, jq, hin, hin1, s=0;
int NTDTINV, nxny, nxnyi, imat, imat1, imat2, IDXI, IDYI, hi, NTST, NTSTI, partest, FREQFILT;
int lsnap, nsnap=0, lsamp=0, buffsize, invtime, invtimer, sws, swstestshot, snapseis, snapseis1, PML;
int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, nsrc_glob=0, ishot, irec, nshots=0, nshots1, Lcount, itest, Lcountsum, itestshot;

float pum, ppim, ppim1, ppim2, thetaf, thetab, e33, e33b, e11, e11b, muss, lamss; 
float memdyn, memmodel, memseismograms, membuffer, memtotal, dngn, fphi, sum, avggrad, beta, betan, betaz, betaLog, betaVp, betaVs, betarho, eps_scale, L2old;
float fac1, fac2, wavefor, waverecipro, dump, dump1, epsilon, gradsign, mun, eps1, gradplastiter, gradglastiter, gradclastiter, betar, sig_max, sig_max1;
float signL1, RMS, opteps_vp, opteps_vs, opteps_rho, Vs, Vp, Vp_avg, C_vp, Vs_avg, C_vs, Cd, rho_avg, C_rho, Vs_sum, Vp_sum, rho_sum, Zp, Zs;
float freqshift, dfreqshift, memfwt, memfwt1, memfwtdata;
char *buff_addr, ext[10], *fileinp;
char wave_forward[225], wave_recipro[225], wave_conv[225], jac[225], jac2[225], jacsum[225], dwavelet[225], vyf[STRING_SIZE];

double time1, time2, time3, time4, time5, time6, time7, time8,
	time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0, 
	time_av_s_exchange=0.0, time_av_timestep=0.0;
	
float L2, L2sum, L2_all_shots, L2sum_all_shots, *L2t, alphanomsum, alphanom, alphadenomsum, alphadenom, scaleamp ,sdummy, lamr; 

float energy, energy_sum, energy_all_shots, energy_sum_all_shots;	

float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
	** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvxdiff=NULL, ** sectionvxdiffold=NULL, ** sectionvydiffold=NULL,
	** sectionvydiff=NULL, ** sectionvydata=NULL, ** sectionpn=NULL, ** sectionread=NULL, ** sectionvy_conv=NULL, ** sectionvy_obs=NULL,** sectionvx_conv=NULL,** sectionvx_obs=NULL,
	* source_time_function=NULL;
float  **  absorb_coeff, ** taper_coeff, * epst1, * epst2,  * epst3, * picked_times;
float  ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, **srcpos_loc_back=NULL, ** signals=NULL, ** signals_rec=NULL, *hc=NULL, ** dsignals=NULL;
int   ** recpos=NULL, ** recpos_loc=NULL;
/*int   ** tracekill=NULL, TRKILL, DTRKILL;*/
int * DTINV_help;
     
float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top; 

/* PML variables */
float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half, * alpha_prime_x_half, * a_x_half, * b_x_half, * d_y, * K_y, * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half, * alpha_prime_y_half, * a_y_half, * b_y_half;
float ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y, ** psi_sxy_x, ** psi_vxx, ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;


/* Variables for viscoelastic modeling */
float **ptaus=NULL, **ptaup=NULL, *etaip=NULL, *etajm=NULL, *peta=NULL, **ptausipjp=NULL, **fipjp=NULL, ***dip=NULL, *bip=NULL, *bjm=NULL;
float *cip=NULL, *cjm=NULL, ***d=NULL, ***e=NULL, ***pr=NULL, ***pp=NULL, ***pq=NULL, **f=NULL, **g=NULL;

/* Variables for step length calculation */
int step1, step2, step3=0, itests, iteste, stepmax, countstep;
float scalefac;

/* Variables for Pseudo-Hessian calculation */
int RECINC, ntr1;
float * jac_rho, * jac_u, * jac_lam_x, * jac_lam_y;
float * temp_TS, * temp_TS1, * temp_TS2, * temp_TS3, * temp_TS4, * temp_TS5, * temp_conv, * temp_conv1, * temp_conv2;
float TSHIFT_back, temp_hess, mulamratio;
float ** hessian, ** hessian_u, ** hessian_rho;
int QUELLART_OLD;

/* Variables of the L-BFGS method */
float *** y_LBFGS_vp, *** s_LBFGS_vp, * rho_LBFGS, * alpha_LBFGS; 
float *** y_LBFGS_vs, *** s_LBFGS_vs;
float *** y_LBFGS_rho, *** s_LBFGS_rho;
int NLBFGS;
float * rho_LBFGS_vp, * rho_LBFGS_vs, * alpha_LBFGS_vp, * alpha_LBFGS_vs;

int * recswitch=NULL;
float ** fulldata=NULL, ** fulldata_vx=NULL, ** fulldata_vy=NULL;


/*vector for abort criterion*/
float * L2_hist=NULL;

/* help variable for MIN_ITER */
int min_iter_help=0;

/* variable for time domain filtering */
float FC;

FILE *fprec, *FP2, *FP3, *FP4, *FP5, *FPL2, *FP6, *FP7;
	
MPI_Request *req_send, *req_rec;
MPI_Status  *send_statuses, *rec_statuses;

/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

setvbuf(stdout, NULL, _IONBF, 0);

if (MYID == 0){
   time1=MPI_Wtime(); 
   clock();
  }

/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open Denise input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
		err(" --- ");
	}
}

if (strstr(fileinp,".json")){
	/* read json formatted input file */
	read_par_json(stdout,fileinp);}
else{
	/* read "old" input file *.inp */
	read_par(FP);}

NPROCX=1;
NPROCY=1;
	
exchange_par();

/* open log-file (each PE is using different file) */
/*	fp=stdout; */
sprintf(ext,".%i",MYID);  
strcat(LOG_FILE,ext);

if ((MYID==0) && (LOG==1)) FP=stdout;
else FP=fopen(LOG_FILE,"w");
fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
	
/* domain decomposition */
initproc();

NT=iround(TIME/DT);  	  /* number of timesteps */
/*ns=iround(NT/NDT);*/           /* number of samples per trace */
ns=NT;	/* in a FWI one has to keep all samples of the forward modeled data
	at the receiver positions to calculate the adjoint sources and to do 
	the backpropagation; look at function saveseis_glob.c to see that every
	NDT sample for the forward modeled wavefield is written to su files*/
lsnap=iround(TSNAP1/DT);      /* first snapshot at this timestep */
lsamp=NDT;


/* output of parameters to log-file or stdout */
if (MYID==0) write_par(FP);

	
/* NXG, NYG denote size of the entire (global) grid */
NXG=NX;
NYG=NY;

/* In the following, NX and NY denote size of the local grid ! */
NX = IENDX;
NY = IENDY;


if (SEISMO){
   recpos=receiver(FP, &ntr);
   recswitch = ivector(1,ntr);
   recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
   ntr_glob=ntr;
   ntr=ntr_loc;
}

fulldata = matrix(1,ntr_glob,1,NT);
fulldata_vx = matrix(1,ntr_glob,1,NT);
fulldata_vy = matrix(1,ntr_glob,1,NT);

/* memory allocation for abort criterion*/
L2_hist = vector(1,1000);

/* estimate memory requirement of the variables in megabytes*/
	
switch (SEISMO){
case 1 : /* particle velocities only */
	nseismograms=2;	
	break;	
case 2 : /* pressure only */
	nseismograms=1;	
	break;	
case 3 : /* curl and div only */
	nseismograms=2;		
	break;	
case 4 : /* everything */
	nseismograms=5;		
	break;
}	
	
/* use only every DTINV time sample for the inversion */
DTINV_help=ivector(1,NT);
NTDTINV=ceil((float)NT/(float)DTINV);		/* round towards next higher integer value */

/* save every IDXI and IDYI spatial point during the forward modelling */
IDXI=1;
IDYI=1;

/*allocate memory for dynamic, static and buffer arrays */
fac1=(NX+FDORDER)*(NY+FDORDER);
fac2=sizeof(float)*pow(2.0,-20.0);

nd = FDORDER/2 + 1;
fdo3 = 2*nd;

if (L){
	memdyn=(5.0+3.0*(float)L)*fac1*fac2;
	memmodel=(12.0+3.0*(float)L)*fac1*fac2;
	
} else {
	memdyn=5.0*fac1*fac2;
	memmodel=6.0*fac1*fac2;
}
memseismograms=nseismograms*ntr*ns*fac2;

memfwt=5.0*((NX/IDXI)+FDORDER)*((NY/IDYI)+FDORDER)*NTDTINV*fac2;
memfwt1=20.0*NX*NY*fac2;
memfwtdata=6.0*ntr*ns*fac2;

membuffer=2.0*fdo3*(NY+NX)*fac2;
buffsize=2.0*2.0*fdo3*(NX +NY)*sizeof(MPI_FLOAT);
memtotal=memdyn+memmodel+memseismograms+memfwt+memfwt1+memfwtdata+membuffer+(buffsize*pow(2.0,-20.0));


if (MYID==0){
   fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
   fprintf(FP," Size of local grids: NX=%d \t NY=%d\n",NX,NY);
   fprintf(FP," Each process is now trying to allocate memory for:\n");
   fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
   fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
   fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
   fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
   fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
   fprintf(FP," ------------------------------------------------ \n");
   fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
   }


/* allocate buffer for buffering messages */
buff_addr=malloc(buffsize);
if (!buff_addr) err("allocation failure for buffer for MPI_Bsend !");
MPI_Buffer_attach(buff_addr,buffsize);

NTST=20;
NTSTI=NTST/DTINV;

nxny=NX*NY;
nxnyi=(NX/IDXI)*(NY/IDYI);

if (ntr>0){
	switch (SEISMO){
	case 1 : /* particle velocities only */
		sectionvx=matrix(1,ntr,1,ns);
		sectionvy=matrix(1,ntr,1,ns);	
		break;	
	case 2 : /* pressure only */
		sectionp=matrix(1,ntr,1,ns);
		sectionpnp1=matrix(1,ntr,1,ns);
		sectionpn=matrix(1,ntr,1,ns);
		break;	
	case 3 : /* curl and div only */
		sectioncurl=matrix(1,ntr,1,ns);
		sectiondiv=matrix(1,ntr,1,ns);	
		break;	
	case 4 : /* everything */
		sectionvx=matrix(1,ntr,1,ns);
		sectionvy=matrix(1,ntr,1,ns);	
		sectioncurl=matrix(1,ntr,1,ns);
		sectiondiv=matrix(1,ntr,1,ns);		
		sectionp=matrix(1,ntr,1,ns);
		break;	
	
	}
}	

/* Memory for seismic data */
sectionvxdata=matrix(1,ntr,1,ns);
sectionread=matrix(1,ntr_glob,1,ns);
sectionvxdiff=matrix(1,ntr,1,ns);
sectionvydata=matrix(1,ntr,1,ns);
sectionvydiff=matrix(1,ntr,1,ns);
sectionvxdiffold=matrix(1,ntr,1,ns);
sectionvydiffold=matrix(1,ntr,1,ns);

picked_times = vector(1,ntr);
	
fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

		
/* Reading source positions from SOURCE_FILE */ 	
srcpos=sources(&nsrc);
nsrc_glob=nsrc;

nshots = nsrc_glob;
SHOTINC=1;

for (ishot=1;ishot<=nshots;ishot+=SHOTINC){

if (MYID==0){
printf("Picking first arrivals for shot %i of %i  \n",ishot,nshots);
}

TIME_FILT=1;
FC=10.0;
ORDER=5;

/* read seismic data from SU file vx */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

	if(INV_STF!=0){
			inseis(fprec,ishot,sectionread,ntr_glob,ns,7,iter);
			}
	else{
			inseis(fprec,ishot,sectionread,ntr_glob,ns,1,iter);
			}

if ((TIME_FILT==1)&&(INV_STF==0)){
  timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);
}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

  stalta(sectionvxdata, ntr, ns, picked_times, ishot);
  /*mer(sectionvxdata, ntr, ns, picked_times, ishot);*/

}

/* read seismic data from SU file vy */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */

	if(INV_STF!=0){
			inseis(fprec,ishot,sectionread,ntr_glob,ns,8,iter);
			}
	else{
			inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);
	}

if ((TIME_FILT==1)&&(INV_STF==0)){
  timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);
}
  
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

  stalta(sectionvydata, ntr, ns, picked_times, ishot);
  /*mer(sectionvydata, ntr, ns, picked_times, ishot);*/  

} /* end QUELLTYPB */

} /* end of loop over shots */

MPI_Finalize();
return 0;	

}/*main*/
