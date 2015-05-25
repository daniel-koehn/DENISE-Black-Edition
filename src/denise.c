/*------------------------------------------------------------------------
 *  DENISE Black Edition: 2D isotropic elastic time domain FWI Code 
 *
 *
 *  Authors:
 *  -----------  
 * 
 *  D. Koehn    (FWI code + updates)
 *  D. De Nil   (FWI code + updates)
 *  L. Rehor    (viscoelastic modelling, Butterworth-filter)
 *  A. Kurzmann (original step length estimation)
 *  M. Schaefer (source wavelet inversion)
 *  S. Heider   (time-windowing)
 *  T. Bohlen   (original FD forward code) 
 *  L. Zhang    (towed streamer, pressure inversion)
 *  
 *  
 *  In case of questions contact the author:
 *	Dr. Daniel Koehn, Kiel University, Institute of Geoscience,
 *	Otto-Hahn-Platz 1, D-24098 Kiel, Germany, ph: +49 431 880 4566,
 *	mailto:dkoehn@geophysik.uni-kiel.de,
 *	Homepage: http://www.geophysik.uni-kiel.de/~dkoehn
 *
 *
 *  DENISE Black Edition is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, version 2.0 of the License only. 
 *  
 *  DENISE Black Edition is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with DENISE Black Edition (see file LICENSE.md) 
 *
 *  If you show modelling/inversion results in a paper or presentation please 
 *  give a reference to the following papers:
 *
 *  Daniel Koehn, Denise De Nil, Andre Kurzmann, Anna Przebindowska and Thomas Bohlen (2012): 
 *  On the influence of model parametrization in elastic full waveform tomography, 
 *  Geophysical Journal International, 191(1), 325-345.
 *
 *  Daniel Koehn (2011): Time Domain 2D Elastic Full Waveform Tomography, PhD-Thesis, Kiel University
 *  Available at: http://nbn-resolving.de/urn:nbn:de:gbv:8-diss-67866 
 * 
 *  
 *  Thank you for your co-operation, 
 *  Daniel Koehn
 * 
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */
#include "cseife.h"

int main(int argc, char **argv){
/* variables in main */
int ns, nseismograms=0, nt, nd, fdo3, j, i, ii, jj, shotid, recid, k, nc, iter, h, iter_true, infoout, SHOTINC, TIMEWIN, test_eps, lq, iq, jq, hin, hin1, s=0;
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

float  **  psxx, **  psxy, **  psyy, ** ux, ** uy, ** uxy, ** uyx, ** Vp0, ** uttx, ** utty, ** Vs0, ** Rho0;
float  **  pvx, **  pvy, **waveconv, **waveconv_lam, **waveconv_mu, **waveconv_rho, **waveconv_rho_s, **waveconv_u, **waveconvtmp, **wcpart, **wavejac;
float **waveconv_shot, **waveconv_u_shot, **waveconv_rho_shot;
float  **  pvxp1, **  pvyp1, **  pvxm1, **  pvym1;
float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u;
float  **  prho,**  prhonp1, **prip=NULL, **prjp=NULL, **pripnp1=NULL, **prjpnp1=NULL, **  ppi, **  pu, **  punp1, **  puipjp, **  ppinp1;
float  **  vpmat, *forward_prop_x, *forward_prop_y, *forward_prop_rho_x, *forward_prop_u, *forward_prop_rho_y;

float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
	** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvxdiff=NULL, ** sectionvxdiffold=NULL, ** sectionvydiffold=NULL,
	** sectionvydiff=NULL, ** sectionvydata=NULL, ** sectionpn=NULL, ** sectionread=NULL, ** sectionvy_conv=NULL, ** sectionvy_obs=NULL,** sectionvx_conv=NULL,** sectionvx_obs=NULL,
	* source_time_function=NULL;
float ** sectionpdata=NULL, ** sectionpdiff=NULL, ** sectionpdiffold=NULL;	
float  **  absorb_coeff, ** taper_coeff, * epst1, * epst2,  * epst3;
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
int step1, step2, step3=0, itests, iteste, countstep;
float eps_true, tmp, tmp1;

/* Variables for Pseudo-Hessian calculation */
int RECINC, ntr1;
float ** hessian_lam=NULL, ** hessian_mu=NULL, ** hessian_rho=NULL;

/* Variables for Laplace-domain inversion */
float  ** forward_propl_x, ** forward_propl_y, ** forward_propl_rho_x, ** forward_propl_u, ** forward_propl_rho_y;
float  ** back_prop_x, ** back_prop_y, ** back_prop_rho_x, ** back_prop_u, ** back_prop_rho_y;
float time;

/* Variables for the L-BFGS method */
float * rho_LBFGS, * alpha_LBFGS, * beta_LBFGS; 
float * y_LBFGS, * s_LBFGS, * q_LBFGS, * r_LBFGS;
int NLBFGS_class, LBFGS_pointer, NLBFGS_vec;

/* Variables for energy weighted gradient */
float ** Ws, **Wr, **We;

/* parameters for FWI-workflow */
int stagemax=0, nstage;

int * recswitch=NULL;
float ** fulldata=NULL, ** fulldata_vx=NULL, ** fulldata_vy=NULL;
float ** fulldata_p=NULL, ** fulldata_curl=NULL,  ** fulldata_div=NULL;

/*vector for abort criterion*/
float * L2_hist=NULL;

/* help variable for MIN_ITER */
int min_iter_help=0;

/* parameters for gravity inversion */
float * gz_mod, * gz_res;
float ** gravpos=NULL, ** rho_grav=NULL, ** rho_grav_ext=NULL;
float ** grad_grav=NULL;
int ngrav=0, nxgrav, nygrav;
float L2_grav, FWImax, GRAVmax, FWImax_all, GRAVmax_all ;
char jac_grav[STRING_SIZE];

/* variable for time domain filtering */
/*float FC;*/

FILE *fprec, *FP2, *FP3, *FP4, *FP5, *FPL2, *FP6, *FP7, *FP_stage, *FP_GRAV;
	
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
		

/* print program name, version etc to stdout*/
if (MYID == 0) info(stdout);

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


/* read input file *.inp */
read_par(FP);

/* define default values */
SPATFILTER=0;
WD_DAMP=1e-4;
SPAT_FILT_SIZE=40;
SPAT_FILT_1=1;
SPAT_FILT_ITER=0;

INV_RHO_ITER=0;

HESSIAN=0;
FC_HESS_START=4.0;
FC_HESS_INC=0.4;
NFREQ=40;

INV_STF=0;
N_STF=10;
N_STF_START=1;

TIME_FILT=1;
FC_START=5.0;
FC_END=80.0;
FC_INCR=2.0;
ORDER=8;

LNORM=2;

PRO=0.0;

TIMEWIN=0;
TWLENGTH_PLUS=0.01;
TWLENGTH_MINUS=0.01;
GAMMA=1e11;

NORMALIZE=0;

INV_VP_ITER=0;
INV_VS_ITER=0;
	
exchange_par();

/* read parameters from workflow-file (stdin) */
fileinp=argv[2];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open Denise workflow input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
		err(" --- ");
	}
}

/* estimate number of lines in FWI-workflow */
i=0;
stagemax=0;
while ((i=fgetc(FP)) != EOF)
if (i=='\n') ++stagemax;
rewind(FP);
stagemax--;
fclose(FP);

if (MYID == 0) note(stdout);

 
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

   recpos=receiver(FP, &ntr, ishot);
   recswitch = ivector(1,ntr);
   recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
   ntr_glob=ntr;
   ntr=ntr_loc;
   
   if(N_STREAMER>0){
     free_imatrix(recpos,1,3,1,ntr_glob);
     if(ntr>0) free_imatrix(recpos_loc,1,3,1,ntr);
     free_ivector(recswitch,1,ntr_glob);
   }
   
}

if(N_STREAMER==0){

   if (ntr>0){
           switch (SEISMO){
           case 1 : /* particle velocities only */
                   sectionvx=matrix(1,ntr,1,ns);
                   sectionvy=matrix(1,ntr,1,ns);
                   break;
           case 2 : /* pressure only */
                   sectionp=matrix(1,ntr,1,ns);
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
           case 5: /* pressure and particle velocities only */
                   sectionvx=matrix(1,ntr,1,ns);
                   sectionvy=matrix(1,ntr,1,ns);
                   sectionp=matrix(1,ntr,1,ns);
                   break;
           }
   }

   /* Memory for seismic data */
   sectionread=matrix(1,ntr_glob,1,ns);
   
   if((QUELLTYPB==1)||(QUELLTYPB==3)){
      sectionvxdata=matrix(1,ntr,1,ns);
      sectionvxdiff=matrix(1,ntr,1,ns);
      sectionvxdiffold=matrix(1,ntr,1,ns);
   }
   
   if((QUELLTYPB==1)||(QUELLTYPB==2)){
      sectionvydata=matrix(1,ntr,1,ns);
      sectionvydiff=matrix(1,ntr,1,ns);
      sectionvydiffold=matrix(1,ntr,1,ns);
   }
   
   if(QUELLTYPB==4){
     sectionpdata=matrix(1,ntr,1,ns);
     sectionpdiff=matrix(1,ntr,1,ns);
     sectionpdiffold=matrix(1,ntr,1,ns);
   }

}

if(SEISMO){
  fulldata = matrix(1,ntr_glob,1,NT);
}

if(SEISMO==1){
  fulldata_vx = matrix(1,ntr_glob,1,NT);
  fulldata_vy = matrix(1,ntr_glob,1,NT);
}

if(SEISMO==2){
  fulldata_p = matrix(1,ntr_glob,1,NT);
}

if(SEISMO==3){
  fulldata_curl = matrix(1,ntr_glob,1,NT); 
  fulldata_div = matrix(1,ntr_glob,1,NT);
}

if(SEISMO==4){
  fulldata_vx = matrix(1,ntr_glob,1,NT);
  fulldata_vy = matrix(1,ntr_glob,1,NT);
  fulldata_p = matrix(1,ntr_glob,1,NT); 
  fulldata_curl = matrix(1,ntr_glob,1,NT);
  fulldata_div = matrix(1,ntr_glob,1,NT); 
}

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

/* allocation for request and status arrays */
req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));

/* memory allocation for dynamic (wavefield) arrays */
psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uy   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uxy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uyx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconvtmp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
wcpart = matrix(1,3,1,3);
wavejac = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
/* memory allocation for static (model) arrays */
prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prhonp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pripnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prjpnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ppinp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
punp1   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
vpmat   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

/* dynamic (wavefield) arrays for viscoelastic modeling */
if (L > 0) {
	pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
}

/* memory allocation for static arrays for viscoelastic modeling */
if (L>0){
	dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	ptausipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	peta =  vector(1,L);
	etaip =  vector(1,L);
	etajm =  vector(1,L);
	bip =  vector(1,L);
	bjm =  vector(1,L);
	cip =  vector(1,L);
	cjm =  vector(1,L);
}

/* Activate energy weighted preconditioning */
EPRECOND=1;

/* FREQFILT = 1 activate frequency filter */
/* FREQFILT = 0 deactivate frequency filter */
FREQFILT=0;

/* If INVMAT==10 deactivate the inversion code and set ITERMAX=1 */
if(INVMAT==10){
 ITERMAX=1;
 TIMELAPSE=0;
 stagemax=1;
 RTM=0;
}

/* Define gradient formulation */
/* GRAD_FORM = 1 - stress-displacement */
/* GRAD_FORM = 2 - stress-velocity in non-conservative form (Castellanos, 2014) */
/* GRAD_FORM = 3 - stress-velocity */
GRAD_FORM = 1;

/* Parameters for gravity modeling/inversion */
/* GRAVITY == 0 no gravity modelling or inversion */
/* GRAVITY == 1 activate gravity modeling */
/* GRAVITY == 2 activate joint FW-gravity inversion */
GRAVITY = 0; 

if(GRAVITY){
  
  sprintf(GRAV_DATA_OUT, "./gravity/grav_mod.dat"); /* output file of gravity data */
  sprintf(GRAV_DATA_IN, "./gravity/grav_field.dat");  /* input file of gravity data */
  sprintf(GRAV_STAT_POS, "./gravity/grav_stat.dat"); /* file with station positions for gravity modelling */   
  NZGRAV = 250;  /* number of grid points in z-direction */
  NGRAVB = 500;  /* thickness of the boundaries [grid points] */

  /* size of the extended gravity model */
  nxgrav = NXG + 2*NGRAVB;
  nygrav = NYG + NGRAVB;

  /* weighting parameter gamma for the gravity objective function */
  GAMMA_GRAV = 5e-1;
  /*GAMMA_GRAV = 0.0;*/

}

if(RTMOD==1){
  RTM=1;
}

/* Activate DC removal from the data residuals */
DC_REMOVE = 0;
FC_DC = 4.0;
ORDER_DC = 6;

if(RTM==1){
  ITERMAX=1;
  TIMELAPSE=0;
  INVMAT=0; 
}
      
lamr=0.01; /* Marquardt factor */

/*nf=4;
nfstart=4;*/

NTST=20;
NTSTI=NTST/DTINV;

nxny=NX*NY;
nxnyi=(NX/IDXI)*(NY/IDYI);

/* Variables for the L-BFGS method */

if(GRAD_METHOD==2){

  NLBFGS_class = 3;                 /* number of parameter classes */ 
  NLBFGS_vec = NLBFGS_class*NX*NY;  /* length of one LBFGS-parameter class */
  LBFGS_pointer = 1;                /* initiate pointer in the cyclic LBFGS-vectors */
  
  y_LBFGS  =  vector(1,NLBFGS_vec*NLBFGS);
  s_LBFGS  =  vector(1,NLBFGS_vec*NLBFGS);

  q_LBFGS  =  vector(1,NLBFGS_vec);
  r_LBFGS  =  vector(1,NLBFGS_vec);

  rho_LBFGS = vector(1,NLBFGS);
  alpha_LBFGS = vector(1,NLBFGS);
  beta_LBFGS = vector(1,NLBFGS);
  
}

if(INVMAT<=1){

  gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  
  gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

  gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

}

if(INVMAT==0){

  forward_prop_x =  vector(1,nxnyi*(NTDTINV));
  forward_prop_y =  vector(1,nxnyi*(NTDTINV));
  forward_prop_rho_x =  vector(1,nxnyi*(NTDTINV));
  forward_prop_rho_y =  vector(1,nxnyi*(NTDTINV));
  forward_prop_u =  vector(1,nxnyi*(NTDTINV));

}

if(INVMAT==1){

  forward_propl_x = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  forward_propl_y = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  forward_propl_rho_x = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  forward_propl_rho_y = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  forward_propl_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

  back_prop_x = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  back_prop_y = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  back_prop_rho_x = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  back_prop_rho_y = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
  back_prop_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

}

if((EPRECOND==1)||(EPRECOND==3)){
  Ws = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the source wavefield */
  Wr = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the receiver wavefield */
  We = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of source and receiver wavefield */
}

if(FW>0){

  d_x = vector(1,2*FW);
  K_x = vector(1,2*FW);
  alpha_prime_x = vector(1,2*FW);
  a_x = vector(1,2*FW);
  b_x = vector(1,2*FW);
  
  d_x_half = vector(1,2*FW);
  K_x_half = vector(1,2*FW);
  alpha_prime_x_half = vector(1,2*FW);
  a_x_half = vector(1,2*FW);
  b_x_half = vector(1,2*FW);

  d_y = vector(1,2*FW);
  K_y = vector(1,2*FW);
  alpha_prime_y = vector(1,2*FW);
  a_y = vector(1,2*FW);
  b_y = vector(1,2*FW);
  
  d_y_half = vector(1,2*FW);
  K_y_half = vector(1,2*FW);
  alpha_prime_y_half = vector(1,2*FW);
  a_y_half = vector(1,2*FW);
  b_y_half = vector(1,2*FW);

  psi_sxx_x =  matrix(1,NY,1,2*FW); 
  psi_syy_y =  matrix(1,2*FW,1,NX);
  psi_sxy_y =  matrix(1,2*FW,1,NX);
  psi_sxy_x =  matrix(1,NY,1,2*FW);
  psi_vxx   =  matrix(1,NY,1,2*FW);
  psi_vxxs  =  matrix(1,NY,1,2*FW); 
  psi_vyy   =  matrix(1,2*FW,1,NX);
  psi_vxy   =  matrix(1,2*FW,1,NX);
  psi_vyx   =  matrix(1,NY,1,2*FW);
   
}

taper_coeff=  matrix(1,NY,1,NX);



/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */
bufferlef_to_rig = matrix(1,NY,1,fdo3);
bufferrig_to_lef = matrix(1,NY,1,fdo3);
buffertop_to_bot = matrix(1,NX,1,fdo3);
bufferbot_to_top = matrix(1,NX,1,fdo3);

/* memory for source position definition */
srcpos1=fmatrix(1,8,1,1);

/* memory of L2 norm */
L2t = vector(1,4);
epst1 = vector(1,3);
epst2 = vector(1,3);
epst3 = vector(1,3);
	
fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

		
/* Holberg coefficients for FD operators*/
hc = holbergcoeff();

MPI_Barrier(MPI_COMM_WORLD);

/* Reading source positions from SOURCE_FILE */ 	
srcpos=sources(&nsrc);
nsrc_glob=nsrc;


/* create model grids */

if(L){
	if (READMOD) readmod(prho,ppi,pu,ptaus,ptaup,peta);
		else model(prho,ppi,pu,ptaus,ptaup,peta);
} else{
	if (READMOD) readmod_elastic(prho,ppi,pu);
    		else model_elastic(prho,ppi,pu);
}

/* check if the FD run will be stable and free of numerical dispersion */
if(L){
	checkfd_ssg_visc(FP,prho,ppi,pu,ptaus,ptaup,peta,hc);
} else{
	checkfd_ssg_elastic(FP,prho,ppi,pu,hc);
}


/* calculate damping coefficients for CPMLs*/
if(FW>0){PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half, 
               d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half);}

MPI_Barrier(MPI_COMM_WORLD);

if (CHECKPTREAD){
	if (MYID==0){
		time3=MPI_Wtime();
 		fprintf(FP," Reading wavefield from check-point file %s \n",CHECKPTFILE);	
	}
	
	read_checkpoint(-1, NX+2, -1, NY+2, pvx, pvy, psxx, psyy, psxy);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0){
		time4=MPI_Wtime();
      		fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
	}
}

if(GRAVITY){
 
  /* read station positions */
  MPI_Barrier(MPI_COMM_WORLD);
  gravpos=read_grav_pos(&ngrav);

  /* define model and residual data vector for gz (z-component of the gravity field) */
  gz_mod = vector(1,ngrav);
  gz_res = vector(1,ngrav);

  /* only model gravity data */
  if(GRAVITY==1){

    /* global density model */
    rho_grav =  matrix(1,NYG,1,NXG);
    rho_grav_ext =  matrix(1,nygrav,1,nxgrav);

    read_density_glob(rho_grav,1);
    extend_mod(rho_grav,rho_grav_ext,nxgrav,nygrav);
    grav_mod(rho_grav_ext,ngrav,gravpos,gz_mod,nxgrav,nygrav,NZGRAV);

    free_matrix(rho_grav,1,NYG,1,NXG);
    free_matrix(rho_grav_ext,1,nygrav,1,nxgrav);

  }

  if(GRAVITY==2){
    grad_grav =  matrix(1,NY,1,NX);
  }

}  
      
snapseis=1;
snapseis1=1;
SHOTINC=1;
RECINC=1;
    
iter_true=1;
/* Begin of FWI-workflow */
for(nstage=1;nstage<=stagemax;nstage++){

/* read workflow input file *.inp */
fileinp=argv[2];
FP_stage=fopen(fileinp,"r");
read_par_inv(FP_stage,nstage,stagemax);

FC=FC_END;

if(INVMAT==10){TIME_FILT=0;}

iter=1;
/* --------------------------------------
 * Begin of Full Waveform iteration loop
 * -------------------------------------- */
while(iter<=ITERMAX){

if(EPRECOND==2){
  hessian_lam  = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* Pseudo-Hessian for Lambda */
  hessian_mu  = matrix(-nd+1,NY+nd,-nd+1,NX+nd);  /* Pseudo-Hessian for Mu */
  hessian_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);  /* Pseudo-Hessian for Rho */
}

if(GRAD_METHOD==2){
  
  /* increase pointer to LBFGS-vector*/
  if(iter>2){
    LBFGS_pointer++;
  }
  
  /* if LBFGS-pointer > NLBFGS -> set LBFGS_pointer=1 */ 
  if(LBFGS_pointer>NLBFGS){LBFGS_pointer=1;}

}


if (MYID==0)
   {
   time2=MPI_Wtime();
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   fprintf(FP,"\n\n\n                   TDFWI ITERATION %d \t of %d \n",iter,ITERMAX);
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   }

/* For the calculation of the material parameters between gridpoints
   they have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
if (L){
	matcopy(prho,ppi,pu,ptaus,ptaup);
} else{
	matcopy_elastic(prho, ppi, pu);
}

MPI_Barrier(MPI_COMM_WORLD);

av_mue(pu,puipjp,prho);
av_rho(prho,prip,prjp);
if (L) av_tau(ptaus,ptausipjp);

/*printf("pu = %f \n",pu[10][10]);
printf("puipjp = %f \n",puipjp[10][10]);
printf("INVMAT1 = %d \n",INVMAT1);*/

/* Preparing memory variables for update_s (viscoelastic) */
if (L) prepare_update_s(etajm,etaip,peta,fipjp,pu,puipjp,ppi,prho,ptaus,ptaup,ptausipjp,f,g,
		bip,bjm,cip,cjm,dip,d,e);


if(iter_true==1){
    for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	
	if(INVMAT1==1){
	
	  Vp0[j][i] = ppi[j][i];
	  Vs0[j][i] = pu[j][i];
	  Rho0[j][i] = prho[j][i];
        }
	  
                 
		 
	if(INVMAT1==2){
        
	  Vp0[j][i] = sqrt((ppi[j][i]+2.0*pu[j][i])*prho[j][i]);
	  Vs0[j][i] = sqrt((pu[j][i])*prho[j][i]);
	  Rho0[j][i] = prho[j][i];
	
	}
	 
	if(INVMAT1==3){
        
	  Vp0[j][i] = ppi[j][i];
	  Vs0[j][i] = pu[j][i];
	  Rho0[j][i] = prho[j][i];
	
	}  
	
    }
    }

/* ----------------------------- */
/* calculate Covariance matrices */
/* ----------------------------- */

	 Vp_avg = 0.0;
	 Vs_avg = 0.0;
	 rho_avg = 0.0;
	 
        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	  
		 /* calculate average Vp, Vs */
                 Vp_avg+=ppi[j][i];
		 Vs_avg+=pu[j][i];
		 
		 /* calculate average rho */
		 rho_avg+=prho[j][i];
	
           }
        }
		
        /* calculate average Vp, Vs and rho of all CPUs*/
        Vp_sum = 0.0;
        MPI_Allreduce(&Vp_avg,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        Vp_avg=Vp_sum;
	
	Vs_sum = 0.0;
        MPI_Allreduce(&Vs_avg,&Vs_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        Vs_avg=Vs_sum;
	
	rho_sum = 0.0;
        MPI_Allreduce(&rho_avg,&rho_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        rho_avg=rho_sum;
	
	Vp_avg /=NXG*NYG; 
	Vs_avg /=NXG*NYG; 
	rho_avg /=NXG*NYG;
	
	if(MYID==0){
           printf("Vp_avg = %e \t Vs_avg = %e \t rho_avg = %e \n ",Vp_avg,Vs_avg,rho_avg);	
	}
	
	C_vp = Vp_avg;
	C_vs = Vs_avg;
	C_rho = rho_avg;


}

/* Open Log File for L2 norm */

if(INVMAT!=10){
  
  if(MYID==0){
    if(iter_true==1){
      FPL2=fopen(MISFIT_LOG_FILE,"w");
    }

    if(iter_true>1){
      FPL2=fopen(MISFIT_LOG_FILE,"a");
    }
  }
}

/* initialization of L2 calculation */
L2=0.0;
Lcount=0;
energy=0.0;
L2_all_shots=0.0;
energy_all_shots=0.0;

EPSILON=0.0;  /* test step length */
exchange_par();

/* set gradient and preconditioning matrices 0 before next iteration*/
for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	
	   if(INVMAT<=1){
	     waveconv[j][i]=0.0;
	     waveconv_rho[j][i]=0.0;
	     waveconv_u[j][i]=0.0;    
           }
          
        }   
}


itestshot=TESTSHOT_START;
swstestshot=0;

if(INVTYPE==2){ 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;

for (ishot=1;ishot<=nshots;ishot+=SHOTINC){
/* for (ishot=1;ishot<=1;ishot+=1){*/

if(N_STREAMER>0){

   if (SEISMO){
      recpos=receiver(FP, &ntr, ishot);
      recswitch = ivector(1,ntr);
      recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
      ntr_glob=ntr;
      ntr=ntr_loc;
   }

   if (ntr>0){
           switch (SEISMO){
           case 1 : /* particle velocities only */
                   sectionvx=matrix(1,ntr,1,ns);
                   sectionvy=matrix(1,ntr,1,ns);
                   break;
           case 2 : /* pressure only */
                   sectionp=matrix(1,ntr,1,ns);
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
           case 5: /* pressure and particle velocities only */
                   sectionvx=matrix(1,ntr,1,ns);
                   sectionvy=matrix(1,ntr,1,ns);
                   sectionp=matrix(1,ntr,1,ns);
                   break;
           }
   }

   /* Memory for seismic data */
   sectionread=matrix(1,ntr_glob,1,ns);
   
   if((QUELLTYPB==1)||(QUELLTYPB==3)){
      sectionvxdata=matrix(1,ntr,1,ns);
      sectionvxdiff=matrix(1,ntr,1,ns);
      sectionvxdiffold=matrix(1,ntr,1,ns);
   }
   
   if((QUELLTYPB==1)||(QUELLTYPB==2)){
      sectionvydata=matrix(1,ntr,1,ns);
      sectionvydiff=matrix(1,ntr,1,ns);
      sectionvydiffold=matrix(1,ntr,1,ns);
   }
   
   if(QUELLTYPB==4){
     sectionpdata=matrix(1,ntr,1,ns);
     sectionpdiff=matrix(1,ntr,1,ns);
     sectionpdiffold=matrix(1,ntr,1,ns);
   }

}

/* estimate source time function by Wiener deconvolution */
if((INV_STF)&&(iter==1)&&(INVMAT<=1)){

        fprintf(FP,"\n==================================================================================\n");
        fprintf(FP,"\n MYID=%d *****  Starting simulation (STF) for shot %d of %d  ********** \n",MYID,ishot,nshots);
	fprintf(FP,"\n==================================================================================\n\n");

                for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
                if (RUN_MULTIPLE_SHOTS){

                                /* find this single source positions on subdomains */
                                if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
                                srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		        }


                else{
                                /* Distribute multiple source positions on subdomains */
                               srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                }

QUELLART = 6;
MPI_Barrier(MPI_COMM_WORLD);
/* calculate wavelet for each source point */
signals=NULL;
signals=wavelet(srcpos_loc,nsrc_loc,ishot);

if (nsrc_loc){if(QUELLART==6){

   /* low pass filtering of spike */
   timedomain_filt(signals,FC_SPIKE_2,ORDER_SPIKE,nsrc_loc,ns,1);

   /* band-pass filtering of spike */
   if(FC_SPIKE_1 > 0.0){
   timedomain_filt(signals,FC_SPIKE_1,ORDER_SPIKE,nsrc_loc,ns,2);
   }

   }
}

/* time domain filtering*/
/*if ((TIME_FILT)&&(INVMAT!=10)){*/

   /*time domain low pass filtering of the source signal */
   /*timedomain_filt(signals,FC,ORDER,nsrc_loc,ns,1);*/

   /*time domain band-pass filtering of the source signal */
   /*if(TIME_FILT==2){
     timedomain_filt(signals,FC_START,ORDER,nsrc_loc,ns,2);
   }*/

/*}*/

/*char  source_signal_file[STRING_SIZE];
sprintf(source_signal_file,"source_signal.%d.su.shot%d.it%d",MYID,ishot,iter);
fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
output_source_signal(fopen(source_signal_file,"w"),signals,NT,3);*/

/* output source signal e.g. for cross-correlation of comparison with analytical solutions */
/*if(RUN_MULTIPLE_SHOTS){

	if(nsrc_loc>0){
        	   char  source_signal_file[STRING_SIZE];
        	   sprintf(source_signal_file,"%s_source_signal.%d.su.shot%d", MFILE, MYID,ishot);
        	   fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
        	   output_source_signal(fopen(source_signal_file,"w"),signals,NT,1);
	}                                
                                
	MPI_Barrier(MPI_COMM_WORLD);
}*/
		    
/* initialize wavefield with zero */
if (L){
	zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pr,pp,pq);
}else{	
	zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	
}                                                         
     
/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%10000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /* update of particle velocities */
      update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                         
	if (MYID==0){
		time4=MPI_Wtime();
		time_av_v_update+=(time4-time3);
		if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
	}
                                                   
	/* exchange of particle velocities between PEs */
	exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);
                                                               
	if (MYID==0){
	  time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
	  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
	}                                                                                         	

    if (L)    /* viscoelastic */
    	update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout,
				pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
				K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
    else
   	update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
         	               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  


    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc,0);
   
   /* moment tensor source */
   if ((!CHECKPTREAD)&&(QUELLTYP==5)) 	
   msource(nt,psxx,psyy,psxy,srcpos_loc,signals,nsrc_loc,0);

   if ((FREE_SURF) && (POS[2]==0)){
   	if (L)    /* viscoelastic */
		surface_PML(1, pvx, pvy, psxx, psyy, psxy, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs);
	else      /* elastic */
   		surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);
   }


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	  time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  

	/* store amplitudes at receivers in section-arrays */
	if (SEISMO){
		seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		/*lsamp+=NDT;*/
	}

      
   if (MYID==0){
      time8=MPI_Wtime();
	  time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (STF estimation) ----------*/
		
   catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);

   /*MPI_Barrier(MPI_COMM_WORLD);   */

   /* estimate STF */	
   if (nsrc_loc>0){
       
      /* read seismic data from SU file vy */
      /* --------------------------------- */
      inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);

      if (TIME_FILT){
         timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

         if (TIME_FILT==2){  /* apply band-pass*/
            timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
         }

      }

      stf(sectionread,fulldata_vy,ntr_glob,ishot,ns,iter,nshots,signals,recpos,srcpos);
      /*saveseis_glob(FP,sectionread,fulldata_vy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);*/
      
   }

   MPI_Barrier(MPI_COMM_WORLD);

   QUELLART=3;

}


 
        fprintf(FP,"\n==================================================================================\n");
        fprintf(FP,"\n MYID=%d *****  Starting simulation (forward model) for shot %d of %d  ********** \n",MYID,ishot,nshots);
	fprintf(FP,"\n==================================================================================\n\n");
		
                for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
                if (RUN_MULTIPLE_SHOTS){

                                /* find this single source positions on subdomains */
                                if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
                                srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		        }


                else{
                                /* Distribute multiple source positions on subdomains */
                               srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                }

MPI_Barrier(MPI_COMM_WORLD);
/* calculate wavelet for each source point */
signals=NULL;
signals=wavelet(srcpos_loc,nsrc_loc,ishot);

if (nsrc_loc){if(QUELLART==6){

   /* low pass filtering of spike */
   timedomain_filt(signals,FC_SPIKE_2,ORDER_SPIKE,nsrc_loc,ns,1);

   /* band-pass filtering of spike */
   if(FC_SPIKE_1 > 0.0){
   timedomain_filt(signals,FC_SPIKE_1,ORDER_SPIKE,nsrc_loc,ns,2);
   }
   
   }
}

/* time domain filtering*/
if ((TIME_FILT)&&(INVMAT!=10)&&(INV_STF==0)){
	
	/* time domain filtering of the source signal */
	timedomain_filt(signals,FC,ORDER,nsrc_loc,ns,1);

        if(TIME_FILT==2){ /* band-pass filter */
          timedomain_filt(signals,FC_START,ORDER,nsrc_loc,ns,2);
        }

}

/*printf("MYID=%d, nsrc_loc = %d \n",MYID,nsrc_loc);*/

/*char  source_signal_file[STRING_SIZE];
sprintf(source_signal_file,"source_signal.%d.su.shot%d.it%d",MYID,ishot,iter);
fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
output_source_signal(fopen(source_signal_file,"w"),signals,NT,3);*/

/* output source signal e.g. for cross-correlation of comparison with analytical solutions */
if(RUN_MULTIPLE_SHOTS){

	if(nsrc_loc>0){
        	   char  source_signal_file[STRING_SIZE];
        	   sprintf(source_signal_file,"%s_source_signal.%d.su.shot%d", MFILE, MYID,ishot);
        	   fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
        	   output_source_signal(fopen(source_signal_file,"w"),signals,NT,1);
	}                                
                                
	MPI_Barrier(MPI_COMM_WORLD);
}
		    
/* initialize wavefield with zero */
if (L){
	zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pr,pp,pq);
}else{	
	zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	
}


/*initialize gradient matrices for each shot with zeros*/
if(INVMAT<=1){

	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_shot[j][i]=0.0;
		}
	}

	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_rho_shot[j][i]=0.0;
		}
	}

	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_u_shot[j][i]=0.0;

                   if((EPRECOND==1)||(EPRECOND==3)){
                     Ws[j][i]=0.0;
                     Wr[j][i]=0.0;
	             We[j][i]=0.0;
                   }
		}
	}
}

/*initialize Laplace-domain wavefields for each shot with zeros*/
if(INVMAT==1){

   for(i=1;i<=NX;i=i+IDX){
      for(j=1;j<=NY;j=j+IDY){
      
         forward_propl_rho_x[j][i]=0.0;
         forward_propl_rho_y[j][i]=0.0;
         forward_propl_u[j][i]=0.0;
         forward_propl_x[j][i]=0.0;
         forward_propl_y[j][i]=0.0;

         back_prop_rho_x[j][i]=0.0;
         back_prop_rho_y[j][i]=0.0;
         back_prop_u[j][i]=0.0;
         back_prop_x[j][i]=0.0;
         back_prop_y[j][i]=0.0;
   
      }
   }

}                                                         
     
/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

hin=1;
hin1=1;

imat=1;
imat1=1;
imat2=1;
hi=1;

if(RTMOD==0){  

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%10000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /* update of particle velocities */
      update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                         
	if (MYID==0){
		time4=MPI_Wtime();
		time_av_v_update+=(time4-time3);
		if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
	}
                                                   
	/* exchange of particle velocities between PEs */
	exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);
                                                               
	if (MYID==0){
	  time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
	  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
	}                                                                                         	

    if (L)    /* viscoelastic */
    	update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout,
				pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
				K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
    else
   	update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
         	               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  


    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc,0);
   
   /* moment tensor source */
   if ((!CHECKPTREAD)&&(QUELLTYP==5)) 	
   msource(nt,psxx,psyy,psxy,srcpos_loc,signals,nsrc_loc,0);

   if ((FREE_SURF) && (POS[2]==0)){
   	if (L)    /* viscoelastic */
		surface_PML(1, pvx, pvy, psxx, psyy, psxy, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs);
	else      /* elastic */
   		surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);
   }


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	  time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  

	/* store amplitudes at receivers in section-arrays */
	if (SEISMO){
		seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		/*lsamp+=NDT;*/
	}


if(nt==hin1){

    /* save forward wavefields for time-domain inversion */
    if(INVMAT==0){
    
        for (i=1;i<=NX;i=i+IDXI){
	    for (j=1;j<=NY;j=j+IDYI){
		 forward_prop_rho_x[imat1]=pvxp1[j][i];
		 forward_prop_rho_y[imat1]=pvyp1[j][i];
                 imat1++;                                   
	    }
        }   

        for (i=1;i<=NX;i=i+IDXI){ 
	    for (j=1;j<=NY;j=j+IDYI){
	    
	        /* gradients with data integration */
                if(GRAD_FORM==1){
	           forward_prop_x[imat]=psxx[j][i];
	           forward_prop_y[imat]=psyy[j][i];
                }
	    
	        /* gradients without data integration */
	        if((GRAD_FORM==2)||(GRAD_FORM==3)){
                   forward_prop_x[imat]=ux[j][i];
	           forward_prop_y[imat]=uy[j][i];
	        }
	    
	        imat++;
            }
         }
    

        for (i=1;i<=NX;i=i+IDXI){ 
	    for (j=1;j<=NY;j=j+IDYI){
	    
	        /* gradients with data integration */
                if(GRAD_FORM==1){
 	          forward_prop_u[imat2]=psxy[j][i];
	        }

	        /* gradients without data integration */
                if((GRAD_FORM==2)||(GRAD_FORM==3)){
 	          forward_prop_u[imat2]=uxy[j][i];
                }

	        imat2++;
	    
            }
        }
	
    }
   
    if((EPRECOND==1)||(EPRECOND==3)){
      eprecond(Ws,pvx,pvy);
    }

    if((EPRECOND==2)&&(iter==1)){
       
      /* calculate Pseudo-Hessians in time-domain */
      for(i=1;i<=NX;i=i+IDX){
	 for(j=1;j<=NY;j=j+IDY){

            if(INVMAT1==1){

              muss = prho[j][i] * pu[j][i] * pu[j][i];
	      lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;

              if(muss>0.0){
	        hessian_lam[j][i]  += (psxx[j][i] + psyy[j][i])*(psxx[j][i] + psyy[j][i]);
	        hessian_mu[j][i]  += (psxy[j][i]*psxy[j][i]/(muss*muss)) + ((psxx[j][i]+psyy[j][i])*(psxx[j][i]+psyy[j][i])/(4.0*(lamss+muss)*(lamss+muss))) + ((psxx[j][i]-psyy[j][i])*(psxx[j][i]-psyy[j][i])/(4.0*muss*muss));
	        hessian_rho[j][i] += (pvxp1[j][i]*pvxp1[j][i]) + (pvyp1[j][i]*pvyp1[j][i]);
              }
 
            }

	 }
      }

    }
 
    hin++;
    hin1=hin1+DTINV;
                                                             

DTINV_help[nt]=1;

}

   /* save forward wavefields for time-Laplace domain inversion */   
   if(INVMAT==1){
     
     time = (float)(nt*DT);
     tmp = exp(-GAMMA*(time+DT)*(time+DT));

        for (i=1;i<=NX;i=i+IDXI){ 
	    for (j=1;j<=NY;j=j+IDYI){

	        forward_propl_rho_x[j][i]+=pvxp1[j][i]*tmp;
		forward_propl_rho_y[j][i]+=pvyp1[j][i]*tmp;

	        /* gradients with data integration */
                if(GRAD_FORM==1){
 	          forward_propl_u[j][i]+=psxy[j][i]*tmp;
		  forward_propl_x[j][i]+=psxx[j][i]*tmp;
	          forward_propl_y[j][i]+=psyy[j][i]*tmp;
	        }

	        /* gradients without data integration */
                if((GRAD_FORM==2)||(GRAD_FORM==3)){
 	          forward_propl_u[j][i]+=uxy[j][i]*tmp;
		  forward_propl_x[j][i]+=ux[j][i]*tmp;
	          forward_propl_y[j][i]+=uy[j][i]*tmp;
                }
		
            }
        }
	
    }


   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }

   
      
   if (MYID==0){
      time8=MPI_Wtime();
	  time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (forward model) ----------*/

} /* end if(RTMOD==0) */	
	
if ((SEISMO==1)&&(INVMAT==10)){
        
        catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
	catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
	
	if (MYID==0){
	   saveseis_glob(FP,fulldata_vx,fulldata_vy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);  
	}
}

if ((SEISMO==2)&&(INVMAT==10)){  

        catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);

        if (MYID==0){
           saveseis_glob(FP,sectionvx,sectionvy,fulldata_p,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
        }
}

if ((SEISMO==3)&&(INVMAT==10)){  

        catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
        catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);

        if (MYID==0){
           saveseis_glob(FP,sectionvx,sectionvy,sectionp,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
        }

}


if ((SEISMO==4)&&(INVMAT==10)){  

        catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
        catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
        catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
        catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
        catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);

        if (MYID==0){
           saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
        }         
 
}


if(INVMAT<=1){

if (MYID==0){
printf("-------------------  \n");
printf("Calculate residuals  \n");
printf("-------------------  \n");
}

/* calculate L2-Norm and energy ? */
if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){swstestshot=1;}

if ((ntr > 0)&&(HESSIAN==0)){

/* read seismic data from SU file vx */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,1,iter);

if (TIME_FILT){

   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
  /* read synthetic seismic data at time step t0 vx */
  inseis(fprec,ishot,sectionread,ntr_glob,ns,9,iter);

  if (TIME_FILT){

     timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

     if(TIME_FILT==2){ /* band-pass */
       timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
     }

  }
  
  /* calculate vx_mod(t1) - vx_mod(t0) */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvx[h][j]=sectionvx[h][j]-sectionread[recpos_loc[3][i]][j];
   }
   h++;
}
                      
} /* end of TIMELAPSE */

  L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);     
  if(swstestshot==1){energy=calc_energy(sectionvxdata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}

  L2_all_shots=calc_misfit(sectionvxdiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
  energy_all_shots=calc_energy(sectionvxdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);

} /* end QUELLTYPB */

/* read seismic data from SU file vy */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);


if (TIME_FILT){

   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);
 
   if(TIME_FILT==2){ /* apply band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
    /* read synthetic seismic data at time step t0 vy */
    inseis(fprec,ishot,sectionread,ntr_glob,ns,10,iter); 

    if (TIME_FILT){
	
       timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

       if(TIME_FILT==2){ /* apply band-pass */
         timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
       }

    }
   
    /* calculate vy_mod(t1) - vy_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
       for(j=1;j<=ns;j++){
	     sectionvy[h][j]=sectionvy[h][j]-sectionread[recpos_loc[3][i]][j];
       }
       h++;
    }
                               
} /* end of TIMELAPSE */
                               
L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);
if(swstestshot==1){energy=calc_energy(sectionvydata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}

L2_all_shots=calc_misfit(sectionvydiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
energy_all_shots=calc_energy(sectionvydata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);	   	    


} /* end QUELLTYPB */

/* read seismic data from SU file p */
/* --------------------------------- */
if(QUELLTYPB==4){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,11,iter);


if (TIME_FILT){

   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* apply band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionpdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){

    /* read synthetic seismic data at time step t0 vy */
    inseis(fprec,ishot,sectionread,ntr_glob,ns,15,iter);

    if (TIME_FILT){

       timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

       if(TIME_FILT==2){ /* apply band-pass */
         timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
       }

    }

    /* calculate p_mod(t1) - p_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
       for(j=1;j<=ns;j++){
             sectionp[h][j]=sectionp[h][j]-sectionread[recpos_loc[3][i]][j];
       }
       h++;
    }

} /* end of TIMELAPSE */

L2=calc_res(sectionpdata,sectionp,sectionpdiff,sectionpdiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);
if(swstestshot==1){energy=calc_energy(sectionpdata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}

L2_all_shots=calc_misfit(sectionpdiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
energy_all_shots=calc_energy(sectionpdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);


} /* end QUELLTYPB == 4*/


} /* end HESSIAN != 1 */

if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){
       swstestshot=0;
       itestshot+=TESTSHOT_INCR;
}

if ((SEISMO)&&(iter==1)&&(INVMAT<=1)&&(ishot==1)){


   if(QUELLTYPB==1){
   
      catseis(sectionvxdiff, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
      catseis(sectionvydiff, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
      
      if (MYID==0){
         saveseis_glob(FP,fulldata_vx,fulldata_vy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
      }
      
   }
   
   if(QUELLTYPB==4){
   
      catseis(sectionpdiff, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
      
      if (MYID==0){
         saveseis_glob(FP,sectionvx,sectionvy,fulldata_p,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
      }
      
   }   
           
}
                  
	    		    
/* --------------------------------------------------------------------------------------------------------------------- */

/*srcpos_loc_back=splitsrc_back(recpos,&nsrc_loc,ntr_glob);*/
  
for (irec=1;irec<=1;irec+=RECINC){ /* loop over shots at receiver positions */

/*for (irec=20;irec<=20;irec+=1){*/

hin=1;
hin1=1;

    if(MYID==0){
    printf("\n==================================================================================\n");
    printf("\n MYID=%d *****  Starting simulation (backward model) for shot %d of %d  ********** \n",MYID,irec,1);
    printf("\n==================================================================================\n\n");}
    
    /* Distribute multiple source positions on subdomains */
    /* define source positions at the receivers */
    srcpos_loc_back = matrix(1,6,1,ntr);
    for (i=1;i<=ntr;i++){
        srcpos_loc_back[1][i] = (recpos_loc[1][i]);
        srcpos_loc_back[2][i] = (recpos_loc[2][i]);
    }
    ntr1=ntr;
                                    
                

/* initialize wavefield with zero */
if (L){
	zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pr,pp,pq);
}else{	
	zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	
}

/*----------------------  loop over timesteps (backpropagation) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

/*invtimer=NT/DTINV;*/
for (nt=1;nt<=NT;nt++){     

	
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%10000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /*if(MYID==0){  
      printf("Calculate zero lag x-correlation of forward and backpropagated wavefields \n");
      printf("------------------------------------------------------------------------- \n");
      }*/

   /* update of particle velocities */
   /*update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp,
   	       srcpos_loc_back,sectionvxdiff,sectionvydiff,ntr,absorb_coeff, hc, infoout,1);*/

   update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc_back, sectionvxdiff,sectionvydiff,ntr1,absorb_coeff,hc,infoout,1, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x); 

   if (MYID==0){
      time4=MPI_Wtime();
      time_av_v_update+=(time4-time3);
     if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
      }

   /* exchange of particle velocities between PEs */
   exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   if (MYID==0){
      time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
      if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
   }

   if (L)    /* viscoelastic */
    	update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout,
				pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
				K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
   else
   	update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
         	               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  

   /* Backpropagate pressure wavefield */	
   psource(nt,psxx,psyy,srcpos_loc_back,sectionpdiff,ntr1,1);

     if ((FREE_SURF) && (POS[2]==0)){
   	if (L)    /* viscoelastic */
		surface_PML(1, pvx, pvy, psxx, psyy, psxy, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs);
	else      /* elastic */
   		surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);
   }


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
   }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
   }

	/* store amplitudes at receivers in section-arrays */
	/*if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionp, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
	} */ 
   
   /* calculate maximum pressure values at each grid point*/
   if(RTMOD==1){
   
     for (i=1;i<=NX;i=i+IDXI){      
         for (j=1;j<=NY;j=j+IDYI){
   
             if(waveconv_shot[j][i]<(sqrt((psxx[j][i] + psyy[j][i])*(psxx[j][i] + psyy[j][i])))){
	           waveconv_shot[j][i]=sqrt((psxx[j][i] + psyy[j][i])*(psxx[j][i] + psyy[j][i]));
             }
         }
     }
         
   } /*end if(RTMOD==1)*/
   
   if((RTMOD==0)&&(DTINV_help[NT-nt+1]==1)){
    
        if(HESSIAN==0){imat=((nxnyi*(NTDTINV)) - hin*nxnyi)+1;}
	
	/* save backpropagated wavefields for time-domain inversion and partially assemble gradients */                 
        if((HESSIAN==0)&&(INVMAT==0)){
	
	    for (i=1;i<=NX;i=i+IDXI){   
	        for (j=1;j<=NY;j=j+IDYI){ 
                                           
		   	waveconv_rho_shot[j][i]+=(pvxp1[j][i]*forward_prop_rho_x[imat])+(pvyp1[j][i]*forward_prop_rho_y[imat]);
                        waveconv_shot[j][i]+= (forward_prop_x[imat]+forward_prop_y[imat])*(psxx[j][i]+psyy[j][i]);
			
		   	/* mu-gradient with data integration */
		   	if(GRAD_FORM==1){			  

                           if(INVMAT1==1){
			       muss = prho[j][i] * pu[j][i] * pu[j][i];
			       lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
			   }
	           
			   if(INVMAT1==3){
			       muss = pu[j][i];
			       lamss = ppi[j][i]; 
			   } 
	                        
			   if(pu[j][i]>0.0){
			      waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[imat] * psxy[j][i])) 
                                  + ((1.0/4.0) * ((forward_prop_x[imat] + forward_prop_y[imat]) * (psxx[j][i] + psyy[j][i])) / ((lamss+muss)*(lamss+muss)))  
                                  + ((1.0/4.0) * ((forward_prop_x[imat] - forward_prop_y[imat]) * (psxx[j][i] - psyy[j][i])) / (muss*muss));
			   }
			   
                        }

			/* Vs-gradient without data integration (stress-velocity in non-conservative form) */
                        if(GRAD_FORM==2){

                           if(INVMAT1==1){
			       muss = prho[j][i] * pu[j][i] * pu[j][i];
			       lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
			   }
	           
			   if(INVMAT1==3){
			       muss = pu[j][i];
			       lamss = ppi[j][i]; 
			   } 
	                        
			   if(pu[j][i]>0.0){
			      waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[imat] * psxy[j][i])) 
                                  + (((6.0*(lamss*lamss)+4.0*(muss*muss)+8.0*lamss*muss) * (forward_prop_x[imat] * psxx[j][i] + forward_prop_y[imat] * psyy[j][i]))  
                                  - ((3.0*(lamss*lamss)+4.0*lamss*muss) * (forward_prop_y[imat] * psxx[j][i] + forward_prop_x[imat] * psyy[j][i]))) / (2.0*muss*muss*(3.0*lamss+2*muss)*(3*lamss+2*muss));
			   } 
                          
                        }
			
			/* mu-gradient without data integration (stress-velocity with elastic tensor) */
			if(GRAD_FORM==3){

		          waveconv_u_shot[j][i]+= 2.0*((forward_prop_x[imat]*psxx[j][i]) + (forward_prop_y[imat]*psyy[j][i])) + (forward_prop_u[imat] * psxy[j][i]);
                          
                        }
                      		                                                                                                             
		   imat++;
		   }
	    }
		
		
	  }
	  
	  
	  if(EPRECOND==1){
	     eprecond(Wr,pvx,pvy);
	  }
                                                                                                                               
    hin++;
    }
    
    	  /* save backpropagated wavefields for time-Laplace domain inversion */   
          if(INVMAT==1){

             for (i=1;i<=NX;i=i+IDXI){ 
	        for (j=1;j<=NY;j=j+IDYI){
	    
	            back_prop_rho_x[j][i]+=pvxp1[j][i];
		    back_prop_rho_y[j][i]+=pvyp1[j][i];

	            /* gradients with data integration */
                    if(GRAD_FORM==1){
 	              back_prop_u[j][i]+=psxy[j][i];
		      back_prop_x[j][i]+=psxx[j][i];
	              back_prop_y[j][i]+=psyy[j][i];
	            }

	            /* gradients without data integration */
                    if((GRAD_FORM==2)||(GRAD_FORM==3)){
 	              back_prop_u[j][i]+=uxy[j][i];
		      back_prop_x[j][i]+=ux[j][i];
	              back_prop_y[j][i]+=uy[j][i];
                    }
		
                 }
             }
	
          }

     /*fclose(FP2);*/
     
   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }
                                                                                                                                               
   
   if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (backpropagation)----------*/


} /* end of loop over shots */
	
/*if ((ntr > 0) && (SEISMO)){
	saveseis(FP,sectionpdiff,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,0);
}*/


/* partially assemble Laplace-domain gradients */                 
if(INVMAT==1){
	
   for (i=1;i<=NX;i=i+IDXI){   
       for (j=1;j<=NY;j=j+IDYI){ 
                                           
	   waveconv_rho_shot[j][i] = ((back_prop_rho_x[j][i]*forward_propl_rho_x[j][i])+(back_prop_rho_y[j][i]*forward_propl_rho_y[j][i]));
           waveconv_shot[j][i] = (forward_propl_x[j][i]+forward_propl_y[j][i])*(back_prop_x[j][i]+back_prop_y[j][i]);
			
           /* mu-gradient with data integration */
	   if(GRAD_FORM==1){			  

               if(INVMAT1==1){
                 muss = prho[j][i] * pu[j][i] * pu[j][i];
	         lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
	       }
	           
	       if(INVMAT1==3){
	         muss = pu[j][i];
	         lamss = ppi[j][i]; 
	       } 
	                        
	       if(pu[j][i]>0.0){
	         waveconv_u_shot[j][i]= (((1.0/(muss*muss))*(forward_propl_u[j][i] * back_prop_u[j][i])) 
                     + ((1.0/4.0) * ((forward_propl_x[j][i] + forward_propl_y[j][i]) * (back_prop_x[j][i] + back_prop_y[j][i])) / ((lamss+muss)*(lamss+muss)))  
                     + ((1.0/4.0) * ((forward_propl_x[j][i] - forward_propl_y[j][i]) * (back_prop_x[j][i] - back_prop_y[j][i])) / (muss*muss)));
	       }
			   
            }

	    /* Vs-gradient without data integration (stress-velocity in non-conservative form) */
            if(GRAD_FORM==2){

                if(INVMAT1==1){
		  muss = prho[j][i] * pu[j][i] * pu[j][i];
	          lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
	        }
	           
	        if(INVMAT1==3){
	          muss = pu[j][i];
	          lamss = ppi[j][i]; 
                } 
	                        
		if(pu[j][i]>0.0){
		  waveconv_u_shot[j][i] = DT*DT*(((1.0/(muss*muss))*(forward_propl_u[j][i] * back_prop_u[j][i])) 
                      + (((6.0*(lamss*lamss)+4.0*(muss*muss)+8.0*lamss*muss) * (forward_propl_x[j][i] * back_prop_x[j][i] + forward_propl_y[j][i] * back_prop_y[j][i]))  
                      - ((3.0*(lamss*lamss)+4.0*lamss*muss) * (forward_propl_y[j][i] * back_prop_x[j][i] + forward_propl_x[j][i] * back_prop_y[j][i]))) / (2.0*muss*muss*(3.0*lamss+2*muss)*(3*lamss+2*muss)));
	         } 
                          
             }
			
	     /* mu-gradient without data integration (stress-velocity with elastic tensor) */
	     if(GRAD_FORM==3){

	       waveconv_u_shot[j][i]+= DT*DT*(2.0*((forward_propl_x[j][i]*back_prop_x[j][i]) + (forward_propl_y[j][i]*back_prop_y[j][i])) + (forward_propl_u[j][i] * back_prop_u[j][i]));
                          
             }
                      		                                                                                                             
       }
   }
				
}

/* calculate gradient for lambda, Vp or Zp */
/* --------------------------------------- */
	
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){

       /* calculate lambda gradient */           
       waveconv_lam[j][i] = - DT * waveconv_shot[j][i];
		 
       if(INVMAT1==3){

         if(GRAD_FORM==1){
	 
	     muss = pu[j][i];
            lamss = ppi[j][i];
            waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * waveconv_lam[j][i];
	    
          }

          waveconv_shot[j][i] = waveconv_lam[j][i];
	}
		 
        if(INVMAT1==1){

	   /* calculate Vp gradient */ 
	   if(GRAD_FORM==1){
	   
              muss = prho[j][i] * pu[j][i] * pu[j][i];
              lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 *  muss;
	      waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * waveconv_lam[j][i];
	      waveconv_shot[j][i] = 2.0 * ppi[j][i] * prho[j][i] * waveconv_lam[j][i];
	      
           }

           if(GRAD_FORM==2){                             	
	   
	      muss = prho[j][i] * pu[j][i] * pu[j][i];
	      lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 *  muss;
              waveconv_lam[j][i] = (1.0/((3.0*lamss+2.0*muss) * (3.0*lamss+2.0*muss))) * waveconv_lam[j][i];
	      waveconv_shot[j][i] = 2.0 * ppi[j][i] * prho[j][i] * waveconv_lam[j][i]; 
	      
            }
		   
	    if(GRAD_FORM==3){                             	
	    
	       waveconv_shot[j][i] = 2.0 * ppi[j][i] * prho[j][i] * waveconv_lam[j][i];
	       
	    } 		   
		 	
            if(RTMOD==1){
	      waveconv_shot[j][i] = -waveconv_lam[j][i];
	    }
		
	}
		 
        if(INVMAT1==2){
	
	   /* calculate Zp gradient */
           waveconv_shot[j][i] = 2.0 * ppi[j][i] * waveconv_lam[j][i];
	   
	}
		                                                                       
   }
}

/* calculate gradient for mu, Vs or Zs */
/* ----------------------------------- */

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
		 
      /* calculate mu gradient */ 
      waveconv_mu[j][i] = - DT * waveconv_u_shot[j][i];
		 
      if(INVMAT1==1){
		 
         /* calculate Vs gradient */		 
         if((GRAD_FORM==1)||(GRAD_FORM==3)){
	    waveconv_u_shot[j][i] = (- 4.0 * prho[j][i] * pu[j][i] * waveconv_lam[j][i]) + 2.0 * prho[j][i] * pu[j][i] * waveconv_mu[j][i];
         }
		    
         if(GRAD_FORM==2){ 
            waveconv_u_shot[j][i] = (- 4.0 * prho[j][i] * pu[j][i] * waveconv_lam[j][i]) + 2.0 * prho[j][i] * pu[j][i] * waveconv_mu[j][i];                   
         }
		 
      }
		 
      if(INVMAT1==2){
        /* calculate Zs gradient */
        waveconv_u_shot[j][i] = (- 4.0 * pu[j][i] * waveconv_lam[j][i]) + (2.0 * pu[j][i] * waveconv_mu[j][i]);
      }
		 
      if(INVMAT1==3){
        /* calculate u gradient */
        waveconv_u_shot[j][i] = waveconv_mu[j][i];
      }
		                                                                       
   }
}

/* calculate gradient for density */
/* ------------------------------ */
for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){

       /* calculate density gradient rho' */
       if(GRAD_FORM==1){
          waveconv_rho_s[j][i]= - DT * waveconv_rho_shot[j][i];
       }

       if(GRAD_FORM==2){
          waveconv_rho_s[j][i]= DT * waveconv_rho_shot[j][i];
       }
		
       if(GRAD_FORM==3){
          waveconv_rho_s[j][i]= DT * waveconv_rho_shot[j][i];
       }
		 
       if(INVMAT1==1){
          /* calculate density gradient */
          waveconv_rho_shot[j][i] = ((((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * waveconv_lam[j][i]) + (pu[j][i] * pu[j][i] * waveconv_mu[j][i]) + waveconv_rho_s[j][i]);
       }
		 
       if(INVMAT1==3){
          /* calculate density gradient */
          waveconv_rho_shot[j][i] = waveconv_rho_s[j][i];
       }
	 
    }
}

if((EPRECOND==1)||(EPRECOND==3)){
  /* calculate energy weights */
  eprecond1(We,Ws,Wr);
      
  /* scale gradient with energy weights*/
  for(i=1;i<=NX;i=i+IDX){
      for(j=1;j<=NY;j=j+IDY){

             waveconv_shot[j][i] = waveconv_shot[j][i]/(We[j][i]*C_vp*C_vp);
	     waveconv_u_shot[j][i] = waveconv_u_shot[j][i]/(We[j][i]*C_vs*C_vs);
             if(C_vs==0.0){waveconv_u_shot[j][i] = 0.0;}
	     waveconv_rho_shot[j][i] = waveconv_rho_shot[j][i]/(We[j][i]*C_rho*C_rho);

      }
  }
}

if (SWS_TAPER_CIRCULAR_PER_SHOT){    /* applying a circular taper at the source position to the gradient of each shot */
	
	/* applying the preconditioning */
	taper_grad_shot(waveconv_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	taper_grad_shot(waveconv_rho_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	taper_grad_shot(waveconv_u_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	
} /* end of SWS_TAPER_CIRCULAR_PER_SHOT == 1 */

for(i=1;i<=NX;i=i+IDX){
	for(j=1;j<=NY;j=j+IDY){
		waveconv[j][i] += waveconv_shot[j][i];
		waveconv_rho[j][i] += waveconv_rho_shot[j][i];
		waveconv_u[j][i] += waveconv_u_shot[j][i];
	}
}
   
} /* end of invtype == 1*/

} /* end of invmat==10 */

if(N_STREAMER>0){

   if (SEISMO) free_imatrix(recpos,1,3,1,ntr_glob);

   if ((ntr>0) && (SEISMO)){

           free_imatrix(recpos_loc,1,3,1,ntr);
           recpos_loc = NULL;
 
           switch (SEISMO){
           case 1 : /* particle velocities only */
                   free_matrix(sectionvx,1,ntr,1,ns);
                   free_matrix(sectionvy,1,ntr,1,ns);
                   sectionvx=NULL;
                   sectionvy=NULL;
                   break;
            case 2 : /* pressure only */
                   free_matrix(sectionp,1,ntr,1,ns);
                   break;
            case 3 : /* curl and div only */
                   free_matrix(sectioncurl,1,ntr,1,ns);
                   free_matrix(sectiondiv,1,ntr,1,ns);
                   break;
            case 4 : /* everything */
                   free_matrix(sectionvx,1,ntr,1,ns);
                   free_matrix(sectionvy,1,ntr,1,ns);
                   free_matrix(sectionp,1,ntr,1,ns);
                   free_matrix(sectioncurl,1,ntr,1,ns);
                   free_matrix(sectiondiv,1,ntr,1,ns);
                   break;

            }

   }

   free_matrix(sectionread,1,ntr_glob,1,ns);
   free_ivector(recswitch,1,ntr);
   
   if((QUELLTYPB==1)||(QUELLTYPB==3)){
      free_matrix(sectionvxdata,1,ntr,1,ns);
      free_matrix(sectionvxdiff,1,ntr,1,ns);
      free_matrix(sectionvxdiffold,1,ntr,1,ns);
   }
   
   if((QUELLTYPB==1)||(QUELLTYPB==2)){   
      free_matrix(sectionvydata,1,ntr,1,ns);
      free_matrix(sectionvydiff,1,ntr,1,ns);
      free_matrix(sectionvydiffold,1,ntr,1,ns);
   }
   
   if(QUELLTYPB==4){   
      free_matrix(sectionpdata,1,ntr,1,ns);
      free_matrix(sectionpdiff,1,ntr,1,ns);
      free_matrix(sectionpdiffold,1,ntr,1,ns);
   }
   
 
}

nsrc_loc=0;

} /* end of loop over shots (forward and backpropagation) */   

/* save Hessian */
if(EPRECOND==2){
  hessian_out(hessian_lam,hessian_mu,hessian_rho,ppi,pu,prho);
}

/* calculate L2 norm of all CPUs*/
L2sum = 0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
energy_sum = 0.0;
MPI_Allreduce(&energy,&energy_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
L2sum_all_shots = 0.0;
MPI_Allreduce(&L2_all_shots,&L2sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
energy_sum_all_shots = 0.0;
MPI_Allreduce(&energy_all_shots,&energy_sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

if((LNORM==2)&&(INVMAT!=1)){
     L2t[1]=L2sum/energy_sum;
     L2t[4]=L2sum/energy_sum;
}
else{L2t[1]=L2sum;
     L2t[4]=L2sum;}

/*if(MYID==0){
	printf("L2sum: %e\n", L2sum);
	printf("energy_sum: %e\n\n", energy_sum);
	printf("L2sum_all_shots: %e\n", L2sum_all_shots);
	printf("energy_sum_all_shots: %e\n\n", energy_sum_all_shots);}*/

if(GRAVITY==2){

  /* save seismic L2-norm of seismic data residuals */
  L2sum = L2t[1];

  /* global density model */
  rho_grav =  matrix(1,NYG,1,NXG);
  rho_grav_ext =  matrix(1,nygrav,1,nxgrav);

  /* model gravity data */
  /* save current density model */
  sprintf(jac_grav,"%s_tmp.rho.%i%i",JACOBIAN,POS[1],POS[2]);
  FP_GRAV=fopen(jac_grav,"wb");

  for (i=1;i<=NX;i=i+IDX){
      for (j=1;j<=NY;j=j+IDY){
          fwrite(&prho[j][i],sizeof(float),1,FP_GRAV);
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

  /* save L2_grav value at iter=1 */  
  if(iter==1){L2_GRAV_IT1 = L2_grav;}

  /* calculate lambda_grav */
  LAM_GRAV = GAMMA_GRAV * (L2sum/L2_GRAV_IT1); 

  /* add gravity penalty term to the seismic objective function */
  L2t[1]+=LAM_GRAV * L2_grav;
  L2t[4]+=LAM_GRAV * L2_grav;

  /* calculate gravity gradient */
  grav_grad(ngrav,gravpos,grad_grav,gz_res);
  
  /* apply taper functions */
  if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(grad_grav,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}
   
  if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(grad_grav,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

  /* output of the gravity gradients after preconditioning */
  sprintf(jac,"%s_grav.%i%i",JACOBIAN,POS[1],POS[2]);
  FP_GRAV=fopen(jac,"wb");       

  for (i=1;i<=NX;i=i+IDX){
      for (j=1;j<=NY;j=j+IDY){
          fwrite(&grad_grav[j][i],sizeof(float),1,FP_GRAV);
      }
  }

  fclose(FP_GRAV);

  MPI_Barrier(MPI_COMM_WORLD);        

  /* merge model file */
  sprintf(jac,"%s_grav",JACOBIAN);          
  if (MYID==0) mergemod(jac,3); 

  /* free memory */
  free_matrix(rho_grav,1,NYG,1,NXG);
  free_matrix(rho_grav_ext,1,nygrav,1,nxgrav);
  

}
  
if(INVMAT<=1){

   /* Interpolate missing spatial gradient values in case IDXI > 1 || IDXY > 1 */
   /* ------------------------------------------------------------------------ */

   if((IDXI>1)||(IDYI>1)){

      interpol(IDXI,IDYI,waveconv,1);
      interpol(IDXI,IDYI,waveconv_u,1);
      interpol(IDXI,IDYI,waveconv_rho,1);

   }


   /* Add gravity gradient to FWI density gradient */
   /* -------------------------------------------- */
	
   if(GRAVITY==2){
		 		 
     /* calculate maximum values of waveconv_rho and grad_grav */
     FWImax = 0.0;
     GRAVmax = 0.0;
	
     for (i=1;i<=NX;i++){
        for (j=1;j<=NY;j++){
		
	    if(fabs(waveconv_rho[j][i])>FWImax){FWImax=fabs(waveconv_rho[j][i]);}
	    if(fabs(grad_grav[j][i])>GRAVmax){GRAVmax=fabs(grad_grav[j][i]);}
		
        }
     }
	
     MPI_Allreduce(&FWImax,&FWImax_all,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
     MPI_Allreduce(&GRAVmax,&GRAVmax_all,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
		
     LAM_GRAV_GRAD = GAMMA_GRAV * FWImax_all/GRAVmax_all;
		 
     /* add gravity gradient to seismic density gradient */
     for (i=1;i<=NX;i++){
        for (j=1;j<=NY;j++){
			
            waveconv_rho[j][i] += LAM_GRAV_GRAD * grad_grav[j][i];
				
        }
     }
		
   }

}

if((HESSIAN==0)&&(GRAD_METHOD==1)&&(INVMAT<=1)){
  PCG(waveconv, taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, C_vp, gradp, nfstart_jac, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho);
}

if((HESSIAN==0)&&(GRAD_METHOD==2)&&(INVMAT<=1)){
  LBFGS1(taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, nfstart_jac, waveconv, C_vp, gradp, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho, y_LBFGS, s_LBFGS, rho_LBFGS, alpha_LBFGS, ppi, pu, prho, nxnyi, q_LBFGS, r_LBFGS, beta_LBFGS, LBFGS_pointer, NLBFGS, NLBFGS_vec);
}

opteps_vp=0.0;
opteps_vs=0.0;
opteps_rho=0.0;

/* ============================================================================================================================*/
/* =============================================== test loop L2 ===============================================================*/
/* ============================================================================================================================*/

if(RTM==0){ /* only if RTM==0 */

if((INVMAT<=1) && (HESSIAN==0)){

/* set min_iter_help to initial global value of MIN_ITER */
if(iter==1){min_iter_help=MIN_ITER;}

/* Estimate optimum step length ... */

/* ... by line search (parabolic fitting) */
eps_scale = step_length_est(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
                nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x,
                a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
                sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,sectionpdata,sectionpdiff,sectionpdiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,
		bufferrig_to_lef,buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec,pr,pp,pq,fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,ptaup,ptaus,etajm,peta,etaip,
                ptausipjp,recpos,&step1,&step3,C_vp,gradg,FC,nxgrav,nygrav,ngrav,gravpos,gz_mod,NZGRAV,recswitch,FP,ntr_loc);

/* ... by line search (parabolic fitting+safeguard) */
/*eps_scale = step_length_est1(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	        nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	        a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	        sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	        buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec,pr,pp,pq,fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,ptaup,ptaus,etajm,peta,etaip,
                ptausipjp,recpos,&step1,&step3,C_vp,gradg,FC);*/


/* no model update due to steplength estimation failed or update with the smallest steplength if the number of iteration is smaller than the minimum number of iteration per
frequency MIN_ITER */
if((iter>min_iter_help)&&(step1==0)){ 
	eps_scale=0.0;
	opteps_vp=0.0;
}
else{
	opteps_vp=eps_scale;
}


if(MYID==0){
printf("MYID = %d \t opteps_vp = %e \t opteps_vs = %e \t opteps_rho = %e \n",MYID,opteps_vp,opteps_vs,opteps_rho);
printf("MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \t L2t[4] = %e \n",MYID,L2t[1],L2t[2],L2t[3],L2t[4]);
printf("MYID = %d \t epst1[1] = %e \t epst1[2] = %e \t epst1[3] = %e \n",MYID,epst1[1],epst1[2],epst1[3]);}

if(MYID==0){
if (TIME_FILT==0){
	fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %d \n",opteps_vp,epst1[1],epst1[2],epst1[3],L2t[1],L2t[2],L2t[3],L2t[4],nstage);}
else{
	fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %f \t %f \t %d \n",opteps_vp,epst1[1],epst1[2],epst1[3],L2t[1],L2t[2],L2t[3],L2t[4],FC_START,FC,nstage);}}


/* saving history of final L2*/
L2_hist[iter]=L2t[4];
s=0;


/* calculate optimal change in the material parameters */
eps_true=calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,INVMAT,eps_scale,0,nfstart);

} /* end of if(INVMAT!=4) */

if ((INVMAT<=1)&&(MODEL_FILTER)){
/* smoothing the velocity models vp and vs */
smooth_model(ppi,pu,prho,iter);
}

} /* only if RTM==0 */

if(MYID==0){	
/*	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"STATISTICS FOR ITERATION STEP %d \n",iter);
	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n");*/
/*	fprintf(FPL2,"Low-pass filter at %e Hz\n",freq);
	fprintf(FPL2,"----------------------------------------------\n");
*/	/*fprintf(FPL2,"L2 at iteration step n = %e \n",L2);*/
/*        fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",EPSILON,EPSILON_u,EPSILON_rho,L2t[4],betaVp,betaVs,betarho,sqrt(C_vp));*/

	/*fprintf(FPL2,"----------------------------------------------\n");*/
/*	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n\n\n");*/
}

if(INVMAT!=10){
if(MYID==0){
  fclose(FPL2);
}}

if(iter==nfstart){
nfstart = nfstart + nf;
}

if(iter==nfstart_jac){
nfstart_jac = nfstart_jac + nf_jac;
}

if (iter>min_iter_help){

float diff=0.0, pro=PRO;

/* calculating differnce of the actual L2 and before two iterations, dividing with L2_hist[iter-2] provide changing in procent*/
diff=fabs((L2_hist[iter-2]-L2_hist[iter])/L2_hist[iter-2]);
	
	if((diff<=pro)||(step3==1)){
        
        	/* output of the model at the end of given corner frequency */
        	model_freq_out(ppi,prho,pu,nstage,FC);
		s=1;
		min_iter_help=0;
		min_iter_help=iter+MIN_ITER;
		iter=0;

        	if(GRAD_METHOD==2){
	  		zero_LBFGS(NLBFGS, NLBFGS_vec, y_LBFGS, s_LBFGS, q_LBFGS, r_LBFGS, alpha_LBFGS, beta_LBFGS, rho_LBFGS);
          		LBFGS_pointer = 1;  
		}

        	if(MYID==0){
			if(step3==1){
			        printf("\n Steplength estimation failed step3=%d \n Changing to next FWI stage \n",step3);
			}
			else{
  				printf("\n Reached the abort criterion of pro=%e and diff=%e \n Changing to next FWI stage \n",pro,diff);
			}
	
		}
		break;
	}
}

iter++;
iter_true++;

if(EPRECOND==2){
  free_matrix(hessian_lam,-nd+1,NY+nd,-nd+1,NX+nd);
  free_matrix(hessian_mu,-nd+1,NY+nd,-nd+1,NX+nd);
  free_matrix(hessian_rho,-nd+1,NY+nd,-nd+1,NX+nd);
}

/* ====================================== */
} /* end of fullwaveform iteration loop*/
/* ====================================== */

} /* End of FWI-workflow loop */

if (CHECKPTWRITE){
	if (MYID==0){
		time3=MPI_Wtime();
 		fprintf(FP," Saving wavefield to check-point file %s \n",CHECKPTFILE);	
	}
	
	save_checkpoint(-1, NX+2, -1, NY+2, pvx, pvy, psxx, psyy, psxy);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0){
		time4=MPI_Wtime();
      		fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
	}
}

/* deallocation of memory */
free_matrix(psxx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psxy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psyy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ux,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uxy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uyx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Vp0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Vs0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Rho0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uttx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(utty,-nd+1,NY+nd,-nd+1,NX+nd);


free_matrix(prho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prhonp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prip,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pripnp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prjpnp1,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(ppi,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ppinp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(punp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(vpmat,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(puipjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_lam,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconvtmp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(wcpart,1,3,1,3);

free_matrix(wavejac,-nd+1,NY+nd,-nd+1,NX+nd);

if(FW>0){

free_vector(d_x,1,2*FW);
free_vector(K_x,1,2*FW);
free_vector(alpha_prime_x,1,2*FW);
free_vector(a_x,1,2*FW);
free_vector(b_x,1,2*FW);

free_vector(d_x_half,1,2*FW);
free_vector(K_x_half,1,2*FW);
free_vector(alpha_prime_x_half,1,2*FW);
free_vector(a_x_half,1,2*FW);
free_vector(b_x_half,1,2*FW);

free_vector(d_y,1,2*FW);
free_vector(K_y,1,2*FW);
free_vector(alpha_prime_y,1,2*FW);
free_vector(a_y,1,2*FW);
free_vector(b_y,1,2*FW);

free_vector(d_y_half,1,2*FW);
free_vector(K_y_half,1,2*FW);
free_vector(alpha_prime_y_half,1,2*FW);
free_vector(a_y_half,1,2*FW);
free_vector(b_y_half,1,2*FW);

free_matrix(psi_sxx_x,1,NY,1,2*FW);
free_matrix(psi_syy_y,1,2*FW,1,NX);
free_matrix(psi_sxy_x,1,NY,1,2*FW);
free_matrix(psi_sxy_y,1,2*FW,1,NX);
free_matrix(psi_vxx,1,NY,1,2*FW);
free_matrix(psi_vxxs,1,NY,1,2*FW);
free_matrix(psi_vyy,1,2*FW,1,NX);
free_matrix(psi_vxy,1,2*FW,1,NX);
free_matrix(psi_vyx,1,NY,1,2*FW);

}

/*absorb_coeff=  matrix(1,NY,1,NX);
free_matrix(absorb_coeff,1,NY,1,NX);*/
free_matrix(taper_coeff,1,NY,1,NX);

free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
free_matrix(buffertop_to_bot,1,NX,1,fdo3);
free_matrix(bufferbot_to_top,1,NX,1,fdo3);

free_vector(hc,0,6);

if(INVMAT<=1){

free_matrix(gradg,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradg_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho_s,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradg_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_mu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_u_shot,-nd+1,NY+nd,-nd+1,NX+nd);

}

if(INVMAT==0){

free_vector(forward_prop_x,1,NY*NX*NT);
free_vector(forward_prop_y,1,NY*NX*NT);
free_vector(forward_prop_rho_x,1,NY*NX*NT);
free_vector(forward_prop_rho_y,1,NY*NX*NT);
free_vector(forward_prop_u,1,NY*NX*NT);

}

if(INVMAT==1){

free_matrix(forward_propl_x,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(forward_propl_y,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(forward_propl_rho_x,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(forward_propl_rho_y,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(forward_propl_u,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(back_prop_x,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(back_prop_y,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(back_prop_rho_x,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(back_prop_rho_y,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(back_prop_u,-nd+1,NY+nd,-nd+1,NX+nd);

}


if (nsrc_loc>0){	
	free_matrix(signals,1,nsrc_loc,1,NT);
	/*dsignals=fmatrix(1,nsrc,1,NT);
	free_matrix(dsignals,1,nsrc_loc,1,NT);*/
	free_matrix(srcpos_loc,1,8,1,nsrc_loc);
	free_matrix(srcpos_loc_back,1,6,1,nsrc_loc);
}		   

 /* free memory for global source positions */
 free_matrix(srcpos,1,8,1,nsrc);

 /* free memory for source position definition */
 free_matrix(srcpos1,1,8,1,1);
 
 /* free memory for abort criterion */
 free_vector(L2_hist,1,1000);
 		
 free_vector(L2t,1,4);
 free_vector(epst1,1,3);
 free_vector(epst2,1,3);
 free_vector(epst3,1,3);

 if(N_STREAMER==0){

    if (SEISMO) free_imatrix(recpos,1,3,1,ntr_glob);

    if ((ntr>0) && (SEISMO)){

            free_imatrix(recpos_loc,1,3,1,ntr);
            recpos_loc = NULL;
 
            switch (SEISMO){
            case 1 : /* particle velocities only */
                    free_matrix(sectionvx,1,ntr,1,ns);
                    free_matrix(sectionvy,1,ntr,1,ns);
                    sectionvx=NULL;
                    sectionvy=NULL;
                    break;
             case 2 : /* pressure only */
                    free_matrix(sectionp,1,ntr,1,ns);
                    break;
             case 3 : /* curl and div only */
                    free_matrix(sectioncurl,1,ntr,1,ns);
                    free_matrix(sectiondiv,1,ntr,1,ns);
                    break;
             case 4 : /* everything */
                    free_matrix(sectionvx,1,ntr,1,ns);
                    free_matrix(sectionvy,1,ntr,1,ns);
                    free_matrix(sectionp,1,ntr,1,ns);
                    free_matrix(sectioncurl,1,ntr,1,ns);
                    free_matrix(sectiondiv,1,ntr,1,ns);
                    break;

             }

    }

    free_matrix(sectionread,1,ntr_glob,1,ns);
    free_ivector(recswitch,1,ntr);
    
    if((QUELLTYPB==1)||(QUELLTYPB==3)){
       free_matrix(sectionvxdata,1,ntr,1,ns);
       free_matrix(sectionvxdiff,1,ntr,1,ns);
       free_matrix(sectionvxdiffold,1,ntr,1,ns);
    }

    if((QUELLTYPB==1)||(QUELLTYPB==2)){    
       free_matrix(sectionvydata,1,ntr,1,ns);
       free_matrix(sectionvydiff,1,ntr,1,ns);
       free_matrix(sectionvydiffold,1,ntr,1,ns);
    }
    
    if(QUELLTYPB==4){    
       free_matrix(sectionpdata,1,ntr,1,ns);
       free_matrix(sectionpdiff,1,ntr,1,ns);
       free_matrix(sectionpdiffold,1,ntr,1,ns);
    }
    
 }

 if(SEISMO){
  free_matrix(fulldata,1,ntr_glob,1,NT); 
 }

 if(SEISMO==1){
  free_matrix(fulldata_vx,1,ntr_glob,1,NT);
  free_matrix(fulldata_vy,1,ntr_glob,1,NT);
 }

 if(SEISMO==2){
  free_matrix(fulldata_p,1,ntr_glob,1,NT);
 } 
 
 if(SEISMO==3){
  free_matrix(fulldata_curl,1,ntr_glob,1,NT);
  free_matrix(fulldata_div,1,ntr_glob,1,NT);
 }

 if(SEISMO==4){
  free_matrix(fulldata_vx,1,ntr_glob,1,NT);
  free_matrix(fulldata_vy,1,ntr_glob,1,NT);
  free_matrix(fulldata_p,1,ntr_glob,1,NT); 
  free_matrix(fulldata_curl,1,ntr_glob,1,NT);
  free_matrix(fulldata_div,1,ntr_glob,1,NT);
 }

 free_ivector(DTINV_help,1,NT);
 
 /* free memory for viscoelastic modeling variables */
 if (L) {
		free_f3tensor(pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_matrix(ptaus,-nd+1,NY+nd,-nd+1,NX+nd);
		free_matrix(ptausipjp,-nd+1,NY+nd,-nd+1,NX+nd);
		free_matrix(ptaup,-nd+1,NY+nd,-nd+1,NX+nd);
		free_vector(peta,1,L);
		free_vector(etaip,1,L);
		free_vector(etajm,1,L);
		free_vector(bip,1,L);
		free_vector(bjm,1,L);
		free_vector(cip,1,L);
		free_vector(cjm,1,L);
		free_matrix(f,-nd+1,NY+nd,-nd+1,NX+nd);
		free_matrix(g,-nd+1,NY+nd,-nd+1,NX+nd);
		free_matrix(fipjp,-nd+1,NY+nd,-nd+1,NX+nd);
		free_f3tensor(dip,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(d,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		free_f3tensor(e,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
}

if(GRAVITY){

  free_matrix(gravpos,1,2,1,ngrav);
  free_vector(gz_mod,1,ngrav);
  free_vector(gz_res,1,ngrav);

  if(GRAVITY==2){
    free_matrix(grad_grav,1,NY,1,NX);
  }

}
 
/* de-allocate buffer for messages */
MPI_Buffer_detach(buff_addr,&buffsize);

/*for (ii=0;ii<=3;ii++){
	MPI_Request_free(&req_send[ii]);
	MPI_Request_free(&req_rec[ii]);}*/
	
	
/* merge snapshot files created by the PEs into one file */
/* if ((SNAP) && (MYID==0)){ 
	snapmerge(nsnap);
}
*/

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0){
	fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
	fprintf(FP," CPU time of program per PE: %li seconds.\n",clock()/CLOCKS_PER_SEC);
	time8=MPI_Wtime();
	fprintf(FP," Total real time of program: %4.2f seconds.\n",time8-time1);
	time_av_v_update=time_av_v_update/(double)NT;
	time_av_s_update=time_av_s_update/(double)NT;
	time_av_v_exchange=time_av_v_exchange/(double)NT;
	time_av_s_exchange=time_av_s_exchange/(double)NT;
	time_av_timestep=time_av_timestep/(double)NT;
	fprintf(FP," Average times for \n");
	fprintf(FP," velocity update:  \t %5.3f seconds  \n",time_av_v_update);
	fprintf(FP," stress update:  \t %5.3f seconds  \n",time_av_s_update);
	fprintf(FP," velocity exchange:  \t %5.3f seconds  \n",time_av_v_exchange);
	fprintf(FP," stress exchange:  \t %5.3f seconds  \n",time_av_s_exchange);
	fprintf(FP," timestep:  \t %5.3f seconds  \n",time_av_timestep);
		
}

fclose(FP);

MPI_Finalize();
return 0;	

}/*main*/
