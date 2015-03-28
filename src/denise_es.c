/*------------------------------------------------------------------------
 *  Parallel 2-D-Viscoelastic Finite Difference Modelling 
 *
 *
 *  Copyright (c) T. Bohlen
 *  
 *
 * 
 * In case of questions contact the author:                       
 *	Dr. Thomas Bohlen, Kiel University, Institute of Geoscience,
 *	Otto-Hahn-Platz 1, D-24098 Kiel, Germany, ph: +49 431 880 4648,
 *	fax: +49 431 880 4432, mailto:tbohlen@geophysik.uni-kiel.de,
 *	Homepage: http://www.geophysik.uni-kiel.de/~tbohlen
 *
 *
 *  Do not freely distribute the program. In case of interest of other
 *  people in obtaining the source code please refer them to the author.
 *  I would like to keep track where the program is used and modified.
 *  If you show modelling results in a paper or presentation please give a reference
 *  to the following paper:
 *  Bohlen, 2002, Parallel 3-D viscoelastic finite-difference seismic modelling,
 *  Computers and Geociences, 28, 887-899.
 *  
 *  Thank you for your co-operation, 
 *  Thomas Bohlen
 

	History of major modifications:

	Febr. 2002	Version 1.0	original parallel implementation
					T. Bohlen
					
	Febr. 2002  Version 1.1
					- fixed bugs which occured when reading
					  source signals from ASCII file
					- ported the coded to CRAY T3E system (changed outseis.c
					   to output correct SU format on CRAY T3E)
					- changed storage of density in field prho[j][i]:
					  it does now store density instead of 1/density ! 
					T. Bohlen

	April 2002  Version 1.5
					- dipping incident plane wave implemented
					T. Bohlen
					
	May 2002  Version 1.6
					- improved merging of snapshots (no memory allocation
					  required), output format of snapshot data was modified
					T. Bohlen

	July 2002  Version 1.7
					- snapshots files can now be merged to one data file using
					  program snapmerge. Usage: 'snapmerge < fdveps.inp' after
					  fdveps has finished. The merged data file can be visualized
					  with the SU program xmovie.
					- version for anisotropic media implemeted and tested.
					T. Bohlen

	June 2003 Version 1.9		- modification of  locations of wavefield and material parameters
	   				  on the standard staggered grid according to Vossen et al.,2002,
					  Geophysics, 67, No. 2, 618-624. In the new implementation only
					  density and shear slowness are averaged. This distribution has
					  higher accuracy in the presence of
					  fluid/solid contrasts (seafloor) compared to the
					  distribution I have used previously (Bohlen, Phd thesis,
					  1998)

	December, 2004, Version 2.0	- Since Revision 2.0 update with CVS

 *  ----------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */


int main(int argc, char **argv){
/* variables in main */
int ns, nseismograms=0, nt, nd, fdo3, j, i, ii, jj, shotid, recid, k, nc, iter, h, infoout, SHOTINC, INVMAT, TIMEWIN, test_eps, lq, iq, jq, hin, hin1;
int DTINV, nxny, nxnyi, imat, imat1, imat2, IDXI, IDYI, hi, NTST, NTSTI, partest, SPATFILTER, FREQFILT;
int lsnap, nsnap=0, lsamp=0, buffsize, invtime, invtimer, sws, snapseis, snapseis1;
int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, ishot, nshots, Lcount, LNORM, itest, Lcountsum;

float memdyn, memmodel, memseismograms, membuffer, memtotal, dngn, fphi, sum, avggrad, beta, betan, betaz, betaLog, betaVp, betaVs, betarho, eps_scale, L2old;
float fac1, fac2, wavefor, waverecipro, dump, dump1, epsilon, gradsign, mun, eps1, gradplastiter, gradglastiter, gradclastiter, betar, sig_max, sig_max1;
float signL1, RMS, opteps_vp, opteps_vs, opteps_rho, Vs, Vp, Vp_avg, C_vp, Vs_avg, C_vs, Cd, rho_avg, C_rho, Vs_sum, Vp_sum, rho_sum, Zp, Zs;
float freqshift, dfreqshift, memfwt, memfwt1, memfwtdata;
char *buff_addr, ext[10], *fileinp;
char wave_forward[225], wave_recipro[225], wave_conv[225], jac[225], jac2[225], jacsum[225], dwavelet[225] ;

double 	time1, time2, time3, time4, time5, time6, time7, time8,
	time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0, 
	time_av_s_exchange=0.0, time_av_timestep=0.0;
	
float L2, L2sum, *L2t, alphanomsum, alphanom, alphadenomsum, alphadenom, scaleamp ,sdummy, lamr; 	

float  **  psxx, **  psxy, **  psyy, ** ux, ** uy, ** uxy, ** uyx, ** Vp0, ** uttx, ** utty, ** Vs0, ** Rho0;
float  **  pvx, **  pvy, **waveconv, **waveconv_lam, **waveconv_mu, **waveconv_rho, **waveconv_rho_s, **waveconv_u, **waveconvtmp, **wcpart, **wavejac;
float  **  pvxp1, **  pvyp1, **  pvxm1, **  pvym1;
float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u, ** hessian;
float  **  prho,**  prhonp1, **prip=NULL, **prjp=NULL, **pripnp1=NULL, **prjpnp1=NULL, **  ppi, **  pu, **  punp1, **  puipjp, **  ppinp1;
float  **  vpmat, *forward_prop_x, *forward_prop_y, *forward_prop_rho_x, *forward_prop_u, *forward_prop_rho_y;

float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
	** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvxdiff=NULL, ** sectionvxdiffold=NULL, ** sectionvydiffold=NULL,
	** sectionvydiff=NULL, ** sectionvydata=NULL, ** sectionpn=NULL, ** sectionread=NULL;
float  **  absorb_coeff, ** taper_coeff, * epst, * picked_times;
float  ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, **srcpos_loc_back=NULL, ** signals=NULL, *hc=NULL, ** dsignals=NULL;
int   **recpos=NULL, ** recpos_loc=NULL;

float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top; 

/* allocate variables for ES */
double drand48(void);
void srand48(long seedval);
float gene, ** matmod, ** matmod1, L2big, minL2;
float minvp, maxvp, minvs, maxvs, minrho, maxrho, dpar1, dpar2, dpar;
int nchild, nparents, ncptot, nparameter;
int is, nmod, nsmod, hs, hes, nrdad, nrmom, minind;
char ES_model[STRING_SIZE];

FILE *fprec, *FP2, *FP3, *FP4, *FP5, *FPL2, *FP6, *FP7 ;
	
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
read_par(FP);
exchange_par();


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
ns=iround(NT/NDT);           /* number of samples per trace */
lsnap=iround(TSNAP1/DT);      /* first snapshot at this timestep */
lsamp=NDT;


/* output of parameters to log-file or stdout */
if (MYID==0)
   write_par(FP);



	
/* NXG, NYG denote size of the entire (global) grid */
NXG=NX;
NYG=NY;

/* In the following, NX and NY denote size of the local grid ! */
NX = IENDX;
NY = IENDY;


if (SEISMO){
   recpos=receiver(FP, &ntr);
   recpos_loc = splitrec(recpos,&ntr_loc, ntr);
   ntr_glob=ntr;
   ntr=ntr_loc;
}



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
	


/*allocate memory for dynamic, static and buffer arrays */
fac1=(NX+FDORDER)*(NY+FDORDER);
fac2=sizeof(float)*pow(2.0,-20.0);

nd = FDORDER/2 + 1;
fdo3 = 2*nd;

memdyn=5.0*fac1*fac2;
memmodel=6.0*fac1*fac2;
memseismograms=nseismograms*ntr*ns*fac2;

memfwt=5.0*(NX/2+FDORDER)*(NY/2+FDORDER)*(NT/2)*fac2;
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

srand48(NXG*NYG*NT);
   

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
wcpart = matrix(1,3,1,3);
wavejac = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
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

/* Inversion of material parameter ? */
/* INVMAT==0 All */
/* INVMAT==1 Paramter 1  Vp, Zp  */
/* INVMAT==2 Parameter 2 density */
/* INVMAT==3 Parameter 3 Vs, Zs  */
/* INVMAT==4 source wavelet      */
INVMAT=10;

/* Time windowing  */
/* TIMEWIN==0 OFF */
/* TIMEWIN==1 ON  */
/* TIMEWIN==2 Time Window with STA/LTA Picker*/
TIMEWIN=0;

/* use only every DTINV time sample for the inversion */
DTINV=2;

/* save every IDXI and IDYI spatial point during the forward modelling */
IDXI=1;
IDYI=1;

/* SPATFILTER = 1 activate spatial wavelengthfilter */
/* SPATFILTER = 0 deactivate spatial wavelengthfilter */
SPATFILTER=0;

/* FREQFILT = 1 activate frequency filter */
/* FREQFILT = 0 deactivate frequency filter */
FREQFILT=0;

lamr=0.01; /* Marquardt factor */

NTST=20;
NTSTI=NTST/DTINV;

nxny=NX*NY;
nxnyi=(NX/IDXI)*(NY/IDYI);

/* define parameters for ES code*/
nchild=45;                  /* number of children */
nparents=10;                /* number of parents */
ncptot=nchild+nparents;    /* sum of children and parents */ 
nparameter=6+1;            /* number of model parameters + value of residual energy */


/* define parents of the 1st generation */
minvp = 950.0;
maxvp = 4500.0;

minvs = 800.0;
maxvs = 4500.0;
maxrho = 2500.0;
 
matmod=matrix(1,ncptot,1,nparameter);
matmod1=matrix(1,ncptot,1,nparameter);

/*matmod[1][1] = 1496.0;
matmod[1][2] = 1300.0;
matmod[1][3] = 1673.0;

matmod[1][4] = 0.3871;
matmod[1][5] = 0.5200;
matmod[1][6] = 0.5328;*/

matmod[1][1] = 800.0;
matmod[1][2] = 800.0;
matmod[1][3] = 1000.0;

matmod[1][4] = 0.2;
matmod[1][5] = 0.2;
matmod[1][6] = 0.2;


/*matmod[1][1] = 1519.390747;
matmod[1][2] = 1293.091187;
matmod[1][3] = 1691.588867;

matmod[1][4] = 0.347294;
matmod[1][5] = 0.312897;
matmod[1][6] = 0.344851;*/


matmod[2][1] = 1000.0;
matmod[2][2] = 1000.0;
matmod[2][3] = 1200.0;

matmod[2][4] = 0.25;
matmod[2][5] = 0.25;
matmod[2][6] = 0.25;

matmod[3][1] = 1200.0;
matmod[3][2] = 1200.0;
matmod[3][3] = 1400.0;

matmod[3][4] = 0.3;
matmod[3][5] = 0.3;
matmod[3][6] = 0.3;


matmod[4][1] = 1400.0;
matmod[4][2] = 1400.0;
matmod[4][3] = 1600.0;

matmod[4][4] = 0.35;
matmod[4][5] = 0.35;
matmod[4][6] = 0.35;

matmod[5][1] = 1600.0;
matmod[5][2] = 1600.0;
matmod[5][3] = 1800.0;

matmod[5][4] = 0.4;
matmod[5][5] = 0.4;
matmod[5][6] = 0.4;

matmod[6][1] = 1800.0;
matmod[6][2] = 1800.0;
matmod[6][3] = 2000.0;

matmod[6][4] = 0.45;
matmod[6][5] = 0.45;
matmod[6][6] = 0.45;

matmod[7][1] = 2000.0;
matmod[7][2] = 2000.0;
matmod[7][3] = 2200.0;

matmod[7][4] = 0.5;
matmod[7][5] = 0.5;
matmod[7][6] = 0.5;

matmod[8][1] = 2200.0;
matmod[8][2] = 2200.0;
matmod[8][3] = 2400.0;

matmod[8][4] = 0.55;
matmod[8][5] = 0.55;
matmod[8][6] = 0.55;

matmod[9][1] = 2400.0;
matmod[9][2] = 2400.0;
matmod[9][3] = 2600.0;

matmod[9][4] = 0.6;
matmod[9][5] = 0.6;
matmod[9][6] = 0.6;

matmod[10][1] = 2600.0;
matmod[10][2] = 2600.0;
matmod[10][3] = 2800.0;

matmod[10][4] = 0.65;
matmod[10][5] = 0.65;
matmod[10][6] = 0.65;

 
L2big=10000;
hs = 47; /* thickness of the water layer in gridpoints */

if((INVMAT==1)||(INVMAT==0)){
/*forward_prop_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
forward_prop_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_x =  vector(1,nxnyi*((NT/DTINV)));
forward_prop_y =  vector(1,nxnyi*((NT/DTINV)));

gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==2)||(INVMAT==0)){
/*forward_prop_rho_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
forward_prop_rho_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_rho_x =  vector(1,nxnyi*((NT/DTINV)));
forward_prop_rho_y =  vector(1,nxnyi*((NT/DTINV)));

gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==3)||(INVMAT==0)){
/*forward_prop_u =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_u =  vector(1,nxnyi*((NT/DTINV)));

gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

absorb_coeff=  matrix(1,NY,1,NX);
taper_coeff=  matrix(1,NY,1,NX);

/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */
bufferlef_to_rig = matrix(1,NY,1,fdo3);
bufferrig_to_lef = matrix(1,NY,1,fdo3);
buffertop_to_bot = matrix(1,NX,1,fdo3);
bufferbot_to_top = matrix(1,NX,1,fdo3);



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

 /* memory for source position definition */
srcpos1=fmatrix(1,6,1,1);

/* memory of L2 norm */
L2t = vector(1,4);
epst = vector(1,3);
picked_times = vector(1,ntr);
	
fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

		
/* Holberg coefficients for FD operators*/
hc = holbergcoeff();

MPI_Barrier(MPI_COMM_WORLD);


/* create model grids */
/*if (READMOD) readmod_elastic(prho,ppi,pu);
    else model_elastic(prho,ppi,pu);*/


/* check if the FD run will be stable and free of numerical dispersion */
/*checkfd_ssg_elastic(FP,prho,ppi,pu,hc);*/


/* calculate 2-D array for exponential damping of reflections
   at the edges of the numerical mesh */
absorb(absorb_coeff);

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
      

/* comunication initialisation for persistent communication */
comm_ini(bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);


/* Reading source positions from SOURCE_FILE */ 	
srcpos=sources(&nsrc);

snapseis=1;
snapseis1=5;
SHOTINC=1;

/* Define used minimization norm */
/* LNORM==1 L1 Norm*/
/* LNORM==2 L2 Norm*/
/* LNORM==3 Cauchy*/
/* LNORM==4 SECH*/
LNORM=2;


for(iter=1;iter<=ITERMAX;iter++){  /* fullwaveform iteration loop */	

/* check if models are numerical stable */

if(iter==1){ 
nsmod=1;
nmod=55;
}

if(iter>1){ 
nsmod=11;
nmod=55;
}


/* generate children */

/* define mutation parameters */

if(iter<=100){
dpar1 = 100.0; /* mutation parameter for vp, vs, rho */
dpar2 = 0.1; /* mutation parameter for gvp, gvs, grho */
}

if(iter>100){
dpar1 = 20.0; /* mutation parameter for vp, vs, rho */
dpar2 = 0.05; /* mutation parameter for gvp, gvs, grho */
}

/* index children hes */
/* index mom i */
/* index dad j */

if(MYID==0){

hes=nparents+1;
for(i=1;i<=nparents-1;i++){
   
   for(j=i+1;j<=nparents;j++){
       
       nrdad=1;
       nrmom=1;
       
      for(k=1;k<=nparameter-1;k++){
      
         gene=(float)drand48();
	 
	 
	 /* mutate model parameter */
	 
	 if(k<=3){
	 dpar = iround((((2.0 * dpar1)*(float)drand48()) - dpar1));}
	 
	 if(k>3){
	 dpar = iround(10000.0 * (((2.0 * dpar2)*(float)drand48()) - dpar2))/10000.0;}
	 
	 /*printf("%e \n",dpar);*/
	 
	 if(nrdad==4){
	    matmod[hes][k] = matmod[i][k] + dpar;
	 }
	 
	 if(nrmom==4){
	    matmod[hes][k] = matmod[j][k] + dpar;
	 }
	 
	 if((gene>0.5)&&(nrmom<=3)){
	    matmod[hes][k] = matmod[i][k] + dpar;
	    nrmom++;
	 }
	 
	 if((gene<=0.5)&&(nrdad<=3)){
	    matmod[hes][k] = matmod[j][k] + dpar;
	    nrdad++;
	 }
	 
      }
      
      /*printf("mom = %d \t dad = %d \n",i,j);*/
    
   hes++;
   }
}

/*printf("%d \n",hes);*/
}

MPI_Barrier(MPI_COMM_WORLD);
exchange_mod_es(matmod,ncptot,nparameter);
MPI_Barrier(MPI_COMM_WORLD);

for(is=nsmod;is<=nmod;is++){ /* start modelling */

/* check if model is stable */

if (((matmod[is][1]-matmod[is][4]*hs*DH)+matmod[is][4]*NYG*DH) > maxvp){matmod[is][7]=L2big;}
if (((matmod[is][2]-matmod[is][5]*hs*DH)+matmod[is][5]*NYG*DH) > maxvs){matmod[is][7]=L2big;}
/*if (((matmod[is][3]-matmod[is][6]*hs*DH)+matmod[is][6]*NYG*DH) > maxrho){matmod[is][7]=L2big;}*/

if (matmod[is][7] < L2big){ /* test if model is stable */


MPI_Barrier(MPI_COMM_WORLD);
readmod_elastic_es(prho,ppi,pu,matmod,is);
MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0)
   {
   time2=MPI_Wtime();
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   fprintf(FP,"\n\n\n                   DENISE ES Generation %d \t of %d \t model %d \n",iter,ITERMAX,is);
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   }

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prho, ppi, pu);
MPI_Barrier(MPI_COMM_WORLD);

av_mue(pu,puipjp,prho);
av_rho(prho,prip,prjp);


/* Open Log File for L2 norm */

if(MYID==0){

sprintf(ES_model,"models_es_gen_%d.dat",iter);
FPL2=fopen(ES_model,"w");

}

/*printf("%d \t %f \t %f \t %f \t %f \t %f \t %f \t %e \n",MYID,matmod[7][1],matmod[7][2],matmod[7][3],matmod[7][4],matmod[7][5],matmod[7][6],matmod[7][7]);*/

if(MYID==0){	

        /*for(i=1;i<=ncptot;i++){
        printf("%f \t %f \t %f \t %f \t %f \t %f \t %e \n",matmod[i][1],matmod[i][2],matmod[i][3],matmod[i][4],matmod[i][5],matmod[i][6],matmod[i][7]);
        }*/

/*fclose(FPL2);*/
}
/*printf("%d \t %f \t %f \t %f \t %f \t %f \t %f \t %e \n",MYID,matmod[5][1],matmod[5][2],matmod[5][3],matmod[5][4],matmod[5][5],matmod[5][6],matmod[5][7]);*/
/*if(iter>1){
FPL2=fopen(,"a");}
}*/

/* initialization of L2 calculation */
L2=0.0;
Lcount=0;

EPSILON=0.0;  /* test step length */
exchange_par();

if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
  
        /*for (ishot=1;ishot<=nshots;ishot+=SHOTINC){*/
	for (ishot=25;ishot<=25;ishot+=SHOTINC){

 
                fprintf(FP,"\n==================================================================================\n");
                fprintf(FP,"\n MYID=%d *****  Starting simulation (forward model) for shot %d of %d  ********** \n",MYID,ishot,nshots);
		fprintf(FP,"\n==================================================================================\n\n");
		
                for (nt=1;nt<=6;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
                if (RUN_MULTIPLE_SHOTS){

                                /* find this single source positions on subdomains */
                                if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
                                srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		           }


                       /* else{*/
                                /* Distribute multiple source positions on subdomains */
                               /* srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                        }*/

/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc);

		    
/* initialize wavefield with zero */
zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,pvxm1,pvym1,pvxp1,pvyp1);	
     
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

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2]))err(" Simulation is unstable !");*/

		
   infoout = !(nt%1000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

   update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout);
   
   /*
   update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, psxx, psyy, psxy, ppi, pu, puipjp,
                       absorb_coeff, hc);
   */

    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);


   if ((FREE_SURF) && (POS[2]==0))
      surface_elastic(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho);	


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


   /* update of particle velocities */
   update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp,
		srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,2); 


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
      

	/* store amplitudes at receivers in section-arrays */
	if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
	}



   /* WRITE SNAPSHOTS TO DISK */
      /*if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }*/
    
      
   if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (forward model) ----------*/

	
	
if ((ntr > 0) && (SEISMO)){
	saveseis(FP,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,0);
}

/*if ((ntr > 0) && (SEISMO) && (iter==snapseis) && (ishot==50)){
	saveseis(FP,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,iter);
	snapseis = snapseis + snapseis1;
}*/


/* define frequency shift */

freqshift=1.0;

if(iter>= 100){
freqshift=2.0;}

/*if(iter>= 200){
freqshift=4.0;}*/

if(iter>= 200){
FREQFILT=0;}

if (MYID==0){
printf("Calculate residuals  \n");
printf("-------------------  \n");
}

if (ntr > 0){
/* read seismic data from SU file vx */
inseis(fprec,ishot,sectionread,ntr_glob,ns,1);

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[i][j];
   }
   h++;
}

L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1);

/* apply frequency filter on the residuals (x-component) */
if(FREQFILT==1){
FFT_filt(sectionvxdiff,freqshift,ntr,ns,1);}

/* read seismic data from SU file vy */
inseis(fprec,ishot,sectionread,ntr_glob,ns,2);

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[i][j];
   }
   h++;
}

L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1);		   	    

/* apply frequency filter on the residuals (y-component) */
if(FREQFILT==1){
FFT_filt(sectionvydiff,freqshift,ntr,ns,1);}


if(TIMEWIN==1){
time_window(sectionvydiff,picked_times,iter,ntr,ns);
time_window(sectionvxdiff,picked_times,iter,ntr,ns);}

if(TIMEWIN==2){
stalta(sectionvy,ntr,ns,picked_times);
time_window(sectionvydiff,picked_times,iter,ntr,ns);

stalta(sectionvx,ntr,ns,picked_times);
time_window(sectionvxdiff,picked_times,iter,ntr,ns);}

}

/*if ((ntr > 0) && (SEISMO)){
	saveseis(FP,sectionvxdiff,sectionvydiff,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,-1);
}*/

	
nsrc_loc=0;
} /* end of loop over shots (forward modeling) */   

/* calculate L2 norm of all CPUs*/
L2sum = 0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
matmod[is][7] = L2sum;

} /* end modelling */	
} /* end of test if model is stable */	
	
        /* calculate average Vp, Vs and rho of all CPUs*/
        /*Lcountsum = 0.0;
        MPI_Allreduce(&Lcount,&Lcountsum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        Lcount=Lcountsum;*/
	
        

if(MYID==0){	

        for(i=1;i<=ncptot;i++){
        fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \n",matmod[i][1],matmod[i][2],matmod[i][3],matmod[i][4],matmod[i][5],matmod[i][6],matmod[i][7]);
        }

}

if(MYID==0){
/* find 4 best fitting models */
for(i=1;i<=nparents;i++){
   
   minL2  = matmod[1][7];
   minind = 1;
   for(j=2;j<=ncptot;j++){
   
     if(matmod[j][7] < minL2){
       minL2 = matmod[j][7];
       minind = j;
     }
   
   }

  /* copy best fitting model parameter to tmp matrix matmod1 */
  for(k=1;k<=nparameter;k++){
    matmod1[i][k] = matmod[minind][k];
  }
  
  /* set L2 = L2big at line minind in matmod */
  matmod[minind][7]=L2big;
  
}

/* replace matmod by matmod1 */

for(i=1;i<=ncptot;i++){
   
 
  for(k=1;k<=nparameter;k++){
    matmod[i][k] = 0.0;
    matmod[i][k] = matmod1[i][k];
  }
   
}

}
MPI_Barrier(MPI_COMM_WORLD);
exchange_mod_es(matmod,ncptot,nparameter);
MPI_Barrier(MPI_COMM_WORLD);

if(MYID==0){
  fclose(FPL2);
}

freqshift += dfreqshift;

/* ====================================== */
} /* end of fullwaveform iteration loop*/
/* ====================================== */

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
free_matrix(wcpart,1,3,1,3);
free_matrix(gradg,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(wavejac,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(absorb_coeff,1,NY,1,NX);
free_matrix(taper_coeff,1,NY,1,NX);

free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
free_matrix(buffertop_to_bot,1,NX,1,fdo3);
free_matrix(bufferbot_to_top,1,NX,1,fdo3);

free_vector(hc,0,6);

if((INVMAT==1)||(INVMAT==0)){
free_vector(forward_prop_x,1,NY*NX*NT);
free_vector(forward_prop_y,1,NY*NX*NT);
}

if((INVMAT==2)||(INVMAT==0)){
free_vector(forward_prop_rho_x,1,NY*NX*NT);
free_vector(forward_prop_rho_y,1,NY*NX*NT);
free_matrix(gradg_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho_s,-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==3)||(INVMAT==0)){
free_vector(forward_prop_u,1,NY*NX*NT);
free_matrix(gradg_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_mu,-nd+1,NY+nd,-nd+1,NX+nd);
}

if (nsrc_loc>0){	
	free_matrix(signals,1,nsrc_loc,1,NT);
	free_matrix(dsignals,1,nsrc_loc,1,NT);
	free_matrix(srcpos_loc,1,6,1,nsrc_loc);
	free_matrix(srcpos_loc_back,1,6,1,nsrc_loc);
}
		   

if (SEISMO) free_imatrix(recpos,1,3,1,ntr_glob);

/* free memory for global source positions */
free_matrix(srcpos,1,6,1,nsrc);

if ((ntr>0) && (SEISMO)){	

      free_imatrix(recpos_loc,1,3,1,ntr);
	switch (SEISMO){
	case 1 : /* particle velocities only */
		free_matrix(sectionvx,1,ntr,1,ns);
		free_matrix(sectionvy,1,ntr,1,ns);		
		break;	
	case 2 : /* pressure only */
		free_matrix(sectionp,1,ntr,1,ns);
		free_matrix(sectionpn,1,ntr,1,ns);
		free_matrix(sectionpnp1,1,ntr,1,ns);
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

 /* free memory for source position definition */
 free_matrix(srcpos1,1,6,1,1);
 
 free_matrix(sectionvxdata,1,ntr,1,ns);
 free_matrix(sectionread,1,ntr_glob,1,ns);
 free_matrix(sectionvxdiff,1,ntr,1,ns);
 free_matrix(sectionvydata,1,ntr,1,ns);
 free_matrix(sectionvydiff,1,ntr,1,ns);	
 free_matrix(sectionvydiffold,1,ntr,1,ns);
 free_matrix(sectionvxdiffold,1,ntr,1,ns);		
 free_vector(L2t,1,4);
 free_vector(epst,1,3);
 free_vector(picked_times,1,ntr);
 
 free_matrix(matmod,1,ncptot,1,nparameter);
 free_matrix(matmod1,1,ncptot,1,nparameter);

  
/* de-allocate buffer for messages */
MPI_Buffer_detach(buff_addr,&buffsize);
	
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
