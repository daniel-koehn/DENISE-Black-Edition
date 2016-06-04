/*------------------------------------------------------------------------
 *   globvar.h - global variables of viscoelastic FD programs   
 *   last update 05/04/02 T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

/* definition of global variables used in the finite difference programs*/
/* generally, for the names of the global variables
   uppercase letters are used */

float XS, YS, DH, TIME, DT, TS, DAMPING, PLANE_WAVE_DEPTH, PHI;
float TSNAP1, TSNAP2, TSNAPINC, *FL, TAU;
float XREC1, XREC2, YREC1, YREC2, ANGLE;
float REFREC[4]={0.0, 0.0, 0.0, 0.0}, FPML;
int   SEISMO, NDT, NGEOPH, NSRC=1, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC, FW=0;
int   NX, NY, NT, QUELLART, QUELLTYP, SNAP, SNAP_FORMAT, LOG, RUN_MULTIPLE_SHOTS, NTRG;
int   L, BOUNDARY, DC, NXG, NYG, IDX, IDY, CHECKPTREAD, CHECKPTWRITE, FDORDER, MAXRELERROR;
char  SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
char  MFILE[STRING_SIZE], REC_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE], LOG_FILE[STRING_SIZE];
char  SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
char  *FILEINP1; 
FILE *FP;

/* Mpi-variables */
int   NP, NPSP, NPROC, NPROCX, NPROCY, MYID, IENDX, IENDY;
int   POS[3], INDEX[5];     
const int TAG1=1,TAG2=2, TAG3=3, TAG4=4, TAG5=5,TAG6=6; 

/* spatial adaptive Code variables*/
int   check_id, cfgt_id, cfgt, jumpid;
float DH1;

/* TDFWI Code DENISE_elastic*/
char JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE];
int TAPER, TAPERLENGTH, INVTYPE;
int GRADT1,GRADT2,GRADT3,GRADT4;
int ITERMAX, REC1, REC2, INVMAT1, INVMAT, QUELLTYPB;
int HESSIAN, GRAD_METHOD, NLBFGS;
float FC_HESS_START, FC_HESS_INC;
int MODEL_FILTER, FILT_SIZE;
float EPSILON, MUN, EPSILON_u, EPSILON_rho;

int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
int SWS_TAPER_FILE;
float SRTRADIUS, EXP_TAPER_GRAD_HOR;

int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
int INV_RHO_ITER, INV_VS_ITER, INV_VP_ITER;
int MIN_ITER, GRAD_FORM, IDXI, IDYI, NTDTINV, NXNYI;

char INV_MODELFILE[STRING_SIZE];

float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;

float npower, k_max_PML;

int INV_STF, N_STF, N_STF_START;
float EPS_STF, OFFSETC_STF;
char PARA[STRING_SIZE];

int TIME_FILT, ORDER, SHOTINC;
float FC_START, FC_END, FC_INCR;

int LNORM, DTINV;

int STEPMAX;
float EPS_SCALE, SCALEFAC;

float PRO;

/* trace kill variables */
int TRKILL;
char TRKILL_FILE[STRING_SIZE];

int TIMEWIN, NORMALIZE;
float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
char PICKS_FILE[STRING_SIZE];
char MISFIT_LOG_FILE[STRING_SIZE]; 

int GRAD_FILTER, FILT_SIZE_GRAD;

/* parameter for RTM */
int RTM;

/* parameter for double difference time-lapse mode */
int TIMELAPSE;
char DATA_DIR_T0[STRING_SIZE];

/* parameters for energy weighted preconditioning */
int EPRECOND;

/* parameters for FD Hessian calculation */
int NFREQ;
float FC_HESS_START, FC_HESS_INC, FC;

/* parameters for wavenumber domain  damping */
float WD_DAMP, WD_DAMP1;

/* Reverse time modelling parameter */
int RTMOD;

/* parameters for offset muting */
float OFFSETC;
int OFFSET_MUTE;

/* scale density update */
float SCALERHO;

/* corner frequencies for spike wavelet */
float FC_SPIKE_1, FC_SPIKE_2; 
int ORDER_SPIKE;

/* parameter for gravity modelling/inversion */
int GRAVITY, NZGRAV, NGRAVB, GRAV_TYPE, BACK_DENSITY;
char GRAV_DATA_OUT[STRING_SIZE], GRAV_DATA_IN[STRING_SIZE], GRAV_STAT_POS[STRING_SIZE], DFILE[STRING_SIZE];
float LAM_GRAV, GAMMA_GRAV, LAM_GRAV_GRAD, L2_GRAV_IT1;

/* parameters for envelope calculation */
int ENV;

/* parameters for towed streamer */
int N_STREAMER;
float REC_INCR_X, REC_INCR_Y;

int PHYSICS;
