/**/
/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 *   last update 29/06/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){

	/* declaration of extern variables */
	extern int   NX, NY, FDORDER, MAXRELERROR, QUELLART, QUELLTYP, SNAP, SNAP_FORMAT, L;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, PLANE_WAVE_DEPTH, PHI;
	extern float XREC1, XREC2, YREC1, YREC2, FPML;
	extern float MUN, EPSILON, EPSILON_u, EPSILON_rho;
	extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC;
	extern int BOUNDARY, LOG, FW;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4], ANGLE;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern char JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE], INV_MODELFILE[STRING_SIZE];
	extern int RUN_MULTIPLE_SHOTS, TAPERLENGTH, INVTYPE;
	extern int NPROC, NPROCX, NPROCY, MYID, IDX, IDY, CHECKPTREAD, CHECKPTWRITE;
	extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT1, INVMAT, QUELLTYPB;
	extern int HESSIAN, GRAD_METHOD, NFREQ;
	extern float FC_HESS_START, FC_HESS_INC; 
	extern int MODEL_FILTER, FILT_SIZE;
	extern int FILT_SIZE_GRAD, GRAD_FILTER;
	
	extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
	extern int SWS_TAPER_FILE;
	extern float SRTRADIUS, WD_DAMP, EXP_TAPER_GRAD_HOR;
	extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
	extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
	extern int MIN_ITER;
	extern int nfstart, nf;
	extern int nfstart_jac, nf_jac;
	extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
	extern float npower, k_max_PML;
	extern int INV_STF, N_STF, N_STF_START;
	extern char PARA[STRING_SIZE];
	extern int TIME_FILT, ORDER;
	extern float FC_START, FC_END, FC_INCR;
	extern int LNORM, DTINV;
	extern int STEPMAX;
	extern float EPS_SCALE, SCALEFAC;
	extern float PRO;
	extern int TRKILL;
	extern char TRKILL_FILE[STRING_SIZE];
	extern int TIMEWIN, NORMALIZE;
	extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
	extern char PICKS_FILE[STRING_SIZE];
	extern char MISFIT_LOG_FILE[STRING_SIZE];
        extern int TIMELAPSE;
        extern char DATA_DIR_T0[STRING_SIZE]; 
        extern float FC_SPIKE_1, FC_SPIKE_2;
        extern int ORDER_SPIKE;
        extern int RTM, NLBFGS;
        extern int N_STREAMER, RTMOD;
        extern float REC_INCR_X, REC_INCR_Y;
        
	/* definition of local variables */
	int idum[NPAR];
	float fdum[NPAR];
	
	
	if (MYID == 0){ 

	fdum[1]  = DH;                                                                                 
        fdum[2]  = TIME;                                                                               
        fdum[3]  = DT;                                                                                 
        fdum[4]  = TS;                                                                                                                                                                  

        fdum[5]  = TAU;                                                                                                                                                                 
        fdum[6]  = TSNAP1;                                                                            
        fdum[7]  = TSNAP2;                                                                            
        fdum[8]  = TSNAPINC;                                                                          
        fdum[9]  = REFREC[1];                                                                         
        fdum[10]  = REFREC[2];                                                                         
        fdum[11]  = PHI;                                                                             

        fdum[12]  = XREC1;                                                                             
        fdum[13]  = YREC1;                                                                             

        fdum[14]  = XREC2;                                                                             
        fdum[15]  = YREC2;                                                                             

        fdum[16]  = DAMPING;                                                                           
        fdum[17]  = REC_INCR_X;                                                                   
        fdum[18]  = REC_INCR_Y;                                                                    
	
	fdum[19]  = MUN;
	fdum[20]  = EPSILON;
	fdum[21]  = EPSILON_u;
	fdum[22]  = EPSILON_rho;
	fdum[23]  = ANGLE;
	fdum[24]  = FPML;
	
	fdum[25]  = SRTRADIUS;
	
	fdum[26]  = VPUPPERLIM;	
	fdum[27]  = VPLOWERLIM;	
	fdum[28]  = VSUPPERLIM;	
	fdum[29]  = VSLOWERLIM;	
	fdum[30]  = RHOUPPERLIM;	
	fdum[31]  = RHOLOWERLIM;
	
	fdum[32]  = npower;
	fdum[33]  = k_max_PML;
	
	fdum[34]  = FC_SPIKE_1;
	
	fdum[35]  = FC_START;
	fdum[36]  = FC_END;
	fdum[37]  = FC_INCR;
	
	fdum[38]  = EPS_SCALE;
	fdum[39]  = SCALEFAC;
	fdum[40]  = PRO;
	fdum[41]  = WD_DAMP;
	fdum[42]  = FC_HESS_START;
	fdum[43]  = FC_HESS_INC;
	
	fdum[44] = FC_SPIKE_2;	
        fdum[45] = EXP_TAPER_GRAD_HOR;
	                                                                                                                                                                                                                                                             
        idum[1]  = NPROCX;                                                                             
        idum[2]  = NPROCY;                                                                             
        idum[3]  = LOG;                                                                             

        idum[4]  = NPROCX*NPROCY;                                                               
        idum[5]  = NX;                                                                                 
        idum[6]  = NY;                                                                                 

        idum[8]  = QUELLART;                                                                           
        idum[9]  = QUELLTYP;                                                                           
        idum[10]  = READMOD;                                                                           
        idum[11]  = L;                                                                                 
        idum[12]  = FREE_SURF;                                                                         
        idum[13]  = SNAP;                                                                                                                                                             

        idum[16]  = BOUNDARY;                                                                                                                                                   
        idum[18]  = SRCREC;                                                                                 
        idum[19]  = IDX;                                                                                 
        idum[20]  = IDY;                                                                                 
        idum[21]  = 0;                                                                                 
        idum[22]  = 0;                                                                                 
        idum[23]  = SNAP_FORMAT;                                                                       
        idum[24]  = SEISMO;                                                                            
        idum[25]  = READREC;                                                                           
        idum[26]  = NGEOPH;                                                                            
        idum[27]  = NDT;                                                                               
        idum[28]  = SEIS_FORMAT;                                                                       
        idum[29]  = CHECKPTREAD;                                                                       
        idum[30]  = CHECKPTWRITE;                                                                       
                                                                                                      
	idum[31]  = FDORDER;
	idum[32]  = MAXRELERROR;
	idum[33]  = RUN_MULTIPLE_SHOTS;	
        idum[34]  = TAPERLENGTH;
	idum[35]  = INVTYPE;
	idum[36]  = GRADT1;
	idum[37]  = GRADT2;
	idum[38]  = GRADT3;
	idum[39]  = GRADT4;
	idum[40]  = ITERMAX;
	idum[41]  = INVMAT1;
	idum[42]  = FW;
	idum[43]  = INVMAT;
	idum[44]  = QUELLTYPB;
	
	idum[45]  = TESTSHOT_START;
	idum[46]  = TESTSHOT_END;
	idum[47]  = TESTSHOT_INCR;
	
	idum[48]  = SWS_TAPER_GRAD_VERT;
	idum[49]  = SWS_TAPER_GRAD_HOR;
	idum[50]  = SWS_TAPER_GRAD_SOURCES;
	idum[51]  = SWS_TAPER_CIRCULAR_PER_SHOT;
	idum[52]  = SRTSHAPE;
	idum[53]  = FILTSIZE;
	
	idum[54]  = SPATFILTER;
	idum[55]  = SPAT_FILT_SIZE;
	idum[56]  = SPAT_FILT_1;
	idum[57]  = SPAT_FILT_ITER;
	
	idum[58]  = INV_RHO_ITER;
	idum[59]  = nfstart;
	idum[60]  = nf;
	
	idum[61]  = nfstart_jac;
	idum[62]  = nf_jac;
	idum[63]  = SWS_TAPER_FILE;
	idum[64]  = HESSIAN;
	idum[65]  = GRAD_METHOD;
	
	idum[66]  = MODEL_FILTER;
	idum[67]  = FILT_SIZE;

	idum[68]  = ORDER_SPIKE;
	
	idum[69]  = INV_STF;
	idum[70]  = N_STF;
	idum[71]  = N_STF_START;
	
	idum[72]  = TIME_FILT;
	idum[73]  = ORDER;
	
	idum[74]  = LNORM;
	idum[75]  = DTINV;
	
	idum[76]  = STEPMAX;
	
	idum[77]  = TRKILL;

	idum[78]  = TIMEWIN;
	fdum[79]  = TWLENGTH_PLUS;
	fdum[80]  = TWLENGTH_MINUS;
	fdum[81]  = GAMMA;

	idum[82]  = NORMALIZE;
	
	idum[83]  = INV_VP_ITER;
	idum[84]  = INV_VS_ITER;
	
	idum[85]  = MIN_ITER;
	
	idum[86]  = GRAD_FILTER;
	idum[87]  = FILT_SIZE_GRAD;
	idum[88]  = TIMELAPSE;
	idum[89]  = NFREQ;
	idum[90]  = RTM;
        idum[91]  = NLBFGS;
        idum[92]  = N_STREAMER;
        idum[93]  = RTMOD;
	
	} /** if (MYID == 0) **/
	
	if (MYID != 0) FL=vector(1,L);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_VX,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_VY,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_CURL,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_DIV,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_P,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&CHECKPTFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&JACOBIAN,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&DATA_DIR,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&INV_MODELFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&PARA,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MISFIT_LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD); 
	MPI_Bcast(&TRKILL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&PICKS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&DATA_DIR_T0,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);

	DH=fdum[1];
	TIME=fdum[2];
	DT=fdum[3];
	TS=fdum[4];

	TAU=fdum[5];
	TSNAP1=fdum[6];
	TSNAP2=fdum[7];
	TSNAPINC=fdum[8];
	REFREC[1]=fdum[9];
	REFREC[2]=fdum[10];
	PHI=fdum[11];
	XREC1=fdum[12];
	YREC1=fdum[13];

	XREC2=fdum[14];
	YREC2=fdum[15];

	DAMPING=fdum[16];
	REC_INCR_X=fdum[17];
	REC_INCR_Y=fdum[18];

        MUN = fdum[19];
	EPSILON = fdum[20];
	EPSILON_u = fdum[21];
	EPSILON_rho = fdum[22];
	ANGLE = fdum[23];
        FPML = fdum[24];
	
	SRTRADIUS = fdum[25];
	
	VPUPPERLIM = fdum[26];	
	VPLOWERLIM = fdum[27];	
	VSUPPERLIM = fdum[28];	
	VSLOWERLIM = fdum[29];	
	RHOUPPERLIM = fdum[30];	
	RHOLOWERLIM = fdum[31];
	
	npower = fdum[32];
	k_max_PML = fdum[33];
	
	FC_SPIKE_1 = fdum[34];
	
	FC_START = fdum[35];
	FC_END = fdum[36];
	FC_INCR = fdum[37];
	
	EPS_SCALE = fdum[38];
	SCALEFAC = fdum[39];
	
	PRO = fdum[40];
	WD_DAMP = fdum[41];
	FC_HESS_START = fdum[42];
	FC_HESS_INC = fdum[43];	
	
	FC_SPIKE_2 = fdum[44];
        EXP_TAPER_GRAD_HOR = fdum[45];

	NPROCX = idum[1];
	NPROCY = idum[2];
	LOG=idum[3];
	NPROC  = idum[4];
	NX = idum[5];
	NY = idum[6];

	QUELLART = idum[8];
	QUELLTYP = idum[9];
	READMOD = idum[10];
	L = idum[11];
	FREE_SURF = idum[12];
	SNAP = idum[13];

	BOUNDARY = idum[16];
	SRCREC = idum[18];
	IDX = idum[19];
	IDY = idum[20];
	
	
	SNAP_FORMAT = idum[23];
	SEISMO = idum[24];
	READREC = idum[25];
	NGEOPH = idum[26];
	NDT = idum[27];
	SEIS_FORMAT = idum[28];
	CHECKPTREAD = idum[29];
	CHECKPTWRITE = idum[30];

	FDORDER = idum[31];
	MAXRELERROR = idum[32];
        RUN_MULTIPLE_SHOTS = idum[33];
        TAPERLENGTH = idum[34];
	INVTYPE = idum[35]; 
	GRADT1 = idum[36];
	GRADT2 = idum[37];
	GRADT3 = idum[38];
	GRADT4 = idum[39];
	ITERMAX = idum[40];
	INVMAT1 = idum[41];
	     FW = idum[42];
	INVMAT  = idum[43];  
	QUELLTYPB = idum[44];
	
	TESTSHOT_START = idum[45];
	TESTSHOT_END = idum[46];
	TESTSHOT_INCR = idum[47];
	
	SWS_TAPER_GRAD_VERT = idum[48];
	SWS_TAPER_GRAD_HOR = idum[49];
	SWS_TAPER_GRAD_SOURCES = idum[50];
	SWS_TAPER_CIRCULAR_PER_SHOT = idum[51];
	SRTSHAPE = idum[52];
	FILTSIZE = idum[53];
	
	SPATFILTER = idum[54];
	SPAT_FILT_SIZE = idum[55];
	SPAT_FILT_1 = idum[56];
	SPAT_FILT_ITER = idum[57];
	
	INV_RHO_ITER = idum[58];
	nfstart = idum[59];
	nf = idum[60];
	
	nfstart_jac = idum[61];
	nf_jac = idum[62];
	SWS_TAPER_FILE = idum[63];
	HESSIAN = idum[64];
	GRAD_METHOD = idum[65];
	
	MODEL_FILTER = idum[66];
	FILT_SIZE = idum[67];

	ORDER_SPIKE = idum[68];
	
	INV_STF = idum[69];
	N_STF = idum[70];
	N_STF_START = idum[71];

	TIME_FILT = idum[72];
	ORDER = idum[73];
	
	LNORM = idum[74];
	DTINV = idum[75];
	
	STEPMAX = idum[76];
	
	TRKILL = idum[77];

	TIMEWIN = idum[78];
	TWLENGTH_PLUS = fdum[79];
	TWLENGTH_MINUS = fdum[80];
	GAMMA = fdum[81];

	NORMALIZE = idum[82];
	
	INV_VP_ITER = idum[83];
	INV_VS_ITER = idum[84];

	MIN_ITER = idum[85];
	
	GRAD_FILTER = idum[86];
	FILT_SIZE_GRAD = idum[87];
	TIMELAPSE = idum[88];
	NFREQ = idum[89];
	RTM = idum[90];

        NLBFGS = idum[91];
        N_STREAMER = idum[92];
        RTMOD = idum[93];	
  
	MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);

}
