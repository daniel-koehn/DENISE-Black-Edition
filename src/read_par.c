/**/
/*------------------------------------------------------------------------
 *   Read FD-Parameters from Stdin                           
 *   last update 05/04/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

/* reading input-parameter from input-file or stdin for
viscoelastic finite-differnce modelling with fdveps */

#include "fd.h"

void read_par(FILE *fp_in){

/* declaration of extern variables */
extern int   NX, NY, FDORDER, MAXRELERROR, QUELLART, SNAP, SNAP_FORMAT, SNAP_SHOT, L;
extern float DH, TIME, DT, TS, *FL, TAU, DAMPING;
extern float FPML;
extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC, RUN_MULTIPLE_SHOTS;
extern int BOUNDARY, LOG, TAPER, TAPERLENGTH, FW, PHYSICS, MODE;
extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE], TFILE[STRING_SIZE];
extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
extern char JACOBIAN[STRING_SIZE],DATA_DIR[STRING_SIZE];
extern int  NPROC, NPROCX, NPROCY, MYID, IDX, IDY; 
extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT1, GRAD_FORM, QUELLTYPB;
extern int HESSIAN, GRAD_METHOD, NFREQ, PCG_BETA, RTM_SHOT;
extern float FC_HESS_START, FC_HESS_INC;
extern int FILT_SIZE, MODEL_FILTER;
extern int FILT_SIZE_GRAD, GRAD_FILTER;

extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
extern int SWS_TAPER_FILE;
extern float SRTRADIUS, EXP_TAPER_GRAD_HOR;
extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
extern int MIN_ITER, RTMOD, INV_MOD_OUT;
extern char INV_MODELFILE[STRING_SIZE];
extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
extern float QSUPPERLIM, QSLOWERLIM;
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

extern float WD_DAMP;

extern float FC_SPIKE_1, FC_SPIKE_2;
extern int ORDER_SPIKE;

extern int N_STREAMER;
extern float REC_INCR_X, REC_INCR_Y;

extern int NLBFGS;

extern int GRAVITY, NGRAVB, NZGRAV, GRAV_TYPE;
extern char DFILE[STRING_SIZE];
extern int BACK_DENSITY;

extern int WRITE_STF, WRITEMOD;

/* definition of local variables */
char s[120];
int  c=0, lineno=0, l;
 
   while ((c=getc(fp_in)) != EOF){
      if ((c=='\n') && (getc(fp_in)!='#')){     
	 lineno++;
	 /* printf(" reading line %d \n",lineno);*/
	 switch (lineno){
         case 1 :
	    fscanf(fp_in,"%s =%i",s,&MODE);
	    break;
         case 2 :
	    fscanf(fp_in,"%s =%i",s,&PHYSICS);
	    break;
	 case 3 :
	    fscanf(fp_in,"%s =%i",s,&NPROCX);
	    break;
	 case 4 :
	    fscanf(fp_in,"%s =%i",s,&NPROCY);
	    break;     		
	 case 5 :
	    fscanf(fp_in,"%s =%i",s,&FDORDER);
	    break;
	 case 6 :
	    fscanf(fp_in,"%s =%i",s,&MAXRELERROR);
	    break;
	 case 7 :
	    fscanf(fp_in,"%s =%i",s,&NX);
	    break;
	 case 8 :
	    fscanf(fp_in,"%s =%i",s,&NY);
	    break;
	 case 9 :
	    fscanf(fp_in,"%s =%f",s,&DH);
	    break;
	 case 10 :
	    fscanf(fp_in,"%s =%f",s,&TIME);
	    break;
	 case 11 :
	    fscanf(fp_in,"%s =%f",s,&DT);
	    break;
	 case 12 :
	    fscanf(fp_in,"%s =%i",s,&QUELLART);
	    break;
	 case 13 :
	    fscanf(fp_in,"%s =%s",s,SIGNAL_FILE);
	    break;
	 case 14 :
	    fscanf(fp_in,"%s =%f",s,&TS);
	    break;
	 case 15 :
	    fscanf(fp_in,"%s =%i",s,&SRCREC);
	    break;
	 case 16 :
	    fscanf(fp_in,"%s =%s",s,SOURCE_FILE);
	    break;
	 case 17 :
            fscanf(fp_in,"\n%s =%i",s,&RUN_MULTIPLE_SHOTS);
            break;
	 case 18 :
            fscanf(fp_in,"%s =%f",s,&FC_SPIKE_1);
            break;
         case 19 :
            fscanf(fp_in,"%s =%f",s,&FC_SPIKE_2);
            break;   
	 case 20 :
            fscanf(fp_in,"%s =%i",s,&ORDER_SPIKE);
            break;
	 case 21 :
            fscanf(fp_in,"%s =%i",s,&WRITE_STF);
            break;
	 case 22 :
	    fscanf(fp_in,"%s =%i",s,&READMOD);
	    break;
	 case 23 :
	    fscanf(fp_in,"%s =%s",s,MFILE);
	    break;
	 case 24 :
            fscanf(fp_in,"%s =%i",s,&WRITEMOD);
            break;
	 case 25 :
	    fscanf(fp_in,"%s =%i",s,&L);
	    FL=vector(1,L);
	    break;
	 case 26 :
	    fscanf(fp_in,"%s =%f",s,&FL[1]);
	    for (l=2;l<=L;l++) fscanf(fp_in,"%f",&FL[l]);
	    break;
	 case 27 :
	    fscanf(fp_in,"%s =%f",s,&TAU);
	    break;
	 case 28 :
	    fscanf(fp_in,"%s =%i",s,&FREE_SURF);
	    break;
	 case 29 :
	    fscanf(fp_in,"%s =%i",s,&FW);
	    if (FW<0) FW=0;
	    break;
	 case 30 :
	    fscanf(fp_in,"%s =%f",s,&DAMPING);
	    break;
	 case 31 :
            fscanf(fp_in,"%s =%f",s,&FPML);
            break;
	 case 32 :
            fscanf(fp_in,"%s =%f",s,&npower);
            break;
	 case 33 :
            fscanf(fp_in,"%s =%f",s,&k_max_PML);
            break;     			
	 case 34 :
	    fscanf(fp_in,"%s =%i",s,&BOUNDARY);
	    break;			
	 case 35 :
	    fscanf(fp_in,"%s =%i",s,&SNAP);
	    break;
	 case 36 :
	    fscanf(fp_in,"%s =%i",s,&SNAP_SHOT);
	    break;
	 case 37 :
	    fscanf(fp_in,"%s =%f",s,&TSNAP1);
	    break;
	 case 38 :
	    fscanf(fp_in,"%s =%f",s,&TSNAP2);
	    break;
	 case 39 :
	    fscanf(fp_in,"%s =%f",s,&TSNAPINC);
	    break;
	 case 40 :
	    fscanf(fp_in,"%s =%i",s,&IDX);
	    break;
	 case 41 :
	    fscanf(fp_in,"%s =%i",s,&IDY);
	    break;
	 case 42 :
	    fscanf(fp_in,"%s =%i",s,&SNAP_FORMAT);
	    break;
	 case 43 :
	    fscanf(fp_in,"%s =%s",s,SNAP_FILE);
	    break;
	 case 44 :
	    fscanf(fp_in,"%s =%i",s,&SEISMO);
	    break;
	 case 45 :
	    fscanf(fp_in,"%s =%i",s,&READREC);
	    break;
	 case 46 :
	    fscanf(fp_in,"%s =%s",s,REC_FILE);
	    break;
	 case 47 :
	    fscanf(fp_in,"%s =%f ,%f",s,&REFREC[1],&REFREC[2]);
	    break;
	 case 48 :
	    fscanf(fp_in,"%s =%i",s,&N_STREAMER);
	    break;
	 case 49 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_X);
	    break;
	 case 50 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_Y);
	    break;
	 case 51 :
	    fscanf(fp_in,"%s =%i",s,&NDT);
	    break;
	 case 52 :
	    fscanf(fp_in,"%s =%i",s,&SEIS_FORMAT);
	    break;
	 case 53 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_VX);
	    break;
	 case 54 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_VY);
	    break;
	 case 55 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_CURL);
	    break;
	 case 56 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_DIV);
	    break;
	 case 57 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_P);
	    break;     			
	 case 58 :
	    fscanf(fp_in,"%s =%s",s,LOG_FILE);
	    break;     			
	 case 59 :
	    fscanf(fp_in,"%s =%i",s,&LOG);
	    break; 
	 case 60 :
	    fscanf(fp_in,"%s =%i",s,&ITERMAX);
	    break;    
	 case 61 :
	    fscanf(fp_in,"%s =%s",s,JACOBIAN);
	    break;   
	 case 62 :
	    fscanf(fp_in,"%s =%s",s,DATA_DIR);
	    break;    
	 case 63 :
	    fscanf(fp_in,"%s =%i",s,&TAPER);
	    break;
	 case 64 :
	    fscanf(fp_in,"%s =%i",s,&TAPERLENGTH);
	    break;        	
	 case 65 :
	    fscanf(fp_in,"%s =%i, %i, %i, %i",s,&GRADT1,&GRADT2,&GRADT3,&GRADT4);
	    break; 
	 case 66 :
	    fscanf(fp_in,"%s =%i",s,&INVMAT1);
	    break;
	 case 67 :
	    fscanf(fp_in,"%s =%i",s,&GRAD_FORM);
	    break;
         case 68 :
            fscanf(fp_in,"%s =%i",s,&QUELLTYPB);
            break; 
	 case 69 :
	    fscanf(fp_in,"%s =%i, %i, %i",s,&TESTSHOT_START,&TESTSHOT_END,&TESTSHOT_INCR); 
	    break; 
	 case 70 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_VERT);
            break;            	        
	 case 71 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_HOR);
            break; 
         case 72 :
            fscanf(fp_in,"%s =%f",s,&EXP_TAPER_GRAD_HOR);
            break;
	 case 73 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_SOURCES);
            break; 
	 case 74 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_CIRCULAR_PER_SHOT);
            break;    
	 case 75 :
            fscanf(fp_in,"%s =%i",s,&SRTSHAPE);
            break;   
	 case 76 :
            fscanf(fp_in,"%s =%f",s,&SRTRADIUS);
            break; 
	 case 77 :
            fscanf(fp_in,"%s =%i",s,&FILTSIZE);
            break;
         case 78 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_FILE);
            break;
	 case 79 :
	    fscanf(fp_in,"%s =%s",s,TFILE);
	    break;
	 case 80 :
            fscanf(fp_in,"%s =%i",s,&INV_MOD_OUT);
            break;                          
	 case 81 :
            fscanf(fp_in,"%s =%s",s,INV_MODELFILE);
            break;
	 case 82 :
            fscanf(fp_in,"%s =%f",s,&VPUPPERLIM);
            break;
	 case 83 :
            fscanf(fp_in,"%s =%f",s,&VPLOWERLIM);
            break; 
	 case 84 :
            fscanf(fp_in,"%s =%f",s,&VSUPPERLIM);
            break; 
	 case 85 :
            fscanf(fp_in,"%s =%f",s,&VSLOWERLIM);
            break;
	 case 86 :
            fscanf(fp_in,"%s =%f",s,&RHOUPPERLIM);
            break;
	 case 87 :
            fscanf(fp_in,"%s =%f",s,&RHOLOWERLIM);
            break;
	 case 88 :
            fscanf(fp_in,"%s =%f",s,&QSUPPERLIM);
            break;
	 case 89 :
            fscanf(fp_in,"%s =%f",s,&QSLOWERLIM);
            break;
	 case 90 :         
            fscanf(fp_in,"%s =%i",s,&GRAD_METHOD);                         
            break;
	 case 91 :         
            fscanf(fp_in,"%s =%i",s,&PCG_BETA);                         
            break;
	 case 92 :         
            fscanf(fp_in,"%s =%i",s,&NLBFGS);                         
            break;   
	 case 93 :         
            fscanf(fp_in,"%s =%i",s,&MODEL_FILTER);                         
            break;  
	 case 94 :         
            fscanf(fp_in,"%s =%i",s,&FILT_SIZE);                         
            break;  
	 case 95 :
	   fscanf(fp_in,"%s =%i",s,&DTINV);                         
            break;
	 case 96 :
	   fscanf(fp_in,"%s =%f",s,&EPS_SCALE);                         
            break;
	 case 97 :
	   fscanf(fp_in,"%s =%i",s,&STEPMAX);                         
            break;
	 case 98 :
	   fscanf(fp_in,"%s =%f",s,&SCALEFAC);                         
            break;
	 case 99 :
	   fscanf(fp_in,"%s =%i",s,&TRKILL);                         
            break;
	 case 100 :
	   fscanf(fp_in,"%s =%s",s,TRKILL_FILE);                         
            break;
	 case 101 :
	   fscanf(fp_in,"%s =%s",s,PICKS_FILE);                         
            break;
	 case 102 :
	   fscanf(fp_in,"%s =%s",s,&MISFIT_LOG_FILE);                         
            break; 
	 case 103 :
	   fscanf(fp_in,"%s =%i",s,&MIN_ITER);                         
            break; 
	 case 104 :
	   fscanf(fp_in,"%s =%i",s,&GRAD_FILTER);                         
            break; 
	 case 105 :
	   fscanf(fp_in,"%s =%i",s,&FILT_SIZE_GRAD);                         
            break;
         case 106 : 
           fscanf(fp_in,"%s =%i",s,&TIMELAPSE);
            break;
         case 107 :                      
           fscanf(fp_in,"%s =%s",s,DATA_DIR_T0);
            break;       
         case 108 :
           fscanf(fp_in,"%s =%i",s,&RTMOD);
           break;
	 case 109 :
	   fscanf(fp_in,"%s =%i",s,&GRAVITY);
	    break;
	 case 110 :
	   fscanf(fp_in,"%s =%i",s,&NGRAVB);
	    break;
	 case 111 :
	   fscanf(fp_in,"%s =%i",s,&NZGRAV);
	    break;
	 case 112 :
	   fscanf(fp_in,"%s =%i",s,&GRAV_TYPE);
	    break;   
	 case 113 :
	   fscanf(fp_in,"%s =%i",s,&BACK_DENSITY);
	    break;
	 case 114 :
	   fscanf(fp_in,"%s =%s",s,DFILE);
	    break;
	 case 115 :
	   fscanf(fp_in,"%s =%i",s,&RTM_SHOT);
	    break;   
	 default:
	    break;
	 }
	 }
      }

/* define total number of MPI processes*/
NPROC = NPROCX*NPROCY;

fclose(fp_in);

}
