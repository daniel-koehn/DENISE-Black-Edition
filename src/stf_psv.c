/*  --------------------------------------------------------------------------
 *   Estime source time function by Wiener deconvolution (PSV problem) 
 *   for a single shot 
 * 
 *   
 *   D. Koehn
 *   Kiel, 24.04.2016
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void stf_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV,
        struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, float *hc, int ishot, int nshots, int nsrc_loc, int nsrc, 
        float ** srcpos_loc,  float ** srcpos, int ** recpos_loc, int ** recpos, float ** signals, int ns, int ntr, int ntr_glob, int iter, float **Ws, float **Wr, int hin, int *DTINV_help, 
        int * recswitch, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
        extern int QUELLART, ORDER, ORDER_SPIKE, TIME_FILT, RUN_MULTIPLE_SHOTS;
        extern float FC_SPIKE_1, FC_SPIKE_2;
        extern float FC, FC_START;

        /* local variables */
        int nt;

        FILE *fprec;

        QUELLART = 6;

	/* calculate wavelet for each source point */
        signals=NULL;
	signals=wavelet(srcpos_loc,nsrc_loc,ishot);

	if (nsrc_loc){if(QUELLART==6){

   		/* time domain filtering of the source signal */
   		apply_tdfilt(signals,nsrc_loc,ns,ORDER_SPIKE,FC_SPIKE_2,FC_SPIKE_1);

	   }
	} 

        /* forward problem */
        psv(&wavePSV,&wavePSV_PML,&matPSV,&fwiPSV,&mpiPSV,&seisPSV,&seisPSVfwi,hc,ishot,nshots,nsrc_loc,srcpos_loc, 
	recpos_loc,signals,ns,ntr,Ws,Wr,hin,DTINV_help,0,req_send,req_rec);
	
        catseis((*seisPSV).sectionvy, (*seisPSV).fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);	   

	/* estimate STF */	
	   if (nsrc_loc>0){
	       
	      /* read seismic data from SU file vy */
	      /* --------------------------------- */
	      inseis(fprec,ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,2,iter);

	      if (TIME_FILT){
                 apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
	      }

	      stf((*seisPSVfwi).sectionread,(*seisPSV).fulldata_vy,ntr_glob,ishot,ns,iter,nshots,signals,recpos,srcpos);
	      /*saveseis_glob(FP,seisPSVfwi.sectionread,seisPSV.fulldata_vy,seisPSV.sectionp,seisPSV.sectioncurl,seisPSV.sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);*/
	      
	   }

	   MPI_Barrier(MPI_COMM_WORLD);

	   QUELLART=3;

       
}
