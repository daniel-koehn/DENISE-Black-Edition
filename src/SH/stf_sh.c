/*  --------------------------------------------------------------------------
 *   Estime source time function by Wiener deconvolution (SH problem) 
 *   for a single shot 
 * 
 *   
 *   D. Koehn
 *   Kiel, 13.12.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void stf_sh(struct waveSH *waveSH, struct waveSH_PML *waveSH_PML, struct matSH *matSH, struct fwiSH *fwiSH, struct mpiPSV *mpiPSV, struct seisSH *seisSH, 
             struct seisSHfwi *seisSHfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, int nsrc, int ns, int ntr, int ntr_glob, int iter, float **Ws, 
             float **Wr, int hin, int *DTINV_help, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
        extern int QUELLART, ORDER, ORDER_SPIKE, TIME_FILT, RUN_MULTIPLE_SHOTS, L;
        extern float FC_SPIKE_1, FC_SPIKE_2;
        extern float FC, FC_START;

        /* local variables */
        int nt;

        FILE *fprec;

        QUELLART = 6;

	/* calculate wavelet for each source point */
        (*acq).signals=NULL;
	(*acq).signals=wavelet((*acq).srcpos_loc,nsrc_loc,ishot);

	if (nsrc_loc){
	   if(QUELLART==6){

   		/* time domain filtering of the source signal with bandlimited spike filter parameters */
   		apply_tdfilt((*acq).signals,nsrc_loc,ns,ORDER_SPIKE,FC_SPIKE_2,FC_SPIKE_1);

	   }
	   
	   if(TIME_FILT){
	   
	        /* time domain filtering of the source signal with FWI filter parameters of current stage */
   		apply_tdfilt((*acq).signals,nsrc_loc,ns,ORDER,FC,FC_START);
	   
	   }
	   
	}
	
	 
        /* forward problem */
	if(L){
           sh_visc(waveSH,waveSH_PML,matSH,fwiSH,mpiPSV,seisSH,seisSHfwi,acq,hc,ishot,nshots,nsrc_loc,ns,ntr,Ws,Wr,hin,DTINV_help,0,req_send,req_rec);
	}else{
           sh(waveSH,waveSH_PML,matSH,fwiSH,mpiPSV,seisSH,seisSHfwi,acq,hc,ishot,nshots,nsrc_loc,ns,ntr,Ws,Wr,hin,DTINV_help,0,req_send,req_rec);
	}

        catseis((*seisSH).sectionvz, (*seisSH).fulldata_vz, (*acq).recswitch, ntr_glob, MPI_COMM_WORLD);	   

	/* estimate STF */	
	   if (nsrc_loc>0){
	       
	      /* read seismic data from SU file vz */
	      /* --------------------------------- */
	      inseis(ishot,(*seisSHfwi).sectionread,ntr_glob,ns,2,iter);

	      if (TIME_FILT){
                 apply_tdfilt((*seisSHfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
	      }

	      stf((*seisSHfwi).sectionread,(*seisSH).fulldata_vz,ntr_glob,ishot,ns,iter,nshots,(*acq).signals,(*acq).recpos,(*acq).srcpos);
	      /*saveseis_glob(FP,seisPSVfwi.sectionread,seisPSV.fulldata_vy,seisPSV.sectionp,seisPSV.sectioncurl,seisPSV.sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);*/
	      
	   }

	   MPI_Barrier(MPI_COMM_WORLD);

	   QUELLART=3;
      
}
