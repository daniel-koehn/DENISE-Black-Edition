/*  --------------------------------------------------------------------------
 *   Calculate objective function for SH problem
 * 
 *   
 *   D. Koehn
 *   Kiel, 22.12.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

double obj_sh(struct waveSH *waveSH, struct waveSH_PML *waveSH_PML, struct matSH *matSH, struct fwiSH *fwiSH, struct mpiPSV *mpiPSV, 
         struct seisSH *seisSH, struct seisSHfwi *seisSHfwi, struct acq *acq, float *hc, int nsrc, int nsrc_loc, int nsrc_glob, int ntr, 
         int ntr_glob, int ns, int itest, int iter, float **Ws, float **Wr, int hin, int *DTINV_help, float eps_scale, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern int RUN_MULTIPLE_SHOTS, TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, N_STREAMER, SEISMO, QUELLART, QUELLTYP, ORDER_SPIKE;
        extern int TIME_FILT, INV_STF, ORDER, L, MYID, LNORM, READREC, QUELLTYPB, LOG, NT;
        extern float FC_SPIKE_2,FC_SPIKE_1, FC, FC_START;

        /* local variables */
	double L2sum, L2_tmp;
        int ntr_loc, nt, ishot, nshots;
        FILE *FP;

	if ((MYID==0) && (LOG==1)) FP=stdout;

        /* initialization of L2 calculation */
	(*seisSHfwi).L2=0.0;
	(*seisSHfwi).energy=0.0;

	/* no differentiation of elastic and viscoelastic modelling because the viscoelastic parameters did not change during the forward modelling */
	matcopy_elastic_SH((*matSH).prho,(*matSH).pu);
	
	MPI_Barrier(MPI_COMM_WORLD);

	av_mu_SH((*matSH).pu, (*matSH).puip, (*matSH).pujp, (*matSH).prho);
	inv_rho_SH((*matSH).prho, (*matSH).prhoi);

	/* Preparing memory variables for update_s (viscoelastic) */
	if (L) prepare_update_s_visc_SH((*matSH).etajm, (*matSH).etaip, (*matSH).peta, (*matSH).fipjp, (*matSH).pujp, (*matSH).puip, 
				(*matSH).prho, (*matSH).ptaus, (*matSH).ptausipjp, (*matSH).f, (*matSH).g, (*matSH).bip, (*matSH).bjm, 
				(*matSH).cip, (*matSH).cjm, (*matSH).dip, (*matSH).d, (*matSH).e); 

		if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;

		for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){	

		if(MYID==0){
		   printf("\n=================================================================================================\n");
		   printf("\n *****  Starting simulation (test-forward model) no. %d for shot %d of %d (rel. step length %.8f) \n",itest,ishot,nshots,eps_scale);
		   printf("\n=================================================================================================\n\n");
		}
		  
		if(READREC==2){

		   if (SEISMO){
		      (*acq).recpos=receiver(FP, &ntr, ishot);
		      (*acq).recswitch = ivector(1,ntr);
		      (*acq).recpos_loc = splitrec((*acq).recpos,&ntr_loc, ntr, (*acq).recswitch);
		      ntr_glob=ntr;
		      ntr=ntr_loc;
		   }

	   	   /* Memory for seismic data */
	   	   alloc_seisSH(ntr,ns,seisSH);
	   
	   	   /* Memory for full data seismograms */
           	   alloc_seisSHfull(seisSH,ntr_glob);

	   	   /* Memory for FWI seismic data */ 
	   	   alloc_seisSHfwi(ntr,ntr_glob,ns,seisSHfwi);

		}

		for (nt=1;nt<=8;nt++) (*acq).srcpos1[nt][1]=(*acq).srcpos[nt][ishot]; 

		/* set QUELLTYP for each shot */
        	QUELLTYP = (*acq).srcpos[8][ishot];

			if (RUN_MULTIPLE_SHOTS){

				/* find this single source positions on subdomains */
				if (nsrc_loc>0) free_matrix((*acq).srcpos_loc,1,8,1,1);
				(*acq).srcpos_loc=splitsrc((*acq).srcpos1,&nsrc_loc, 1);
			}


			else{
				/* Distribute multiple source positions on subdomains */
				(*acq).srcpos_loc = splitsrc((*acq).srcpos,&nsrc_loc, nsrc);
			}

		MPI_Barrier(MPI_COMM_WORLD);

		/*==================================================================================
		           Starting simulation (forward model)
		==================================================================================*/
		
		/* calculate wavelet for each source point */
		(*acq).signals=NULL;
		(*acq).signals=wavelet((*acq).srcpos_loc,nsrc_loc,ishot);

		if (nsrc_loc){if(QUELLART==6){

		   /* time domain filtering of the source signal */
		   apply_tdfilt((*acq).signals,nsrc_loc,ns,ORDER_SPIKE,FC_SPIKE_2,FC_SPIKE_1);
		   
		   }
		}

		/* time domain filtering*/
		if ((TIME_FILT)&&(INV_STF==0)){
	
		   /* time domain filtering of the source signal */
		   apply_tdfilt((*acq).signals,nsrc_loc,ns,ORDER,FC,FC_START);

		}
						                                              
		/* solve forward problem */
		sh(waveSH,waveSH_PML,matSH,fwiSH,mpiPSV,seisSH,seisSHfwi,acq,hc,ishot,nshots,nsrc_loc,ns,ntr,Ws,Wr,hin,DTINV_help,2,req_send,req_rec);

		/* ===============================================
		   Calculate objective function and data residuals
		   =============================================== */
		if (ntr > 0){
		    calc_res_SH(seisSH,seisSHfwi,(*acq).recswitch,(*acq).recpos,(*acq).recpos_loc,ntr_glob,ntr,nsrc_glob,(*acq).srcpos,ishot,ns,iter,1);
		}				
		
	   if(READREC==2){

	     if (SEISMO) free_imatrix((*acq).recpos,1,3,1,ntr_glob);

	     if ((ntr>0) && (SEISMO)){

		   free_imatrix((*acq).recpos_loc,1,3,1,ntr);
		   (*acq).recpos_loc = NULL;
	 
		   switch (SEISMO){
		   case 1 : /* particle velocities only */
		           free_matrix((*seisSH).sectionvz,1,ntr,1,ns);
		           (*seisSH).sectionvz=NULL;
		           break;
		    }

	   }

	   free_matrix((*seisSHfwi).sectionread,1,ntr_glob,1,ns);
	   free_ivector((*acq).recswitch,1,ntr);	   
	   
	   if(QUELLTYPB){   
	      free_matrix((*seisSHfwi).sectionvzdata,1,ntr,1,ns);
	      free_matrix((*seisSHfwi).sectionvzdiff,1,ntr,1,ns);
	      free_matrix((*seisSHfwi).sectionvzdiffold,1,ntr,1,ns);
	   }

    	   if(SEISMO){
              free_matrix((*seisSH).fulldata_vz,1,ntr_glob,1,NT); 
           }
	   
	   ntr=0;
	   ntr_glob=0;
	 
	}

	nsrc_loc=0;

	} /* end of loop over shots */

	/* calculate L2 norm of all CPUs*/
	L2sum = 0.0;
        L2_tmp = (*seisSHfwi).L2;
	MPI_Allreduce(&L2_tmp,&L2sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   

        
        return L2sum;
	
}
