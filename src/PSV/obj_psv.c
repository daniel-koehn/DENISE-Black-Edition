/*  --------------------------------------------------------------------------
 *   Calculate objective function for PSV problem
 * 
 *   
 *   D. Koehn
 *   Kiel, 03.06.2016
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

float obj_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int nsrc, int nsrc_loc, int nsrc_glob, int ntr, 
         int ntr_glob, int ns, int itest, int iter, float **Ws, float **Wr, int hin, int *DTINV_help, float eps_scale, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern int RUN_MULTIPLE_SHOTS, TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, N_STREAMER, SEISMO, QUELLART, QUELLTYP, ORDER_SPIKE;
        extern int TIME_FILT, INV_STF, ORDER, L, MYID, LNORM;
        extern float FC_SPIKE_2,FC_SPIKE_1, FC, FC_START;

        /* local variables */
        float L2sum, L2_all_shots, energy_all_shots, energy_tmp, L2_tmp;
        int ntr_loc, nt, ishot, nshots;
        FILE *FP;

        /* initialization of L2 calculation */
	(*seisPSVfwi).L2=0.0;
	(*seisPSVfwi).energy=0.0;
	L2_all_shots=0.0;
	energy_all_shots=0.0;

	/* no differentiation of elastic and viscoelastic modelling because the viscoelastic parameters did not change during the forward modelling */
	matcopy_elastic((*matPSV).prho,(*matPSV).ppi,(*matPSV).pu);
	
	MPI_Barrier(MPI_COMM_WORLD);

	av_mue((*matPSV).pu,(*matPSV).puipjp,(*matPSV).prho);
	av_rho((*matPSV).prho,(*matPSV).prip,(*matPSV).prjp);

	/* Preparing memory variables for update_s (viscoelastic) */
	if (L) prepare_update_s((*matPSV).etajm,(*matPSV).etaip,(*matPSV).peta,(*matPSV).fipjp,(*matPSV).pu,(*matPSV).puipjp,(*matPSV).ppi,(*matPSV).prho,(*matPSV).ptaus,(*matPSV).ptaup,(*matPSV).ptausipjp,(*matPSV).f,(*matPSV).g,
				(*matPSV).bip,(*matPSV).bjm,(*matPSV).cip,(*matPSV).cjm,(*matPSV).dip,(*matPSV).d,(*matPSV).e);

		if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;

		for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){		

		if(MYID==0){
		   printf("\n=================================================================================================\n");
		   printf("\n *****  Starting simulation (test-forward model) no. %d for shot %d of %d (rel. step length %.8f) \n",itest,ishot,nshots,eps_scale);
		   printf("\n=================================================================================================\n\n");
		}
		  
		if(N_STREAMER>0){

		   if (SEISMO){
		      (*acq).recpos=receiver(FP, &ntr, ishot);
		      (*acq).recswitch = ivector(1,ntr);
		      (*acq).recpos_loc = splitrec((*acq).recpos,&ntr_loc, ntr, (*acq).recswitch);
		      ntr_glob=ntr;
		      ntr=ntr_loc;
		   }

		   /* Memory for seismic data */
		   alloc_seisPSV(ntr,ns,seisPSV);

		   /* Memory for FWI seismic data */ 
		   alloc_seisPSVfwi(ntr,ntr_glob,ns,seisPSVfwi);

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
		psv(wavePSV,wavePSV_PML,matPSV,fwiPSV,mpiPSV,seisPSV,seisPSVfwi,acq,hc,ishot,nshots,nsrc_loc,ns,ntr,Ws,Wr,hin,DTINV_help,2,req_send,req_rec);

		/* ===============================================
		   Calculate objective function and data residuals
		   =============================================== */
		if (ntr > 0){
		   calc_res_PSV(seisPSV,seisPSVfwi,(*acq).recswitch,(*acq).recpos,(*acq).recpos_loc,ntr_glob,ntr,nsrc_glob,(*acq).srcpos,ishot,ns,iter,1);
		}

	}

	/* calculate L2 norm of all CPUs*/
	L2sum = 0.0;
        L2_tmp = (*seisPSVfwi).L2;
	MPI_Allreduce(&L2_tmp,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	/* calculate L2 norm of all CPUs*/
	energy_all_shots = 0.0;
        energy_tmp = (*seisPSVfwi).energy;
	MPI_Allreduce(&energy_tmp,&energy_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	/* if(MYID==0){
		printf("L2sum: %e\n", L2sum);
		printf("energy_sum: %e\n\n", energy_all_shots);
	}*/

	if(LNORM==2){
	     L2sum = L2sum/energy_all_shots;
	}
	else{L2sum=L2sum;}    

        
        return L2sum;
	
}
