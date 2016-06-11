/*------------------------------------------------------------------------
 * Calculate gradient and objective function for PSV problem
 *
 *
 * D. Koehn
 * Kiel, 05.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"


float grad_obj_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int iter, int nsrc, int ns, int ntr, int ntr_glob, int nsrc_glob, 
         int nsrc_loc, int ntr_loc, int nstage, float **We, float **Ws, float **Wr, float ** taper_coeff, int hin, int *DTINV_help, 
         MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern int MYID, TIME_FILT, IDX, IDY, NX, NY, NT, RUN_MULTIPLE_SHOTS, INV_STF, QUELLART;
        extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, SEISMO, EPRECOND, LNORM;
        extern int N_STREAMER, SWS_TAPER_CIRCULAR_PER_SHOT, QUELLTYPB, QUELLTYP, LOG;
        extern int ORDER_SPIKE, ORDER, SHOTINC;
        extern float EPSILON, FC, FC_START, FC_SPIKE_1, FC_SPIKE_2;
        extern float C_vp, C_vs, C_rho;
        extern char MFILE[STRING_SIZE];

        /* local variables */
	int i, j, nshots, ishot, nt, lsnap, itestshot, swstestshot;
        float L2sum, L2_all_shots, energy_all_shots, energy_tmp, L2_tmp;
        char source_signal_file[STRING_SIZE];

	FILE *FP;

	if ((MYID==0) && (LOG==1)) FP=stdout;

	/* initialization of L2 calculation */
	(*seisPSVfwi).L2=0.0;
	(*seisPSVfwi).energy=0.0;
	L2_all_shots=0.0;
	energy_all_shots=0.0;

	EPSILON=0.0;  /* test step length */

	/* set gradient and preconditioning matrices 0 before next iteration*/
	init_grad((*fwiPSV).waveconv);
	init_grad((*fwiPSV).waveconv_rho);
	init_grad((*fwiPSV).waveconv_u);

	itestshot=TESTSHOT_START;
	swstestshot=0;
        SHOTINC=1;

	if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;

	for (ishot=1;ishot<=nshots;ishot+=SHOTINC){
	/*for (ishot=1;ishot<=1;ishot+=1){*/

	/*initialize gradient matrices for each shot with zeros*/
	init_grad((*fwiPSV).waveconv_shot);
	init_grad((*fwiPSV).waveconv_u_shot);
	init_grad((*fwiPSV).waveconv_rho_shot);

	if((EPRECOND==1)||(EPRECOND==3)){
	   init_grad(Ws);
	   init_grad(Wr);
	   init_grad(We);
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
		        Estimate source time function by Wiener deconvolution
	==================================================================================*/

	if((INV_STF)&&(iter==1)){
	  stf_psv(wavePSV,wavePSV_PML,matPSV,fwiPSV,mpiPSV,seisPSV,seisPSVfwi,acq,hc,ishot,nshots,nsrc_loc,nsrc,ns,ntr,ntr_glob,iter,Ws,Wr,hin,DTINV_help,req_send,req_rec);
	}
	 
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

	/*printf("MYID=%d, nsrc_loc = %d \n",MYID,nsrc_loc);*/

	/*char  source_signal_file[STRING_SIZE];
	sprintf(source_signal_file,"source_signal.%d.su.shot%d.it%d",MYID,ishot,iter);
	fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
	output_source_signal(fopen(source_signal_file,"w"),signals,NT,3);*/

	/* output source signal e.g. for cross-correlation of comparison with analytical solutions */
	if(RUN_MULTIPLE_SHOTS){

		if(nsrc_loc>0){
			   sprintf(source_signal_file,"%s_source_signal.%d.su.shot%d", MFILE, MYID,ishot);
			   fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
			   output_source_signal(fopen(source_signal_file,"w"),(*acq).signals,NT,1);
		}                                
		                        
		MPI_Barrier(MPI_COMM_WORLD);
	}
				                                                      
	/* solve forward problem */
	psv(wavePSV,wavePSV_PML,matPSV,fwiPSV,mpiPSV,seisPSV,seisPSVfwi,acq,hc,ishot,nshots,nsrc_loc,ns,ntr,Ws,Wr,hin,DTINV_help,0,req_send,req_rec);

	/* ===============================================
	   Calculate objective function and data residuals
	   =============================================== */

	if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){swstestshot=1;}

	if (ntr > 0){
	   calc_res_PSV(seisPSV,seisPSVfwi,(*acq).recswitch,(*acq).recpos,(*acq).recpos_loc,ntr_glob,ntr,nsrc_glob,(*acq).srcpos,ishot,ns,iter,swstestshot);
	}

	if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){
	       swstestshot=0;
	       itestshot+=TESTSHOT_INCR;
	}

	/* output of time reversed residual seismograms */
	if ((SEISMO)&&(iter==1)&&(ishot==1)){
	   outseis_PSVres(seisPSV,seisPSVfwi,(*acq).recswitch,(*acq).recpos,(*acq).recpos_loc,ntr_glob,(*acq).srcpos,ishot,ns,nstage,FP);       
	}
		          	    		    
	/*================================================================================
		        Starting simulation (backward model)
	==================================================================================*/
	    
	    /* Distribute multiple source positions on subdomains */
	    /* define source positions at the receivers */
	    (*acq).srcpos_loc_back = matrix(1,6,1,ntr);
	    for (i=1;i<=ntr;i++){
		(*acq).srcpos_loc_back[1][i] = ((*acq).recpos_loc[1][i]);
		(*acq).srcpos_loc_back[2][i] = ((*acq).recpos_loc[2][i]);
	    }
		                            
	   /* solve adjoint problem */
	   psv(wavePSV,wavePSV_PML,matPSV,fwiPSV,mpiPSV,seisPSV,seisPSVfwi,acq,hc,ishot,nshots,ntr,ns,ntr,Ws,Wr,hin,DTINV_help,1,req_send,req_rec);               

	   /* assemble PSV gradients for each shot */
	   ass_gradPSV(fwiPSV,matPSV,iter);

	if((EPRECOND==1)||(EPRECOND==3)){
	  /* calculate energy weights */
	  eprecond1(We,Ws,Wr);
	      
	  /* scale gradient with energy weights*/
	  for(i=1;i<=NX;i=i+IDX){
	      for(j=1;j<=NY;j=j+IDY){

		     (*fwiPSV).waveconv_shot[j][i] = (*fwiPSV).waveconv_shot[j][i]/(We[j][i]*C_vp*C_vp);
		     (*fwiPSV).waveconv_u_shot[j][i] = (*fwiPSV).waveconv_u_shot[j][i]/(We[j][i]*C_vs*C_vs);
		     if(C_vs==0.0){(*fwiPSV).waveconv_u_shot[j][i] = 0.0;}
		     (*fwiPSV).waveconv_rho_shot[j][i] = (*fwiPSV).waveconv_rho_shot[j][i]/(We[j][i]*C_rho*C_rho);

	      }
	  }
	}

	if (SWS_TAPER_CIRCULAR_PER_SHOT){    /* applying a circular taper at the source position to the gradient of each shot */
	
		/* applying the preconditioning */
		taper_grad_shot((*fwiPSV).waveconv_shot,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,ishot);
		taper_grad_shot((*fwiPSV).waveconv_rho_shot,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,ishot);
		taper_grad_shot((*fwiPSV).waveconv_u_shot,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,ishot);
	
	} /* end of SWS_TAPER_CIRCULAR_PER_SHOT == 1 */

	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			(*fwiPSV).waveconv[j][i] += (*fwiPSV).waveconv_shot[j][i];
			(*fwiPSV).waveconv_rho[j][i] += (*fwiPSV).waveconv_rho_shot[j][i];
			(*fwiPSV).waveconv_u[j][i] += (*fwiPSV).waveconv_u_shot[j][i];
		}
	}

	if(N_STREAMER>0){

	   if (SEISMO) free_imatrix((*acq).recpos,1,3,1,ntr_glob);

	   if ((ntr>0) && (SEISMO)){

		   free_imatrix((*acq).recpos_loc,1,3,1,ntr);
		   (*acq).recpos_loc = NULL;
	 
		   switch (SEISMO){
		   case 1 : /* particle velocities only */
		           free_matrix((*seisPSV).sectionvx,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectionvy,1,ntr,1,ns);
		           (*seisPSV).sectionvx=NULL;
		           (*seisPSV).sectionvy=NULL;
		           break;
		    case 2 : /* pressure only */
		           free_matrix((*seisPSV).sectionp,1,ntr,1,ns);
		           break;
		    case 3 : /* curl and div only */
		           free_matrix((*seisPSV).sectioncurl,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectiondiv,1,ntr,1,ns);
		           break;
		    case 4 : /* everything */
		           free_matrix((*seisPSV).sectionvx,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectionvy,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectionp,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectioncurl,1,ntr,1,ns);
		           free_matrix((*seisPSV).sectiondiv,1,ntr,1,ns);
		           break;

		    }

	   }

	   free_matrix((*seisPSVfwi).sectionread,1,ntr_glob,1,ns);
	   free_ivector((*acq).recswitch,1,ntr);
	   
	   if((QUELLTYPB==1)||(QUELLTYPB==3)){
	      free_matrix((*seisPSVfwi).sectionvxdata,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionvxdiff,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionvxdiffold,1,ntr,1,ns);
	   }
	   
	   if((QUELLTYPB==1)||(QUELLTYPB==2)){   
	      free_matrix((*seisPSVfwi).sectionvydata,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionvydiff,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionvydiffold,1,ntr,1,ns);
	   }
	   
	   if(QUELLTYPB==4){   
	      free_matrix((*seisPSVfwi).sectionpdata,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionpdiff,1,ntr,1,ns);
	      free_matrix((*seisPSVfwi).sectionpdiffold,1,ntr,1,ns);
	   }
	   
	 
	}

	nsrc_loc=0;

	} /* end of loop over shots (forward and backpropagation) */   

	/* calculate L2 norm of all CPUs*/
	L2sum = 0.0;
        L2_tmp = (*seisPSVfwi).L2;
	MPI_Allreduce(&L2_tmp,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	/* calculate L2 norm of all CPUs*/
	energy_all_shots = 0.0;
        energy_tmp = (*seisPSVfwi).energy;
	MPI_Allreduce(&energy_tmp,&energy_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	/*if(MYID==0){
		printf("L2sum: %e\n", L2sum);
		printf("energy_sum: %e\n\n", energy_all_shots);

	}*/

	if(LNORM==2){
	     L2sum = L2sum/energy_all_shots;
	}
	else{L2sum=L2sum;}

return L2sum;

}

