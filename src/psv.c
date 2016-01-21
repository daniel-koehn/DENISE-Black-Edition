/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D PSV-forward problem by finite-differences 
 *   for a single shot 
 * 
 *   
 *   D. Koehn
 *   Kiel, 28.12.2015
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, float ** ppi, float ** pu, 
        float ** puipjp, float **prho, float  **prip, float **prjp, 
        float *hc, int infoout, float **fipjp, float **f, float **g, float *bip, float *bjm, 
        float *cip, float *cjm, float ***d, float ***e,  float ***dip, float **ptaup, float **ptaus, 
        float *etajm, float *peta, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
        float ** buffertop_to_bot, float ** bufferbot_to_top, int ishot, int nshots, int nsrc_loc, 
        float ** srcpos_loc, int ** recpos_loc, float ** signals, int ns, int ntr, float **sectionp, 
        float **sectionvx, float **sectionvy, float **sectiondiv, float **sectioncurl, 
        float *forward_prop_rho_x, float *forward_prop_rho_y, float *forward_prop_x, float *forward_prop_y, 
        float *forward_prop_u, float **waveconv_shot, float **waveconv_u_shot, float **waveconv_rho_shot, 
        float **Ws, float **Wr, float **sectionvxdiff, float **sectionvydiff, int hin, int *DTINV_help, 
        int mode, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, FDORDER, FW, L, GRAD_FORM, FC_SPIKE_1, FC_SPIKE_2, ORDER_SPIKE;
        extern int NX, NY, FREE_SURF, BOUNDARY, INVMAT, QUELLTYP, QUELLART, FDORDER;
	extern int NPROCX, NPROCY, POS[3], NDT, SEISMO, IDXI, IDYI, GRAD_FORM, DTINV;
        extern int SNAP, INVMAT1, INV_STF, EPRECOND, NTDTINV, NXNYI, NT;
	extern FILE *FP;

        /* local variables */
	int i,j,nt,lsamp,lsnap,nsnap, nd, hin1, imat, imat1, imat2;
        float tmp, tmp1, muss, lamss;

        nd = FDORDER/2 + 1;

	/*MPI_Barrier(MPI_COMM_WORLD);*/
        
        if(MYID==0){

		if((INV_STF==0)&&(mode==0)){

		   fprintf(FP,"\n==================================================================================\n");
		   fprintf(FP,"\n *****  Starting simulation (forward model) for shot %d of %d  ********** \n",ishot,nshots);
		   fprintf(FP,"\n==================================================================================\n\n");

		}


		if((INV_STF==1)&&(mode==0)){

		   fprintf(FP,"\n==================================================================================\n");
		   fprintf(FP,"\n *****  Starting simulation (STF) for shot %d of %d  ********** \n",ishot,nshots);
		   fprintf(FP,"\n==================================================================================\n\n");

		}

		if(mode==1){

		   fprintf(FP,"\n==================================================================================\n");
		   fprintf(FP,"\n *****  Starting simulation (adjoint wavefield)  ********** \n");
		   fprintf(FP,"\n==================================================================================\n\n");

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
			    
	/* initialize PSV wavefields with zero */
	if (L){
		zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd, (*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                                 (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,(*wavePSV_PML).psi_sxy_x,
                                 (*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,(*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,
                                 (*wavePSV_PML).psi_vxxs,(*wavePSV).pr,(*wavePSV).pp,(*wavePSV).pq);
	}else{	
		zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                            (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,
                            (*wavePSV_PML).psi_sxy_x,(*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,
                            (*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,(*wavePSV_PML).psi_vxxs);	
	}                                                         
	     
	/*----------------------  loop over timesteps (forward model) ------------------*/

	lsnap=iround(TSNAP1/DT);  
	lsamp=NDT;
	nsnap=0;

        if(mode==0){
	   hin=1;
	   hin1=1;
	   imat=1;
	   imat1=1;
	   imat2=1;
        }

        if(mode==1){
	   hin=1;
	   hin1=1;
        }

	for (nt=1;nt<=NT;nt++){     
		        
		/* Check if simulation is still stable */
		/*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
		if (isnan((*wavePSV).pvy[NY/2][NX/2])) {
		   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,(*wavePSV).pvy[NY/2][NX/2]);
		   err(" Simulation is unstable !");}

		
	   infoout = !(nt%10000);

	   if (MYID==0){
	      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
	      /*time3=MPI_Wtime();*/
	   }

	      /* update of particle velocities */
              if(mode==0){
	         update_v_PML(1, NX, 1, NY, nt, (*wavePSV).pvx, (*wavePSV).pvxp1, (*wavePSV).pvxm1, (*wavePSV).pvy, (*wavePSV).pvyp1, (*wavePSV).pvym1, (*wavePSV).uttx, (*wavePSV).utty, (*wavePSV).psxx, (*wavePSV).psyy,       
                              (*wavePSV).psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,(*wavePSV_PML).absorb_coeff,hc,infoout, mode, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, 
                              (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_sxx_x, (*wavePSV_PML).psi_syy_y, 
                              (*wavePSV_PML).psi_sxy_y, (*wavePSV_PML).psi_sxy_x);
              }

              if(mode==1){
	         update_v_PML(1, NX, 1, NY, nt, (*wavePSV).pvx, (*wavePSV).pvxp1, (*wavePSV).pvxm1, (*wavePSV).pvy, (*wavePSV).pvyp1, (*wavePSV).pvym1, (*wavePSV).uttx, (*wavePSV).utty, (*wavePSV).psxx, (*wavePSV).psyy, 
                              (*wavePSV).psxy, prip, prjp, srcpos_loc, sectionvxdiff, sectionvydiff,nsrc_loc,(*wavePSV_PML).absorb_coeff,hc,infoout, mode, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half,  
                              (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_sxx_x, 
                              (*wavePSV_PML).psi_syy_y, (*wavePSV_PML).psi_sxy_y, (*wavePSV_PML).psi_sxy_x);
              }
		                 
		/*if (MYID==0){
			time4=MPI_Wtime();
			time_av_v_update+=(time4-time3);
			if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
		}*/
		                                           
		/* exchange of particle velocities between PEs */
		exchange_v((*wavePSV).pvx,(*wavePSV).pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);
		                                                       
		/*if (MYID==0){
		  time5=MPI_Wtime();
		  time_av_v_exchange+=(time5-time4);
		  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
		}*/                                                                                      	

	    if (L)    /* viscoelastic */
	    	update_s_visc_PML(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, ppi, pu, puipjp, prho, hc,
                                  infoout, (*wavePSV).pr, (*wavePSV).pp, (*wavePSV).pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, 
                                  (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx, (*wavePSV_PML).psi_vyy, 
                                  (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx);
	    else
	   	update_s_elastic_PML(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, ppi, pu, puipjp, 
                                     (*wavePSV_PML).absorb_coeff, prho, hc, infoout, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, 
                                     (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx, (*wavePSV_PML).psi_vyy, (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx, mode);  


	    /* explosive source */
	   if (QUELLTYP==1) 	
	   psource(nt,(*wavePSV).psxx,(*wavePSV).psyy,srcpos_loc,signals,nsrc_loc,0);
	   
	   /* moment tensor source */
	   if (QUELLTYP==5) 	
	   msource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,srcpos_loc,signals,nsrc_loc,0);

	   if ((FREE_SURF) && (POS[2]==0)){
	   	if (L)    /* viscoelastic */
			surface_PML(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*wavePSV).pp, (*wavePSV).pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, (*wavePSV_PML).K_x, 
                                    (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
		else      /* elastic */
	   		surface_elastic_PML(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, ppi, pu, prho, hc, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
	   }


	   /*if (MYID==0){
	      time6=MPI_Wtime();
		  time_av_s_update+=(time6-time5);
	      if (infoout)  fprintf(FP," stress exchange between PEs ...");
	      }*/


	   /* stress exchange between PEs */
	    exchange_s((*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy, 
	      bufferlef_to_rig, bufferrig_to_lef, 
	      buffertop_to_bot, bufferbot_to_top,
	      req_send, req_rec);

	   /*if (MYID==0){
	      time7=MPI_Wtime();
	 	  time_av_s_exchange+=(time7-time6);
	     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
	      }  */

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO){
			seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy, 
				sectionp, sectioncurl, sectiondiv, 
				(*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, ppi, pu, prho, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

	      snap(FP,nt,++nsnap,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,pu,ppi,hc);

	      lsnap=lsnap+iround(TSNAPINC/DT);
	   }

	      
	   /*if (MYID==0){
	      time8=MPI_Wtime();
		  time_av_timestep+=(time8-time3);
	      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
	      } */  		

	if((nt==hin1)&&(mode==0)){

	    /* save forward wavefields for time-domain inversion */
            /* ------------------------------------------------- */
	    if(INVMAT==0){
	    
		for (i=1;i<=NX;i=i+IDXI){
		    for (j=1;j<=NY;j=j+IDYI){
			 forward_prop_rho_x[imat1]=(*wavePSV).pvxp1[j][i];
			 forward_prop_rho_y[imat1]=(*wavePSV).pvyp1[j][i];
		         imat1++;                                   
		    }
		}   

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){
		    
			/* gradients with data integration */
		        if(GRAD_FORM==1){
			   forward_prop_x[imat]=(*wavePSV).psxx[j][i];
			   forward_prop_y[imat]=(*wavePSV).psyy[j][i];
		        }
		    
			/* gradients without data integration */
			if(GRAD_FORM==2){
		           forward_prop_x[imat]=(*wavePSV).ux[j][i];
			   forward_prop_y[imat]=(*wavePSV).uy[j][i];
			}
		    
			imat++;
		    }
		 }
	    

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){
		    
			/* gradients with data integration */
		        if(GRAD_FORM==1){
	 	          forward_prop_u[imat2]=(*wavePSV).psxy[j][i];
			}

			/* gradients without data integration */
		        if(GRAD_FORM==2){
	 	          forward_prop_u[imat2]=(*wavePSV).uxy[j][i];
		        }

			imat2++;
		    
		    }
		}
	
	    }
	   
	    if((EPRECOND==1)||(EPRECOND==3)){
	      eprecond(Ws,(*wavePSV).pvx,(*wavePSV).pvy);
	    }
	 
	    hin++;
	    hin1=hin1+DTINV;

            DTINV_help[nt]=1;
		                                                     
	}

	/* save backpropagated wavefields for time-domain inversion and partially assemble gradients */ 
        /* ----------------------------------------------------------------------------------------- */
	   if((mode==1)&&(DTINV_help[NT-nt+1]==1)){
	    
		imat=((NXNYI*(NTDTINV)) - hin*NXNYI)+1;
	                
		    for (i=1;i<=NX;i=i+IDXI){   
			for (j=1;j<=NY;j=j+IDYI){ 
		                                   
			   	waveconv_rho_shot[j][i]+=((*wavePSV).pvxp1[j][i]*forward_prop_rho_x[imat])+((*wavePSV).pvyp1[j][i]*forward_prop_rho_y[imat]);
			
			   	/* mu-gradient with data integration */
			   	if(GRAD_FORM==1){
			
				   waveconv_shot[j][i]+= (forward_prop_x[imat]+forward_prop_y[imat])*((*wavePSV).psxx[j][i]+(*wavePSV).psyy[j][i]);			  

		                   if(INVMAT1==1){
				       muss = prho[j][i] * pu[j][i] * pu[j][i];
				       lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
				   }
			   
				   if(INVMAT1==3){
				       muss = pu[j][i];
				       lamss = ppi[j][i]; 
				   } 
			                
				   if(muss>0.0){
				      waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[imat] * (*wavePSV).psxy[j][i])) 
		                          + ((1.0/4.0) * ((forward_prop_x[imat] + forward_prop_y[imat]) * ((*wavePSV).psxx[j][i] + (*wavePSV).psyy[j][i])) / ((lamss+muss)*(lamss+muss)))  
		                          + ((1.0/4.0) * ((forward_prop_x[imat] - forward_prop_y[imat]) * ((*wavePSV).psxx[j][i] - (*wavePSV).psyy[j][i])) / (muss*muss));
				   }
				   
		                }

				/* Vs-gradient without data integration (stress-velocity in non-conservative form) */
		                if(GRAD_FORM==2){
			
				   waveconv_shot[j][i]+= (forward_prop_x[imat]+forward_prop_y[imat])*((*wavePSV).psxx[j][i]+(*wavePSV).psyy[j][i]);

		                   if(INVMAT1==1){
				       muss = prho[j][i] * pu[j][i] * pu[j][i];
				       lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
				   }
			   
				   if(INVMAT1==3){
				       muss = pu[j][i];
				       lamss = ppi[j][i]; 
				   } 
			                
				   if(muss>0.0){
								
				      tmp = (1.0/(4.0*(lamss+muss)*(lamss+muss))) - (1.0/(4.0*muss*muss));
				      tmp1 = (1.0/(4.0*(lamss+muss)*(lamss+muss))) + (1.0/(4.0*muss*muss));
			
				      waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[imat] * (*wavePSV).psxy[j][i])) 
		                                           + ( tmp1 * (forward_prop_x[imat] * (*wavePSV).psxx[j][i] + forward_prop_y[imat] * (*wavePSV).psyy[j][i]))  
		                                           + ( tmp  * (forward_prop_x[imat] * (*wavePSV).psyy[j][i] + forward_prop_y[imat] * (*wavePSV).psxx[j][i]));	  
					  
				   } 
		                  
		                }			
						                                                                                                     
			   imat++;
			   }
		    }  
		  
		  if(EPRECOND==1){
		     eprecond(Wr,(*wavePSV).pvx,(*wavePSV).pvy);
		  }
		                                                                                                                       
	    hin++;
	    }

	   /* save forward wavefields for time-Laplace domain inversion */   
	   /*if(INVMAT==1){
	     
	     time = (float)(nt*DT);
	     tmp = exp(-GAMMA*(time+DT)*(time+DT));

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){

			forward_propl_rho_x[j][i]+=pvxp1[j][i]*tmp;
			forward_propl_rho_y[j][i]+=pvyp1[j][i]*tmp;

		        if(GRAD_FORM==1){
	 	          forward_propl_u[j][i]+=psxy[j][i]*tmp;
			  forward_propl_x[j][i]+=psxx[j][i]*tmp;
			  forward_propl_y[j][i]+=psyy[j][i]*tmp;
			}

		        if(GRAD_FORM==2){
	 	          forward_propl_u[j][i]+=uxy[j][i]*tmp;
			  forward_propl_x[j][i]+=ux[j][i]*tmp;
			  forward_propl_y[j][i]+=uy[j][i]*tmp;
		        }
		
		    }
		}
	
	    }*/

    	  /* save backpropagated wavefields for time-Laplace domain inversion */   
          /*if(INVMAT==1){

             for (i=1;i<=NX;i=i+IDXI){ 
	        for (j=1;j<=NY;j=j+IDYI){
	    
	            back_prop_rho_x[j][i]+=pvxp1[j][i];
		    back_prop_rho_y[j][i]+=pvyp1[j][i];

                    if(GRAD_FORM==1){
 	              back_prop_u[j][i]+=psxy[j][i];
		      back_prop_x[j][i]+=psxx[j][i];
	              back_prop_y[j][i]+=psyy[j][i];
	            }

                    if(GRAD_FORM==2){
 	              back_prop_u[j][i]+=uxy[j][i];
		      back_prop_x[j][i]+=ux[j][i];
	              back_prop_y[j][i]+=uy[j][i];
                    }
		
                 }
             }
	
          }*/

	   }/*--------------------  End  of loop over timesteps ----------*/		

}
