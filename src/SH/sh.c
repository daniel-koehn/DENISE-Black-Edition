/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D SH-forward problem by finite-differences 
 *   for a single shot 
 *
 *   mode = 0 - forward modelling only, STF estimation or FWI gradient calculation
 *   mode = 1 - adjoint modelling
 *   mode = 2 - evaluation objective function for step length estimation  
 * 
 *   
 *   D. Koehn
 *   Kiel, 13.08.2018
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void sh(struct waveSH *waveSH, struct waveSH_PML *waveSH_PML, struct matSH *matSH, struct fwiSH *fwiSH, struct mpiPSV *mpiPSV, 
         struct seisSH *seisSH, struct seisSHfwi *seisSHfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, int hin, int *DTINV_help, int mode, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, FDORDER, FW, L, GRAD_FORM, FC_SPIKE_1, FC_SPIKE_2, ORDER_SPIKE;
        extern int NX, NY, FREE_SURF, BOUNDARY, MODE, QUELLTYP, QUELLTYPB, QUELLART, FDORDER;
	extern int NPROCX, NPROCY, POS[3], NDT, SEISMO, IDXI, IDYI, GRAD_FORM, DTINV;
        extern int SNAP, INVMAT1, INV_STF, EPRECOND, NTDTINV, NXNYI, NT, ADJ_SIGN;
	extern FILE *FP;

        /* local variables */
	int i,j,l,nt,lsamp,lsnap,nsnap, nd, hin1, imat, imat1, imat2, infoout;
        float tmp, tmp1, muss, lamss, P3, P5;
        float hess_rho, hess_mu, hess_ts;
	float SUMr, SUMq;

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
		    
	/* initialize SH wavefields with zero */
	if (L){

                zero_denise_visc_SH(-nd+1,NY+nd,-nd+1,NX+nd, (*waveSH).pvz, (*waveSH).psxz, (*waveSH).psyz, (*waveSH).pvzm1, (*waveSH).pvzp1, 
				   (*waveSH_PML).psi_sxz_x, (*waveSH_PML).psi_syz_y, (*waveSH_PML).psi_vzx,  (*waveSH_PML).psi_vzy, 
                                   (*waveSH).pr, (*waveSH).pp, (*waveSH).pq, (*fwiSH).Rxz, (*fwiSH).Ryz);

	}else{	

                zero_denise_elast_SH(-nd+1,NY+nd,-nd+1,NX+nd, (*waveSH).pvz, (*waveSH).psxz, (*waveSH).psyz, (*waveSH).pvzm1, (*waveSH).pvzp1, 
				   (*waveSH_PML).psi_sxz_x, (*waveSH_PML).psi_syz_y, (*waveSH_PML).psi_vzx,  (*waveSH_PML).psi_vzy);
	
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
	   ADJ_SIGN=1;
        }

	if(mode==2){
	   ADJ_SIGN=1;
        }

        if(mode==1){
	   hin=1;
	   hin1=1;
	   ADJ_SIGN=-1;
        }

	for (nt=1;nt<=NT;nt++){     
		        
		/* Check if simulation is still stable */
		if (isnan((*waveSH).pvz[NY/2][NX/2])) {
		   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,(*waveSH).pvz[NY/2][NX/2]);
		   err(" Simulation is unstable !");}

	   infoout = !(nt%10000);

	   if (MYID==0){
	      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
	      /*time3=MPI_Wtime();*/
	   }

	      /* update of particle velocities */
              if(mode==0 || mode==2){
		update_v_PML_SH(1, NX, 1, NY, nt, (*waveSH).pvz, (*waveSH).pvzp1, (*waveSH).pvzm1, (*waveSH).uttz, (*waveSH).psxz, (*waveSH).psyz, (*matSH).prho, (*matSH).prhoi, (*acq).srcpos_loc, (*acq).signals, nsrc_loc, 
			       (*waveSH_PML).absorb_coeff,hc, infoout, 0, (*waveSH_PML).K_x, (*waveSH_PML).a_x, (*waveSH_PML).b_x, (*waveSH_PML).K_x_half, (*waveSH_PML).a_x_half, (*waveSH_PML).b_x_half, (*waveSH_PML).K_y, 
			       (*waveSH_PML).a_y, (*waveSH_PML).b_y, (*waveSH_PML).K_y_half, (*waveSH_PML).a_y_half, (*waveSH_PML).b_y_half, (*waveSH_PML).psi_sxz_x, (*waveSH_PML).psi_syz_y);
              }

              if(mode==1){
		update_v_PML_SH(1, NX, 1, NY, nt, (*waveSH).pvz, (*waveSH).pvzp1, (*waveSH).pvzm1, (*waveSH).uttz, (*waveSH).psxz, (*waveSH).psyz, (*matSH).prho, (*matSH).prhoi, (*acq).srcpos_loc_back, (*seisSHfwi).sectionvzdiff, ntr, 
			       (*waveSH_PML).absorb_coeff,hc, infoout, 1, (*waveSH_PML).K_x, (*waveSH_PML).a_x, (*waveSH_PML).b_x, (*waveSH_PML).K_x_half, (*waveSH_PML).a_x_half, (*waveSH_PML).b_x_half, (*waveSH_PML).K_y, 
			       (*waveSH_PML).a_y, (*waveSH_PML).b_y, (*waveSH_PML).K_y_half, (*waveSH_PML).a_y_half, (*waveSH_PML).b_y_half, (*waveSH_PML).psi_sxz_x, (*waveSH_PML).psi_syz_y);
              }
		                 
		/*if (MYID==0){
			time4=MPI_Wtime();
			time_av_v_update+=(time4-time3);
			if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
		}*/
		                                           
		/* exchange of particle velocities between PEs */
		exchange_v_SH((*waveSH).pvz, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);
		                                                       
		/*if (MYID==0){
		  time5=MPI_Wtime();
		  time_av_v_exchange+=(time5-time4);
		  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
		}*/                                                                                      	

	    if (L)    /* viscoelastic */

		update_s_visc_PML_SH(1, NX, 1, NY, (*waveSH).pvz, (*waveSH).uz, (*waveSH).uzx, (*waveSH).psyz, (*waveSH).psxz, (*matSH).pujp, (*matSH).puip, (*matSH).prho, hc, infoout,
				    (*waveSH).pr, (*waveSH).pp, (*waveSH).pq, (*matSH).fipjp, (*matSH).f, (*matSH).g, (*matSH).bip, (*matSH).bjm, (*matSH).cip, (*matSH).cjm, (*matSH).d, 
			            (*matSH).e, (*matSH).dip, (*waveSH_PML).K_x, (*waveSH_PML).a_x, (*waveSH_PML).b_x, (*waveSH_PML).K_x_half, (*waveSH_PML).a_x_half, (*waveSH_PML).b_x_half,
        			    (*waveSH_PML).K_y, (*waveSH_PML).a_y, (*waveSH_PML).b_y, (*waveSH_PML).K_y_half, (*waveSH_PML).a_y_half, (*waveSH_PML).b_y_half,
        			    (*waveSH_PML).psi_vzx, (*waveSH_PML).psi_vzy, fwiSH, mode);

	    else

		update_s_elastic_PML_SH(1, NX, 1, NY, (*waveSH).pvz, (*waveSH).uz, (*waveSH).uzx, (*waveSH).psyz, (*waveSH).psxz, (*matSH).pujp, (*matSH).puip, (*matSH).prho, hc,infoout,
				     (*waveSH_PML).K_x, (*waveSH_PML).a_x, (*waveSH_PML).b_x, (*waveSH_PML).K_x_half, (*waveSH_PML).a_x_half, (*waveSH_PML).b_x_half,
        			     (*waveSH_PML).K_y, (*waveSH_PML).a_y, (*waveSH_PML).b_y, (*waveSH_PML).K_y_half, (*waveSH_PML).a_y_half, (*waveSH_PML).b_y_half,
        			     (*waveSH_PML).psi_vzx, (*waveSH_PML).psi_vzy, mode);


	   /*if ((FREE_SURF) && (POS[2]==0)){
	   	if (L)   
			surface_visc_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).ppi, 
					    (*matPSV).pu, (*matPSV).prho, (*matPSV).ptaup, (*matPSV).ptaus, (*matPSV).etajm, (*matPSV).peta, hc, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
					    (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
		else  
	   		surface_elastic_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matPSV).ppi, (*matPSV).pu, (*matPSV).prho, hc, 
					       (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
	   }*/

	   /*if (MYID==0){
	      time6=MPI_Wtime();
		  time_av_s_update+=(time6-time5);
	      if (infoout)  fprintf(FP," stress exchange between PEs ...");
	      }*/


	   /* stress exchange between PEs */
	   exchange_s_SH((*waveSH).psxz,(*waveSH).psyz, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, 
	      (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);


	   /*if (MYID==0){
	      time7=MPI_Wtime();
	 	  time_av_s_exchange+=(time7-time6);
	     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
	      }  */

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO && (mode==0 || mode==2)){
			seismo_ssg(nt, ntr, (*acq).recpos_loc, (*seisSH).sectionvz, (*seisSH).sectionvz, (*seisSH).sectionvz, (*seisSH).sectionvz, (*seisSH).sectionvz, 
				(*waveSH).pvz, (*waveSH).pvz, (*waveSH).pvz, (*waveSH).pvz, (*matSH).pu, (*matSH).pu, (*matSH).prho, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (nt==lsnap) && (nt<=iround(TSNAP2/DT))){

	      snap(FP,nt,++nsnap,(*waveSH).pvz,(*waveSH).pvz,(*waveSH).pvz,(*waveSH).pvz,(*matSH).pu,(*matSH).pu,hc);

	      lsnap=lsnap+iround(TSNAPINC/DT);
	   }

	      
	   /*if (MYID==0){
	      time8=MPI_Wtime();
		  time_av_timestep+=(time8-time3);
	      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
	      } */ 

	if((nt==hin1)&&(mode==0)&&(MODE==1||MODE==2)){

	  /* store forward wavefields for time-domain inversion and RTM */
          /* ---------------------------------------------------------- */	    
	    
		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){

			/* store forward wavefield for density gradient */
		        (*fwiSH).forward_prop_rho_z[imat2] = (*waveSH).pvzp1[j][i];

			/* store forward wavefield for mu and Taus gradient */
	 	        (*fwiSH).forward_prop_sxz[imat2] = (*waveSH).uz[j][i];
	 	        (*fwiSH).forward_prop_syz[imat2] = (*waveSH).uzx[j][i];

			for (l=1;l<=L;l++){
			    (*fwiSH).forward_prop_rxz[imat2][l] = (*waveSH).pr[j][i][l] + (*fwiSH).tausl[l] * (*fwiSH).rxzt[j][i][l];
	 	            (*fwiSH).forward_prop_ryz[imat2][l] = (*waveSH).pq[j][i][l] + (*fwiSH).tausl[l] * (*fwiSH).ryzt[j][i][l];
			}			

			/* Compute Pseudo-Hessian main diagonal approximation */
			if(EPRECOND==4){

			   /* Pseudo-Hessian density */
			   hess_rho = (*waveSH).pvz[j][i] * (*waveSH).pvzp1[j][i];			  

			   /* Pseudo-Hessian Vs and Taus*/
			   muss = (*matSH).prho[j][i] * (*matSH).pu[j][i] * (*matSH).pu[j][i];
			   
			   SUMr=SUMq=0.0;
                           for (l=1;l<=L;l++){
			       SUMr += (*fwiSH).Rxz[j][i][l];
                               SUMq += (*fwiSH).Ryz[j][i][l];
			   }
			   
			   P3 = (((*waveSH).psxz[j][i]-SUMr) * (*waveSH).uz[j][i]) + (((*waveSH).psyz[j][i]-SUMq) * (*waveSH).uzx[j][i]);
			   P5 = 0.0;  			            
		           for (l=1;l<=L;l++){
			       P5 += ((*fwiSH).forward_prop_rxz[imat2][l] * (*fwiSH).Rxz[j][i][l]) + ((*fwiSH).forward_prop_ryz[imat2][l] * (*fwiSH).Ryz[j][i][l]); 
			   }

			   /* Pseudo-Hessian main-diagonal and non-diagonal elements */
			   /* ------------------------------------------------------ */
			   (*fwiSH).hess_rho2[j][i] += hess_rho * hess_rho;			   

		           if(muss>0.0){
			   
			      hess_mu = - (*fwiSH).c1mu[j][i] * P3 + (*fwiSH).c4mu[j][i] * P5;
			      hess_ts = - (*fwiSH).c1ts[j][i] * P3 + (*fwiSH).c4ts[j][i] * P5;
			   			
		                (*fwiSH).hess_mu2[j][i] += hess_mu * hess_mu;
			        (*fwiSH).hess_ts2[j][i] += hess_ts * hess_ts;
			       (*fwiSH).hess_muts[j][i] += hess_mu * hess_ts;
			      (*fwiSH).hess_murho[j][i] += hess_mu * hess_rho;
			      (*fwiSH).hess_tsrho[j][i] += hess_rho * hess_ts; 
			       
		           }			   			   			    		 

			 }
			 
			 imat2++;			 			 
		    
		    }
		}
	   
	    if((EPRECOND==1)||(EPRECOND==3)){
	      eprecond_SH(Ws,(*waveSH).pvz);
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
		                                   
			   	(*fwiSH).waveconv_rho_shot[j][i] += (*waveSH).pvz[j][i] * (*fwiSH).forward_prop_rho_z[imat];
			
				/* Vs/Taus-gradient (symmetrized stress-velocity formulation) according to Fabien-Ouellet et al. (2017) */
		                if(GRAD_FORM==2){			

		           	    if(INVMAT1==1){
		               	        muss = (*matSH).prho[j][i] * (*matSH).pu[j][i] * (*matSH).pu[j][i];
	                   	    }
	           
		           	    if(INVMAT1==3){
		               	        muss = (*matSH).pu[j][i];
				    }

				    /* correlate forward and adjoint wavefields */
				    P3 = ((*fwiSH).forward_prop_sxz[imat] * (*waveSH).uz[j][i]) + ((*fwiSH).forward_prop_syz[imat] * (*waveSH).uzx[j][i]);
				    
				    P5 = 0.0;  			            
				    for (l=1;l<=L;l++){
					P5 += ((*fwiSH).forward_prop_rxz[imat][l] * (*fwiSH).Rxz[j][i][l]) + ((*fwiSH).forward_prop_ryz[imat][l] * (*fwiSH).Ryz[j][i][l]); 
				    }    

				    if(muss>0.0){			
				         (*fwiSH).waveconv_u_shot[j][i] += - (*fwiSH).c1mu[j][i] * P3 + (*fwiSH).c4mu[j][i] * P5;
			  		//(*fwiSH).waveconv_ts_shot[j][i] += - (*fwiSH).c1ts[j][i] * P3 + (*fwiSH).c4ts[j][i] * P5;
					(*fwiSH).waveconv_ts_shot[j][i] = 0.0;
				    } 		                  
		                }			
						                                                                                                     
			imat++;
			}
		    }  
		  
		  if(EPRECOND==1){
		     eprecond_SH(Wr,(*waveSH).pvz);
		  }
		                                                                                                                       
	    hin++;
	    }

	   }/*--------------------  End  of loop over timesteps ----------*/		

}
