/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D PSV TTI-forward problem by finite-differences 
 *   for a single shot 
 *
 *   mode = 0 - forward modelling only, STF estimation or FWI gradient calculation
 *   mode = 1 - backpropagation of data residuals
 *   mode = 2 - evaluation objective function for step length estimation  
 * 
 *   
 *   D. Koehn
 *   Kiel, 05.02.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void TTI(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matTTI *matTTI, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, int hin, int *DTINV_help, int mode, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, FDORDER, FW, L, GRAD_FORM, FC_SPIKE_1, FC_SPIKE_2, ORDER_SPIKE;
        extern int NX, NY, FREE_SURF, BOUNDARY, MODE, QUELLTYP, QUELLTYPB, QUELLART, FDORDER;
	extern int NPROCX, NPROCY, POS[3], NDT, SEISMO, IDXI, IDYI, GRAD_FORM, DTINV;
        extern int SNAP, INVMAT1, INV_STF, EPRECOND, NTDTINV, NXNYI, NT;
	extern FILE *FP;

        /* local variables */
	int i,j,nt,lsamp,lsnap,nsnap, nd, hin1, imat, imat1, imat2, infoout;
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
			    
	/* initialize PSV wavefields with zero */
	if (L){
		zero_denise_visc_PSV(-nd+1,NY+nd,-nd+1,NX+nd, (*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                                 (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,(*wavePSV_PML).psi_sxy_x,
                                 (*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,(*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,
                                 (*wavePSV_PML).psi_vxxs,(*wavePSV).pr,(*wavePSV).pp,(*wavePSV).pq);
	}else{	
		zero_denise_elast_PSV(-nd+1,NY+nd,-nd+1,NX+nd,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
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
              if(mode==0 || mode==2){
	         update_v_PML_PSV(1, NX, 1, NY, nt, (*wavePSV).pvx, (*wavePSV).pvxp1, (*wavePSV).pvxm1, (*wavePSV).pvy, (*wavePSV).pvyp1, (*wavePSV).pvym1, (*wavePSV).uttx, (*wavePSV).utty, (*wavePSV).psxx, (*wavePSV).psyy,       
                              (*wavePSV).psxy, (*matTTI).prip, (*matTTI).prjp, (*acq).srcpos_loc,(*acq).signals,(*acq).signals,nsrc_loc,(*wavePSV_PML).absorb_coeff,hc,infoout, 0, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
                              (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, 
                              (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_sxx_x, (*wavePSV_PML).psi_syy_y, (*wavePSV_PML).psi_sxy_y, (*wavePSV_PML).psi_sxy_x);
              }

              if(mode==1){
	         update_v_PML_PSV(1, NX, 1, NY, nt, (*wavePSV).pvx, (*wavePSV).pvxp1, (*wavePSV).pvxm1, (*wavePSV).pvy, (*wavePSV).pvyp1, (*wavePSV).pvym1, (*wavePSV).uttx, (*wavePSV).utty, (*wavePSV).psxx, (*wavePSV).psyy, 
                              (*wavePSV).psxy, (*matTTI).prip, (*matTTI).prjp, (*acq).srcpos_loc_back, (*seisPSVfwi).sectionvxdiff, (*seisPSVfwi).sectionvydiff,ntr,(*wavePSV_PML).absorb_coeff,hc,infoout, 1, (*wavePSV_PML).K_x,
 	                      (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, 
                              (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_sxx_x, (*wavePSV_PML).psi_syy_y, (*wavePSV_PML).psi_sxy_y, (*wavePSV_PML).psi_sxy_x);
              }
		                 
		/*if (MYID==0){
			time4=MPI_Wtime();
			time_av_v_update+=(time4-time3);
			if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
		}*/
		                                           
		/* exchange of particle velocities between PEs */
		exchange_v_PSV((*wavePSV).pvx,(*wavePSV).pvy, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);
		                                                       
		/*if (MYID==0){
		  time5=MPI_Wtime();
		  time_av_v_exchange+=(time5-time4);
		  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
		}*/                                                                                      	

	    /* update stress components */

	    /* viscoelastic case */
	    /* if (L)    
	    	update_s_visc_PML_PSV(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matPSV).ppi, (*matPSV).pu, 
                                  (*matPSV).puipjp, (*matPSV).prho, hc, infoout, (*wavePSV).pr, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).fipjp, (*matPSV).f, (*matPSV).g, (*matPSV).bip, (*matPSV).bjm, (*matPSV).cip, (*matPSV).cjm, 
                                  (*matPSV).d, (*matPSV).e, (*matPSV).dip, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, 
                                  (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx, 
                                  (*wavePSV_PML).psi_vyy, (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx, mode);
	    else*/

            /* elastic case */
	    update_s_elastic_PML_TTI(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matTTI).d11, (*matTTI).d13, 
                                     (*matTTI).d15, (*matTTI).d15h, (*matTTI).d33, (*matTTI).d35, (*matTTI).d35h, (*matTTI).d55h, (*wavePSV_PML).absorb_coeff, hc, infoout, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, 				     (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, 		
				     (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx, (*wavePSV_PML).psi_vyy, (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx, mode);  

	   /* explosive source */
	   if (QUELLTYP==1){
	   
	       if(mode==0 || mode==2){
	         psource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);
	       }
	       
           }

	   
	   /* adjoint explosive source */
           if((QUELLTYPB>=4)&&(mode==1)){ 	
	       psource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*acq).srcpos_loc_back,(*seisPSVfwi).sectionpdiff,nsrc_loc,1);
           }
 
	   /* moment tensor source */
	   if (QUELLTYP==5) 	
	   msource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);

	   if ((FREE_SURF) && (POS[2]==0)){

		/* viscoelastic */
	   	/* if (L)   
			surface_visc_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).ppi, 
					    (*matPSV).pu, (*matPSV).prho, (*matPSV).ptaup, (*matPSV).ptaus, (*matPSV).etajm, (*matPSV).peta, hc, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
					    (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs); */
		/* elastic */
		/* else      
	   		surface_elastic_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matPSV).ppi, (*matPSV).pu, (*matPSV).prho, hc, 
					       (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs); */
	   }

	   /*if (MYID==0){
	      time6=MPI_Wtime();
		  time_av_s_update+=(time6-time5);
	      if (infoout)  fprintf(FP," stress exchange between PEs ...");
	      }*/


	   /* stress exchange between PEs */
	    exchange_s_PSV((*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy, 
	      (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, 
	      (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top,
	      req_send, req_rec);

	   /*if (MYID==0){
	      time7=MPI_Wtime();
	 	  time_av_s_exchange+=(time7-time6);
	     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
	      }  */

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO && (mode==0 || mode==2)){
			seismo_ssg_VTI(nt, ntr, (*acq).recpos_loc, (*seisPSV).sectionvx, (*seisPSV).sectionvy, 
				(*seisPSV).sectionp, (*seisPSV).sectioncurl, (*seisPSV).sectiondiv, 
				(*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (nt==lsnap) && (nt<=iround(TSNAP2/DT))){

	      snap_VTI(FP,nt,++nsnap,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,hc);
	      lsnap=lsnap+iround(TSNAPINC/DT);

	   }

	      
	   /*if (MYID==0){
	      time8=MPI_Wtime();
		  time_av_timestep+=(time8-time3);
	      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
	      } */  		

	if((nt==hin1)&&(mode==0)&&(MODE>0)){

	    /* store forward wavefields for time-domain inversion and RTM */
            /* ---------------------------------------------------------- */
	    
		for (i=1;i<=NX;i=i+IDXI){
		    for (j=1;j<=NY;j=j+IDYI){
			 (*fwiPSV).forward_prop_rho_x[imat1]=(*wavePSV).pvxp1[j][i];
			 (*fwiPSV).forward_prop_rho_y[imat1]=(*wavePSV).pvyp1[j][i];
		         imat1++;                                   
		    }
		}   

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){
		    
			/* gradients with data integration */
		        if(GRAD_FORM==1){
			   (*fwiPSV).forward_prop_x[imat]=(*wavePSV).psxx[j][i];
			   (*fwiPSV).forward_prop_y[imat]=(*wavePSV).psyy[j][i];
		        }
		    
			/* gradients without data integration */
			if(GRAD_FORM==2){
		           (*fwiPSV).forward_prop_x[imat]=(*wavePSV).ux[j][i];
			   (*fwiPSV).forward_prop_y[imat]=(*wavePSV).uy[j][i];
			}
		    
			imat++;
		    }
		 }
	    

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){
		    
			/* gradients with data integration */
		        if(GRAD_FORM==1){
	 	          (*fwiPSV).forward_prop_u[imat2]=(*wavePSV).psxy[j][i];
			}

			/* gradients without data integration */
		        if(GRAD_FORM==2){
	 	          (*fwiPSV).forward_prop_u[imat2]=(*wavePSV).uxy[j][i];
		        }

			imat2++;
		    
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
		                                   
			   	(*fwiPSV).waveconv_rho_shot[j][i]+=((*wavePSV).pvxp1[j][i]*(*fwiPSV).forward_prop_rho_x[imat])+((*wavePSV).pvyp1[j][i]*(*fwiPSV).forward_prop_rho_y[imat]);
			
			   	/* mu-gradient with data integration */
			   	if(GRAD_FORM==1){
			
				   (*fwiPSV).waveconv_shot[j][i]+= ((*fwiPSV).forward_prop_x[imat]+(*fwiPSV).forward_prop_y[imat])*((*wavePSV).psxx[j][i]+(*wavePSV).psyy[j][i]);			  

		                   /*if(INVMAT1==1){
				       muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
				       lamss = (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] - 2.0 * muss;
				   }
			   
				   if(INVMAT1==3){
				       muss = (*matPSV).pu[j][i];
				       lamss = (*matPSV).ppi[j][i]; 
				   } 
			                
				   if(muss>0.0){
				      (*fwiPSV).waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*((*fwiPSV).forward_prop_u[imat] * (*wavePSV).psxy[j][i])) 
		                          + ((1.0/4.0) * (((*fwiPSV).forward_prop_x[imat] + (*fwiPSV).forward_prop_y[imat]) * ((*wavePSV).psxx[j][i] + (*wavePSV).psyy[j][i])) / ((lamss+muss)*(lamss+muss)))  
		                          + ((1.0/4.0) * (((*fwiPSV).forward_prop_x[imat] - (*fwiPSV).forward_prop_y[imat]) * ((*wavePSV).psxx[j][i] - (*wavePSV).psyy[j][i])) / (muss*muss));
				   }
				   */
				   
		                }

				/* Vs-gradient without data integration (stress-velocity in non-conservative form) */
		                if(GRAD_FORM==2){
			
				   (*fwiPSV).waveconv_shot[j][i]+= ((*fwiPSV).forward_prop_x[imat]+(*fwiPSV).forward_prop_y[imat])*((*wavePSV).psxx[j][i]+(*wavePSV).psyy[j][i]);

		                   /*if(INVMAT1==1){
				       muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
				       lamss = (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] - 2.0 * muss;
				   }
			   
				   if(INVMAT1==3){
				       muss = (*matPSV).pu[j][i];
				       lamss = (*matPSV).ppi[j][i]; 
				   } 
			                
				   if(muss>0.0){
								
				      tmp = (1.0/(4.0*(lamss+muss)*(lamss+muss))) - (1.0/(4.0*muss*muss));
				      tmp1 = (1.0/(4.0*(lamss+muss)*(lamss+muss))) + (1.0/(4.0*muss*muss));
			
				      (*fwiPSV).waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*((*fwiPSV).forward_prop_u[imat] * (*wavePSV).psxy[j][i])) 
		                                           + ( tmp1 * ((*fwiPSV).forward_prop_x[imat] * (*wavePSV).psxx[j][i] + (*fwiPSV).forward_prop_y[imat] * (*wavePSV).psyy[j][i]))  
		                                           + ( tmp  * ((*fwiPSV).forward_prop_x[imat] * (*wavePSV).psyy[j][i] + (*fwiPSV).forward_prop_y[imat] * (*wavePSV).psxx[j][i]));	  
					  
				   } */
		                  
		                }			
						                                                                                                     
			   imat++;
			   }
		    }  
		  
		  if(EPRECOND==1){
		     eprecond(Wr,(*wavePSV).pvx,(*wavePSV).pvy);
		  }
		                                                                                                                       
	    hin++;
	    }

	   }/*--------------------  End  of loop over timesteps ----------*/		

}
