/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D acoustic forward problem by finite-differences 
 *   for a single shot 
 *
 *   mode = 0 - forward modelling only, STF estimation or FWI gradient calculation
 *   mode = 1 - backpropagation of data residuals
 *   mode = 2 - evaluation objective function for step length estimation  
 * 
 *   
 *   D. Koehn
 *   Kiel, 10.06.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void ac(struct waveAC *waveAC, struct waveAC_PML *waveAC_PML, struct matAC *matAC, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, int hin, int *DTINV_help, int mode, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, MYID_SHOT, FDORDER, FW, L, GRAD_FORM, FC_SPIKE_1, FC_SPIKE_2, ORDER_SPIKE;
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
			    
	/* initialize AC wavefields with zero */
	/*if (L){
		zero_denise_visc_PSV(-nd+1,NY+nd,-nd+1,NX+nd, (*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                                 (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,(*wavePSV_PML).psi_sxy_x,
                                 (*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,(*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,
                                 (*wavePSV_PML).psi_vxxs,(*wavePSV).pr,(*wavePSV).pp,(*wavePSV).pq);
	}else{*/	
		zero_denise_acoustic_AC(-nd+1,NY+nd,-nd+1,NX+nd,(*waveAC).pvx,(*waveAC).pvy,(*waveAC).p,
                            (*waveAC).ux,(*waveAC).pvxp1, (*waveAC).pvyp1,(*waveAC_PML).psi_p_x,
                            (*waveAC_PML).psi_vxx,(*waveAC_PML).psi_p_y,(*waveAC_PML).psi_vyy,(*waveAC_PML).psi_vxxs);	
	/*} */                                                        
	     
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
		if (isnan((*waveAC).pvy[NY/2][NX/2])) {
		   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,(*waveAC).pvy[NY/2][NX/2]);
		   err(" Simulation is unstable !");}

	   infoout = !(nt%10000);

	   if (MYID==0){
	      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
	      /*time3=MPI_Wtime();*/
	   }

	      /* update of particle velocities */
              if(mode==0 || mode==2){
	         update_v_PML_AC(1, NX, 1, NY, nt, (*waveAC).pvx, (*waveAC).pvxp1, (*waveAC).pvxm1, (*waveAC).pvy, (*waveAC).pvyp1, (*waveAC).pvym1, (*waveAC).uttx, (*waveAC).utty, (*waveAC).p,       
                              (*matAC).prip, (*matAC).prjp, (*acq).srcpos_loc,(*acq).signals,(*acq).signals,nsrc_loc,(*waveAC_PML).absorb_coeff,hc,infoout, 0, (*waveAC_PML).K_x, (*waveAC_PML).a_x, 
                              (*waveAC_PML).b_x, (*waveAC_PML).K_x_half, (*waveAC_PML).a_x_half, (*waveAC_PML).b_x_half, (*waveAC_PML).K_y, (*waveAC_PML).a_y, (*waveAC_PML).b_y, (*waveAC_PML).K_y_half, 
                              (*waveAC_PML).a_y_half, (*waveAC_PML).b_y_half, (*waveAC_PML).psi_p_x, (*waveAC_PML).psi_p_y);
              }

              if(mode==1){
	         update_v_PML_AC(1, NX, 1, NY, nt, (*waveAC).pvx, (*waveAC).pvxp1, (*waveAC).pvxm1, (*waveAC).pvy, (*waveAC).pvyp1, (*waveAC).pvym1, (*waveAC).uttx, (*waveAC).utty, (*waveAC).p, 
                              (*matAC).prip, (*matAC).prjp, (*acq).srcpos_loc_back, (*seisPSVfwi).sectionvxdiff, (*seisPSVfwi).sectionvydiff,ntr,(*waveAC_PML).absorb_coeff,hc,infoout, 1, (*waveAC_PML).K_x,
 	                      (*waveAC_PML).a_x, (*waveAC_PML).b_x, (*waveAC_PML).K_x_half, (*waveAC_PML).a_x_half, (*waveAC_PML).b_x_half, (*waveAC_PML).K_y, (*waveAC_PML).a_y, (*waveAC_PML).b_y, (*waveAC_PML).K_y_half, 
                              (*waveAC_PML).a_y_half, (*waveAC_PML).b_y_half, (*waveAC_PML).psi_p_x, (*waveAC_PML).psi_p_y);
              }
		                 
		/*if (MYID==0){
			time4=MPI_Wtime();
			time_av_v_update+=(time4-time3);
			if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
		}*/
		                                           
		/* exchange of particle velocities between PEs */
		exchange_v_AC((*waveAC).pvx,(*waveAC).pvy, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);
		                                                       
		/*if (MYID==0){
		  time5=MPI_Wtime();
		  time_av_v_exchange+=(time5-time4);
		  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
		}*/                                                                                      	

	    /*if (L)  */  /* viscoelastic */
	    	/*update_s_visc_PML_PSV(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matPSV).ppi, (*matPSV).pu, 
                                  (*matPSV).puipjp, (*matPSV).prho, hc, infoout, (*wavePSV).pr, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).fipjp, (*matPSV).f, (*matPSV).g, (*matPSV).bip, (*matPSV).bjm, (*matPSV).cip, (*matPSV).cjm, 
                                  (*matPSV).d, (*matPSV).e, (*matPSV).dip, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, 
                                  (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx, 
                                  (*wavePSV_PML).psi_vyy, (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx, mode);
	    else*/

	   	update_s_acoustic_PML_AC(1, NX, 1, NY, (*waveAC).pvx, (*waveAC).pvy, (*waveAC).ux, (*waveAC).p, (*matAC).ppi, 
                                     (*waveAC_PML).absorb_coeff, (*matAC).prho, hc, infoout, (*waveAC_PML).K_x, (*waveAC_PML).a_x, (*waveAC_PML).b_x, (*waveAC_PML).K_x_half, (*waveAC_PML).a_x_half, 
                                     (*waveAC_PML).b_x_half, (*waveAC_PML).K_y, (*waveAC_PML).a_y, (*waveAC_PML).b_y, (*waveAC_PML).K_y_half, (*waveAC_PML).a_y_half, (*waveAC_PML).b_y_half, (*waveAC_PML).psi_vxx,  
                                     (*waveAC_PML).psi_vyy, mode);  


	   /* explosive source */
	   if (QUELLTYP==1){
	   
	       if(mode==0 || mode==2){
	         psource_AC(nt,(*waveAC).p,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);
	       }
	       
           }

	   
	   /* adjoint explosive source */
           if((QUELLTYPB>=4)&&(mode==1)){ 	
	       psource_AC(nt,(*waveAC).p,(*acq).srcpos_loc_back,(*seisPSVfwi).sectionpdiff,ntr,1);
           }
 	 
	   if ((FREE_SURF) && (POS[2]==0)){
	   	/* if (L) */   /* visco-acoustic */
			/*surface_visc_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).ppi, 
					    (*matPSV).pu, (*matPSV).prho, (*matPSV).ptaup, (*matPSV).ptaus, (*matPSV).etajm, (*matPSV).peta, hc, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
					    (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
		else */     /* acoustic */
	   		
                surface_acoustic_PML_AC(1, (*waveAC).p);
	   }

	   /*if (MYID==0){
	      time6=MPI_Wtime();
		  time_av_s_update+=(time6-time5);
	      if (infoout)  fprintf(FP," stress exchange between PEs ...");
	      }*/


	   /* stress exchange between PEs */
	   exchange_p_AC((*waveAC).p, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);

	   /*if (MYID==0){
	      time7=MPI_Wtime();
	 	  time_av_s_exchange+=(time7-time6);
	     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
	      }  */

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO && (mode==0 || mode==2)){
			seismo_AC(nt, ntr, (*acq).recpos_loc, (*seisPSV).sectionvx, (*seisPSV).sectionvy, 
				(*seisPSV).sectionp, (*seisPSV).sectioncurl, (*seisPSV).sectiondiv, 
				(*waveAC).pvx, (*waveAC).pvy, (*waveAC).p, (*matAC).ppi, (*matAC).ppi, (*matAC).prho, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (nt==lsnap) && (nt<=iround(TSNAP2/DT))){

	      snap_AC(FP,nt,++nsnap,(*waveAC).pvx,(*waveAC).pvy,(*waveAC).p,(*matAC).ppi,(*matAC).ppi,hc);

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
			 (*fwiPSV).forward_prop_rho_x[imat1]=(*waveAC).pvxp1[j][i];
			 (*fwiPSV).forward_prop_rho_y[imat1]=(*waveAC).pvyp1[j][i];
		         imat1++;                                   
		    }
		}   

		for (i=1;i<=NX;i=i+IDXI){ 
		    for (j=1;j<=NY;j=j+IDYI){
		    
			/* gradients with data integration */
		        if(GRAD_FORM==1){
			   (*fwiPSV).forward_prop_x[imat]=(*waveAC).p[j][i];
		        }
			
			/* gradients without data integration */
		        if(GRAD_FORM==2){
			   (*fwiPSV).forward_prop_x[imat]=(*waveAC).ux[j][i];
		        }
					    		    
			imat++;
		    }
		 }
	    
	   
	    if((EPRECOND==1)||(EPRECOND==3)){
	      eprecond(Ws,(*waveAC).pvx,(*waveAC).pvy);
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
		                
				if(GRAD_FORM==1){
				   (*fwiPSV).waveconv_shot[j][i] += (*fwiPSV).forward_prop_x[imat] * (*waveAC).p[j][i];                   
			   	   (*fwiPSV).waveconv_rho_shot[j][i] += ((*waveAC).pvxp1[j][i]*(*fwiPSV).forward_prop_rho_x[imat])+((*waveAC).pvyp1[j][i]*(*fwiPSV).forward_prop_rho_y[imat]);
				}
				
				if(GRAD_FORM==2){
				   (*fwiPSV).waveconv_shot[j][i] += (*fwiPSV).forward_prop_x[imat] * (*waveAC).p[j][i];                   
			   	   (*fwiPSV).waveconv_rho_shot[j][i] += ((*waveAC).pvxp1[j][i]*(*fwiPSV).forward_prop_rho_x[imat])+((*waveAC).pvyp1[j][i]*(*fwiPSV).forward_prop_rho_y[imat]);
				}
				   
																                                                                                                     
			   imat++;
			   }
		    }  
		  
		  if(EPRECOND==1){
		     eprecond(Wr,(*waveAC).pvx,(*waveAC).pvy);
		  }
		                                                                                                                       
	    hin++;
	    }

	   }/*--------------------  End  of loop over timesteps ----------*/		

}
