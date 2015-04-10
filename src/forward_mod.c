/*------------------------------------------------------------------------
 * Pure forward modelling code
 *
 * D. Koehn
 * Kiel, 3.2.2012
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void forward_mod(FILE *fprec, float ** waveconv, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, float ** ppi, float ** ppinp1, int iter, float eps_scale, int nfstart,
	int nsrc, float ** puipjp, float ** prip, float ** prjp, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvx, float ** pvy, float ** psxx, float ** psyy, float ** psxy, float ** ux, float ** uy, float ** pvxp1, float ** pvyp1, float ** psi_sxx_x, float ** psi_sxy_x,
	float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ** pvxm1, float ** pvym1, float ** uttx, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uxy, float ** uyx, int ntr, int **recpos_loc, float **sectionvx, float **sectionvy, float **sectionp, float **sectioncurl, 
	float **sectiondiv, float **sectionread, int ntr_glob, float ** sectionvxdata, float ** sectionvxdiff, float ** sectionvxdiffold, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, float ** sectionpdata, float ** sectionpdiff, float ** sectionpdiffold, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, 
        float ** bufferrig_to_lef, float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, int itest, int nsrc_glob, int nsrc_loc, MPI_Request * req_send, MPI_Request * req_rec, float ***pr, 
        float ***pp, float ***pq, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, float **ptaup, float **ptaus, 
        float *etajm, float *peta, float *etaip, float **ptausipjp, int **recpos, float FC, int * recswitch, FILE *FP, int ntr_loc){

	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY, FDORDER, INVMAT, MYID, QUELLART, QUELLTYP, SEISMO, QUELLTYPB, LNORM;
	extern int RUN_MULTIPLE_SHOTS, TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, NT, NDT, FREE_SURF, TRKILL, REC1, REC2;
	extern int FREQFILT, L, ORDER_SPIKE, TIME_FILT, INV_STF, ORDER, CHECKPTREAD, TIMELAPSE, ENV, N_STREAMER;
	
	extern float TSNAP1, DT, FC_SPIKE_1, FC_SPIKE_2, FC_START;

	int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;
	float alphanom, alphadenom, tmp;


/* calculate change in the material parameters */
tmp=calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,INVMAT,eps_scale,1,nfstart);

/*char modfile[STRING_SIZE];

sprintf(modfile,"%s_vp_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,ppinp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_countstep%d.bin",INV_MODELFILE,countstep);

writemod(modfile,punp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_rho_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,prhonp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);*/


/*if(MYID==0){printf("EPSILON = %e \n",EPSILON);}*/

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prhonp1, ppinp1, punp1);	/* no differentiation of elastic and viscoelastic modelling because the 
						viscoelastic parameters did not change during the forward modelling */
MPI_Barrier(MPI_COMM_WORLD);

av_mue(punp1,puipjp,prhonp1);
av_rho(prhonp1,prip,prjp);


/* Preparing memory variables for update_s (viscoelastic) */
if (L) prepare_update_s(etajm,etaip,peta,fipjp,punp1,puipjp,ppinp1,prhonp1,ptaus,ptaup,ptausipjp,f,g,
		bip,bjm,cip,cjm,dip,d,e);
		
/* initialization of L2 calculation */
L2=0.0;

alphanom = 0.0;
alphadenom = 0.0;

exchange_par();
 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
    
        for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){
 
           if(N_STREAMER>0){ 
	   
              if (SEISMO){
                  recpos=receiver(FP, &ntr, ishot);
                  recswitch = ivector(1,ntr);
                  recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
                  ntr_glob=ntr;
                  ntr=ntr_loc;
              }

              if (ntr>0){

                 switch (SEISMO){
                 case 1 : /* particle velocities only */
                        sectionvx=matrix(1,ntr,1,ns);
                        sectionvy=matrix(1,ntr,1,ns);
                        break;
                 case 2 : /* pressure only */
                        sectionp=matrix(1,ntr,1,ns);
                        break;
                 case 3 : /* curl and div only */
                        sectioncurl=matrix(1,ntr,1,ns);
                        sectiondiv=matrix(1,ntr,1,ns);
                        break;
                 case 4 : /* everything */
                        sectionvx=matrix(1,ntr,1,ns);
                        sectionvy=matrix(1,ntr,1,ns);
                        sectioncurl=matrix(1,ntr,1,ns);
                        sectiondiv=matrix(1,ntr,1,ns);
                        sectionp=matrix(1,ntr,1,ns);
                        break;
                 }

              }

              /* Memory for seismic data */
              sectionread=matrix(1,ntr_glob,1,ns);

              if((QUELLTYPB==1)||(QUELLTYPB==3)){
                sectionvxdata=matrix(1,ntr,1,ns);
                sectionvxdiff=matrix(1,ntr,1,ns);
                sectionvxdiffold=matrix(1,ntr,1,ns); 
              }
              
              if((QUELLTYPB==1)||(QUELLTYPB==2)){
                sectionvydata=matrix(1,ntr,1,ns);
                sectionvydiff=matrix(1,ntr,1,ns);
                sectionvydiffold=matrix(1,ntr,1,ns);
	      }

              if(QUELLTYPB==4){
                sectionpdata=matrix(1,ntr,1,ns);
                sectionpdiff=matrix(1,ntr,1,ns);
                sectionpdiffold=matrix(1,ntr,1,ns);
              }

          }       
 
        if(MYID==0){
          printf("\n=================================================================================================\n");
          printf("\n *****  Starting simulation (test-forward model) no. %d for shot %d of %d (rel. step length %.8f) \n",itest,ishot,nshots,eps_scale);
	  printf("\n=================================================================================================\n\n");
	}
		
        for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
        if (RUN_MULTIPLE_SHOTS){

	    /* find this single source positions on subdomains */
           if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
           srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		}
		           
		else{
		/* Distribute multiple source positions on subdomains */
		   srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
		}

/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc,ishot);

if (nsrc_loc){if(QUELLART==6){

   /* low pass filtering of spike */
   timedomain_filt(signals,FC_SPIKE_2,ORDER_SPIKE,nsrc_loc,ns,1);

   /* band-pass filtering of spike */
   if(FC_SPIKE_1 > 0.0){
   timedomain_filt(signals,FC_SPIKE_1,ORDER_SPIKE,nsrc_loc,ns,2);
   }

  }

}

/* calculate envelope (Chi, Dong & Liu, 2014) */      
if (ENV==1){ 

   calc_envelope(signals,signals,ns,nsrc_loc); 

} 

if((TIME_FILT)&&(INVMAT!=10)&&(INV_STF==0)){

  /*time domain low-pass filtering of the source signal */
  timedomain_filt(signals,FC,ORDER,nsrc_loc,ns,1);

  if(TIME_FILT==2){ /* band-pass filter */
    timedomain_filt(signals,FC_START,ORDER,nsrc_loc,ns,2);
  }

}
		    
/* initialize wavefield with zero */
if (L){
	zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pr,pp,pq);
}else{	
	zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	
}     

/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");

		
   infoout = !(nt%10000);



   /* update of particle velocities */
   update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);


   /* exchange of particle velocities between PEs */
   exchange_v(pvx, pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);


   if (L)    /* viscoelastic */
    	update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, prhonp1, hc, infoout,
				pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
				K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
    else
   	update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, absorb_coeff, prhonp1, hc, infoout,
         	               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  

 
    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc,0);

   if ((FREE_SURF) && (POS[2]==0)){
   	if (L)    /* viscoelastic */
		surface_PML(1, pvx, pvy, psxx, psyy, psxy, pp, pq, ppinp1, punp1, prhonp1, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs);
	else      /* elastic */
   		surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppinp1, punp1, prhonp1, hc, K_x, a_x, b_x, psi_vxxs);
   }

   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec); 

	/* store amplitudes at receivers in section-arrays */
	if (SEISMO){
		seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppinp1, punp1, hc);
		/*lsamp+=NDT;*/
	}
		
   }/*--------------------  End  of loop over timesteps (test forward model) ----------*/


if (MYID==0){
printf("Calculate residuals between test forward model m - mu * dm and actual model m \n");
printf("----------------------------------------------------------------------------- \n");
}

if (ntr > 0){

/* read seismic data from SU file vx */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,1,iter);

/* calculate envelope (Chi, Dong & Liu, 2014) */      
if (ENV==1){ 

   calc_envelope(sectionread,sectionread,ns,ntr_glob); 
   calc_envelope(sectionvx,sectionvx,ns,ntr);  

} 
			
if (TIME_FILT){

   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}


/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
    /* read synthetic seismic data at time step t0 vx */
    inseis(fprec,ishot,sectionread,ntr_glob,ns,9,iter);

    if (TIME_FILT){

       timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

       if(TIME_FILT==2){ /* band-pass */
         timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
       }

    }
       
    /* calculate vx_mod(t1) - vx_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
         for(j=1;j<=ns;j++){   
	        sectionvx[h][j]=sectionvx[h][j]-sectionread[recpos_loc[3][i]][j];
	 }
	 h++;
    } 
                              
} /* end of TIMELAPSE */

L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);

} /* end QUELLTYPB*/

/* read seismic data from SU file vy */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);

/* calculate envelope (Chi, Dong & Liu, 2014) */      
if (ENV==1){ 

   calc_envelope(sectionread,sectionread,ns,ntr_glob); 
   calc_envelope(sectionvy,sectionvy,ns,ntr);  

} 

if (TIME_FILT){
   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}

/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
      /* read synthetic seismic data at time step t0 vy */
      inseis(fprec,ishot,sectionread,ntr_glob,ns,10,iter);

      if (TIME_FILT){
	timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

        if(TIME_FILT==2){ /* band-pass */
          timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
        }
      } 
            
      /* calculate vy_mod(t1) - vy_mod(t0) */
      h=1;
      for(i=1;i<=ntr;i++){
	 for(j=1;j<=ns;j++){
	      sectionvy[h][j]=sectionvy[h][j]-sectionread[recpos_loc[3][i]][j];
	 }
	 h++;
      }
                                                                                                   
} /* end of TIMELAPSE */

L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);			   	    

} /* end QUELLTYPB */

/* read seismic data from SU file p */
/* --------------------------------- */

if(QUELLTYPB==4){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,11,iter);

/* calculate envelope (Chi, Dong & Liu, 2014) */      
if (ENV==1){ 

   calc_envelope(sectionread,sectionread,ns,ntr_glob); 
   calc_envelope(sectionvx,sectionvx,ns,ntr);  

} 
			
if (TIME_FILT){

   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }

}


/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionpdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate p_mod(t1) - p_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
    /* read synthetic seismic data at time step t0 vx */
    inseis(fprec,ishot,sectionread,ntr_glob,ns,15,iter);

    if (TIME_FILT){

       timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

       if(TIME_FILT==2){ /* band-pass */
         timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
       }

    }
       
    /* calculate vx_mod(t1) - vx_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
         for(j=1;j<=ns;j++){   
	        sectionp[h][j]=sectionp[h][j]-sectionread[recpos_loc[3][i]][j];
	 }
	 h++;
    } 
                              
} /* end of TIMELAPSE */

L2=calc_res(sectionpdata,sectionp,sectionpdiff,sectionpdiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);

} /* end QUELLTYPB*/


}

if(N_STREAMER>0){

	if (SEISMO) free_imatrix(recpos,1,3,1,ntr_glob);

	if ((ntr>0) && (SEISMO)){

        	free_imatrix(recpos_loc,1,3,1,ntr);

        	switch (SEISMO){
        	case 1 : /* particle velocities only */
                	free_matrix(sectionvx,1,ntr,1,ns);
              		free_matrix(sectionvy,1,ntr,1,ns);
                	break;
        	case 2 : /* pressure only */
                	free_matrix(sectionp,1,ntr,1,ns);
                	break;
        	case 3 : /* curl and div only */
                	free_matrix(sectioncurl,1,ntr,1,ns);
                	free_matrix(sectiondiv,1,ntr,1,ns);
                	break;
        	case 4 : /* everything */
                	free_matrix(sectionvx,1,ntr,1,ns);
                	free_matrix(sectionvy,1,ntr,1,ns);
                	free_matrix(sectionp,1,ntr,1,ns);
                	free_matrix(sectioncurl,1,ntr,1,ns);
                	free_matrix(sectiondiv,1,ntr,1,ns);
                	break;

        	}

	}
	

 	free_ivector(recswitch,1,ntr);
 	free_matrix(sectionread,1,ntr_glob,1,ns);

        if((QUELLTYPB==1)||(QUELLTYPB==3)){
 	  free_matrix(sectionvxdata,1,ntr,1,ns);
 	  free_matrix(sectionvxdiff,1,ntr,1,ns);
          free_matrix(sectionvxdiffold,1,ntr,1,ns);
        }

        if((QUELLTYPB==1)||(QUELLTYPB==2)){
 	  free_matrix(sectionvydata,1,ntr,1,ns);
 	  free_matrix(sectionvydiff,1,ntr,1,ns);
 	  free_matrix(sectionvydiffold,1,ntr,1,ns);
        }

        if(QUELLTYPB==4){
          free_matrix(sectionvydata,1,ntr,1,ns);
          free_matrix(sectionvydiff,1,ntr,1,ns);
          free_matrix(sectionvydiffold,1,ntr,1,ns);
        }

}

 nsrc_loc=0;


} /* ===========================================================================================================================*/   
/* ==================================== end of loop over shots (test forward) ==================================================*/
/* =============================================================================================================================*/
epst1[itest]=eps_scale;
epst1[1] = 0.0;    

L2sum=0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
if(LNORM==2){   
  L2t[itest] = L2sum/energy_sum;
}
else {L2t[itest] = L2sum;}
      
}

