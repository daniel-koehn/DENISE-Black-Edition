/*------------------------------------------------------------------------
 *  Calculate data residuals and objective function (PSV problem)
 * 
 *  Daniel Koehn
 *  Kiel, 25/04/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_res_PSV(struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob,  int ntr, int nsrc_glob, float ** srcpos, int ishot, int ns, int iter, int swstestshot){ 
		
/* global variables */
extern int QUELLTYPB, TIMELAPSE, TIME_FILT, ORDER, LNORM, MYID, COMP_WEIGHT;
extern float FC, FC_START;	
	
/* local variables */
int i, j, h;

if (MYID==0){
   printf("-------------------  \n");
   printf("Calculate residuals  \n");
   printf("-------------------  \n");
}

/* read seismic data from SU file vx */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==3)||(QUELLTYPB==5)||(QUELLTYPB==7)){ /* if QUELLTYPB */

inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,1,iter);

if (TIME_FILT){
   apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           (*seisPSVfwi).sectionvxdata[h][j]=(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
  /* read synthetic seismic data at time step t0 vx */
  inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,9,iter);

  if (TIME_FILT){
   apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
  }
  
  /* calculate vx_mod(t1) - vx_mod(t0) */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           (*seisPSV).sectionvx[h][j]=(*seisPSV).sectionvx[h][j]-(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
   }
   h++;
}
                      
} /* end of TIMELAPSE */

  COMP_WEIGHT = 0;
  if((QUELLTYPB==5)||(QUELLTYPB==7)){
      COMP_WEIGHT = 1;
  }

  (*seisPSVfwi).L2=calc_res((*seisPSVfwi).sectionvxdata,(*seisPSV).sectionvx,(*seisPSVfwi).sectionvxdiff,(*seisPSVfwi).sectionvxdiffold,ntr,ns,LNORM,(*seisPSVfwi).L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);     

  if(swstestshot==1){(*seisPSVfwi).energy=calc_energy((*seisPSVfwi).sectionvxdata,ntr,ns,(*seisPSVfwi).energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}
  /*L2_all_shots=calc_misfit((*seisPSVfwi).sectionvxdiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
  energy_all_shots=calc_energy((*seisPSVfwi).sectionvxdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);*/

} /* end QUELLTYPB */

/* read seismic data from SU file vy */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==2)||(QUELLTYPB==6)||(QUELLTYPB==7)){ /* if QUELLTYPB */

inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,2,iter);

if (TIME_FILT){
   apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           (*seisPSVfwi).sectionvydata[h][j]=(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
    /* read synthetic seismic data at time step t0 vy */
    inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,10,iter); 

    if (TIME_FILT){	
       apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
    }
   
    /* calculate vy_mod(t1) - vy_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
       for(j=1;j<=ns;j++){
	     (*seisPSV).sectionvy[h][j]=(*seisPSV).sectionvy[h][j]-(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
       }
       h++;
    }
                               
} /* end of TIMELAPSE */

  COMP_WEIGHT = 0;
  if((QUELLTYPB==6)||(QUELLTYPB==7)){
      COMP_WEIGHT = 2;
  }
                               
(*seisPSVfwi).L2=calc_res((*seisPSVfwi).sectionvydata,(*seisPSV).sectionvy,(*seisPSVfwi).sectionvydiff,(*seisPSVfwi).sectionvydiffold,ntr,ns,LNORM,(*seisPSVfwi).L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);

if(swstestshot==1){(*seisPSVfwi).energy=calc_energy((*seisPSVfwi).sectionvydata,ntr,ns,(*seisPSVfwi).energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}
/*L2_all_shots=calc_misfit((*seisPSVfwi).sectionvydiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
energy_all_shots=calc_energy((*seisPSVfwi).sectionvydata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);	*/   	    


} /* end QUELLTYPB */

/* read seismic data from SU file p */
/* --------------------------------- */
if(QUELLTYPB>=4){ /* if QUELLTYPB */

inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,11,iter);

if (TIME_FILT){
   apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
}

h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           (*seisPSVfwi).sectionpdata[h][j]=(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){

    /* read synthetic seismic data at time step t0 vy */
    inseis(ishot,(*seisPSVfwi).sectionread,ntr_glob,ns,15,iter);

    if (TIME_FILT){
       apply_tdfilt((*seisPSVfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
    }

    /* calculate p_mod(t1) - p_mod(t0) */
    h=1;
    for(i=1;i<=ntr;i++){
       for(j=1;j<=ns;j++){
             (*seisPSV).sectionp[h][j]=(*seisPSV).sectionp[h][j]-(*seisPSVfwi).sectionread[recpos_loc[3][i]][j];
       }
       h++;
    }

} /* end of TIMELAPSE */

COMP_WEIGHT = 0;

(*seisPSVfwi).L2=calc_res((*seisPSVfwi).sectionpdata,(*seisPSV).sectionp,(*seisPSVfwi).sectionpdiff,(*seisPSVfwi).sectionpdiffold,ntr,ns,LNORM,(*seisPSVfwi).L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);

if(swstestshot==1){(*seisPSVfwi).energy=calc_energy((*seisPSVfwi).sectionpdata,ntr,ns,(*seisPSVfwi).energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}

/*L2_all_shots=calc_misfit((*seisPSVfwi).sectionpdiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
energy_all_shots=calc_energy((*seisPSVfwi).sectionpdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);*/


} /* end QUELLTYPB >= 4*/			

}
