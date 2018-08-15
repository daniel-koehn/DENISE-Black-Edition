/*------------------------------------------------------------------------
 *  Calculate data residuals and objective function (SH problem)
 * 
 *  Daniel Koehn
 *  Kiel, 13/12/2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_res_SH(struct seisSH *seisSH, struct seisSHfwi *seisSHfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob,  int ntr, int nsrc_glob, float ** srcpos, int ishot, int ns, int iter, int swstestshot){ 
		
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

/* read seismic data from SU file vz */
/* --------------------------------- */
if(QUELLTYPB){ /* if QUELLTYPB */

    inseis(ishot,(*seisSHfwi).sectionread,ntr_glob,ns,2,iter);

    if (TIME_FILT){
        apply_tdfilt((*seisSHfwi).sectionread,ntr_glob,ns,ORDER,FC,FC_START);
    }

    h=1;
    for(i=1;i<=ntr;i++){
         for(j=1;j<=ns;j++){
             (*seisSHfwi).sectionvzdata[h][j]=(*seisSHfwi).sectionread[recpos_loc[3][i]][j];
         }
         h++;
    }

    COMP_WEIGHT = 0;  

    (*seisSHfwi).L2=calc_res((*seisSHfwi).sectionvzdata,(*seisSH).sectionvz,(*seisSHfwi).sectionvzdiff,(*seisSHfwi).sectionvzdiffold,ntr,ns,LNORM,(*seisSHfwi).L2,0,1,swstestshot,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);     

    if(swstestshot==1){(*seisSHfwi).energy=calc_energy((*seisSHfwi).sectionvzdata,ntr,ns,(*seisSHfwi).energy, ntr_glob, recpos_loc, nsrc_glob, ishot);}
    /*L2_all_shots=calc_misfit((*seisPSVfwi).sectionvxdiff,ntr,ns,LNORM,L2_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);
    energy_all_shots=calc_energy((*seisPSVfwi).sectionvxdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot);*/

} /* end QUELLTYPB */


}
