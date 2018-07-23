/*------------------------------------------------------------------------
 *  Output of time reversed residual seismograms (PSV problem)
 * 
 *  Daniel Koehn
 *  Kiel, 25/04/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void outseis_PSVres(struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int nstage, FILE *FP){ 
		
        /* global variables */
	extern int QUELLTYPB, MYID_SHOT;
        extern MPI_Comm SHOT_COMM;	
	
        /* local variables */

        if(QUELLTYPB==1){
   
           catseis((*seisPSVfwi).sectionvxdiff, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
           catseis((*seisPSVfwi).sectionvydiff, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
      
           if (MYID_SHOT==0){
	      saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).sectionp,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
           }
      
        }
   
        if(QUELLTYPB==2){
   
           catseis((*seisPSVfwi).sectionvydiff, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
      
           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).fulldata_vy,(*seisPSV).fulldata_vy,(*seisPSV).sectionvy,(*seisPSV).sectionvy,(*seisPSV).sectionvy,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
           }
      
        }
   
        if(QUELLTYPB==3){
   
           catseis((*seisPSVfwi).sectionvxdiff, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
      
           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).sectionp,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
           }
      
        }
   
        if(QUELLTYPB==4){
   
           catseis((*seisPSVfwi).sectionpdiff, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);
      
           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).sectionvx,(*seisPSV).sectionvy,(*seisPSV).fulldata_p,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
           }
      
        }

        if(QUELLTYPB==5){

           catseis((*seisPSVfwi).sectionvxdiff, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
           catseis((*seisPSVfwi).sectionpdiff, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);
 
           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).fulldata_p,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage);
           }

        }

        if(QUELLTYPB==6){

           catseis((*seisPSVfwi).sectionvydiff, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
           catseis((*seisPSVfwi).sectionpdiff, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);

           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).fulldata_p,(*seisPSV).sectionvy,(*seisPSV).sectionvy,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage);
           }

        }

        if(QUELLTYPB==7){

           catseis((*seisPSVfwi).sectionvxdiff, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
           catseis((*seisPSVfwi).sectionvydiff, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
           catseis((*seisPSVfwi).sectionpdiff, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);
           if (MYID_SHOT==0){
              saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).fulldata_p,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage);
           }

        } 
		
}
