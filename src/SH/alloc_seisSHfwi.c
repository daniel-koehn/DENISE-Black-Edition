/*
 * Allocate memory for seismograms used for FWI (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

void alloc_seisSHfwi(int ntr, int ntr_glob, int ns, struct seisSHfwi *seisSHfwi){

        /* global variables */
	extern int QUELLTYPB;

	/* local variables */	

   	(*seisSHfwi).sectionread=matrix(1,ntr_glob,1,ns);
   
	/* FWI of z-component */
   	if(QUELLTYPB){
      	   (*seisSHfwi).sectionvzdata=matrix(1,ntr,1,ns);
      	   (*seisSHfwi).sectionvzdiff=matrix(1,ntr,1,ns);
      	   (*seisSHfwi).sectionvzdiffold=matrix(1,ntr,1,ns);
   	}   

}



