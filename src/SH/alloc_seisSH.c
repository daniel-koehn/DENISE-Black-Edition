/*
 * Allocate memory for seismograms (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 03/12/2017
 */

#include "fd.h"

void alloc_seisSH(int ntr, int ns, struct seisSH *seisSH){

        /* global variables */
	extern int SEISMO;

	/* local variables */	

   if (ntr>0){
           switch (SEISMO){
           case 1 : /* particle velocities only */
                   (*seisSH).sectionvz=matrix(1,ntr,1,ns);                   
                   break;           
           }
   }


}



