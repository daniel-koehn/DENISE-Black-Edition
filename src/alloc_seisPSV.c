/*
 * Allocate memory for seismograms (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void alloc_seisPSV(int ntr, int ns, struct seisPSV *seisPSV){

        /* global variables */
	extern int SEISMO;

	/* local variables */	

   if (ntr>0){
           switch (SEISMO){
           case 1 : /* particle velocities only */
                   (*seisPSV).sectionvx=matrix(1,ntr,1,ns);
                   (*seisPSV).sectionvy=matrix(1,ntr,1,ns);
                   break;
           case 2 : /* pressure only */
                   (*seisPSV).sectionp=matrix(1,ntr,1,ns);
                   break;
           case 3 : /* curl and div only */
                   (*seisPSV).sectioncurl=matrix(1,ntr,1,ns);
                   (*seisPSV).sectiondiv=matrix(1,ntr,1,ns);
                   break;
           case 4 : /* everything */
                   (*seisPSV).sectionvx=matrix(1,ntr,1,ns);
                   (*seisPSV).sectionvy=matrix(1,ntr,1,ns);
                   (*seisPSV).sectioncurl=matrix(1,ntr,1,ns);
                   (*seisPSV).sectiondiv=matrix(1,ntr,1,ns);
                   (*seisPSV).sectionp=matrix(1,ntr,1,ns);
                   break;
           case 5: /* pressure and particle velocities only */
                   (*seisPSV).sectionvx=matrix(1,ntr,1,ns);
                   (*seisPSV).sectionvy=matrix(1,ntr,1,ns);
                   (*seisPSV).sectionp=matrix(1,ntr,1,ns);
                   break;
           }
   }


}



