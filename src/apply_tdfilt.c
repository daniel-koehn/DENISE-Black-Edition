/*
 * Apply time-domain filter
 *
 * Daniel Koehn
 * Kiel, 23/04/2016
 */

#include "fd.h"

void apply_tdfilt(float **section, int ntr, int ns, int order, float fc2, float fc1){

     /* global variables */
     extern int TIME_FILT;

     /* local variables */

     timedomain_filt(section,fc2,order,ntr,ns,1);

     if(TIME_FILT==2){ /* band-pass */
       timedomain_filt(section,fc1,order,ntr,ns,2);
     } 
       	
}



