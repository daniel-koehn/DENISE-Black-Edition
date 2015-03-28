/*------------------------------------------------------------------------
 *   Write one single amplitude on disk                                   
 *   last update 07/10/98
 *
 *  Copyright (c) 1997, 1998 T. Bohlen
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


void writedsk(FILE *fp_out, float amp, int format){

#ifdef _CRAY1     
	float ampc;
	int type=2, num=1, bitoff=0, stride=1, ierr;
#endif


	switch(format){
	case 1 : /* SU*/

		err(" Sorry, SU-format for snapshots not implemented yet. \n");
		break;
	case 2 :  /*ASCII*/
		fprintf(fp_out,"%e\n", amp);
		break;
	case 3 :   /* BINARY */
#ifdef _CRAY1     
		ierr=CRAY2IEG(&type, &num, &ampc, &bitoff, 
		    &amp,&stride);
		fwrite(&ampc,4,1,fp_out);
#else 
		fwrite(&amp, sizeof(float), 1, fp_out);
#endif
		break;

	default :
		printf(" Don't know the format for the snapshot-data !\n");
		err(" No output was written. ");
	}
}
