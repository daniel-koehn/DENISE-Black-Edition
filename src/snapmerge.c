/*------------------------------------------------------------------------
 *   loop over snapshotfiles which have to be merged.                                   
 *   last update 25/05/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */

int main(int argc, char **argv){

int nsnap;
char *fileinp;

/* read parameters from parameter-file (stdin) */
/*read_par(stdin);*/ 


/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open Denise input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
		err(" --- ");
	}
}

/* read input file */
read_par(FP);
	
NXG=NX;
NYG=NY;	
NX = NXG/NPROCX;
NY = NYG/NPROCY;

nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);

FP=stdout;


switch(SNAP){
case 1 : /*particle velocity*/
   merge(nsnap,1);
   merge(nsnap,2);
   break;
case 2 : /*pressure */
   merge(nsnap,6);
   break;
case 4 : /*particle velocity*/
   merge(nsnap,1);
   merge(nsnap,2);
case 3 :
   merge(nsnap,4);
   merge(nsnap,5);
   break;
default :
   warning(" snapmerge: cannot identify content of snapshot !");
   break;

}	
return 0;	

}
