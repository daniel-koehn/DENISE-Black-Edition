/* 
   Reading number of source positions.
   
   Daniel Koehn 
   Kiel, 10.03.2019
*/

#include "fd.h"

void count_src(){

	/* declaration of extern variables */
	extern char SOURCE_FILE[STRING_SIZE];
	extern int MYID, NSHOTS;
	FILE *fpsrc;

	if (MYID==0){	

	   fpsrc=fopen(SOURCE_FILE,"r");	
	   if (fpsrc==NULL) err(" Source file could not be opened !");
	   NSHOTS=0;

	   /* read number of source positions */	
	   fscanf(fpsrc,"%d",&NSHOTS);

	   fclose(fpsrc);
						
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&NSHOTS,1,MPI_INT,0,MPI_COMM_WORLD);
	
}
