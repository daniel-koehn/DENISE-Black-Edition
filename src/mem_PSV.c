/*
 * Calculate memory requirements for PSV problem 
 *
 * Daniel Koehn
 * Kiel, 23/04/2016
 */

#include "fd.h"

void mem_PSV(int nseismograms,int ntr, int ns, int fdo3, int nd, float buffsize){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER, MYID, IDXI, IDYI, NTDTINV;

	/* local variables */
        float fac1, fac2, memdyn, memmodel, memseismograms, memfwt, memfwt1, memfwtdata;
        float membuffer, memtotal;

	/*allocate memory for dynamic, static and buffer arrays */
	fac1=(NX+FDORDER)*(NY+FDORDER);
	fac2=sizeof(float)*pow(2.0,-20.0);

	if (L){
		memdyn=(5.0+3.0*(float)L)*fac1*fac2;
		memmodel=(12.0+3.0*(float)L)*fac1*fac2;
	
	} else {
		memdyn=5.0*fac1*fac2;
		memmodel=6.0*fac1*fac2;
	}

	memseismograms=nseismograms*ntr*ns*fac2;

	membuffer=2.0*fdo3*(NY+NX)*fac2;
	buffsize=2.0*2.0*fdo3*(NX +NY)*sizeof(MPI_FLOAT);
	memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0));

	if (MYID==0){
	   printf("\n **Message from main (printed by PE %d):\n",MYID);
	   printf(" Size of local grids: NX=%d \t NY=%d\n",NX,NY);
	   printf(" Each process is now trying to allocate memory for:\n");
	   printf(" Dynamic variables: \t\t %6.2f MB\n", memdyn);
	   printf(" Static variables: \t\t %6.2f MB\n", memmodel);
	   printf(" Seismograms: \t\t\t %6.2f MB\n", memseismograms);
	   printf(" Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
	   printf(" Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
	   printf(" ------------------------------------------------ \n");
	   printf(" Total memory required: \t %6.2f MB.\n\n", memtotal);
	}
        	
}



