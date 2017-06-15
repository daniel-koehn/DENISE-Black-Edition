/*
 * Calculate memory requirements for acoustic problem 
 *
 * Daniel Koehn
 * Kiel, 12/06/2017
 */

#include "fd.h"

void mem_fwiAC(int nseismograms, int ntr, int ns, int fdo3, int nd, float buffsize, int ntr_glob){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER, MYID, IDXI, IDYI, NTDTINV, QUELLTYPB;

	/* local variables */
        float fac1, fac2, memdyn, memmodel, memseismograms, memfwt, memfwt1, memfwtdata;
        float memseismograms1, membuffer, memtotal;

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

        /* storage of forward modelled wavefield */
	memfwt=5.0*(NX/IDXI)*(NY/IDYI)*NTDTINV*fac2;
        
        /* gradients, updated model parameters   */
	memfwt1=23.0*NX*NY*fac2;					  

        /* residual seismograms*/
	memfwtdata=3.0*ntr*ns*fac2;
        if(QUELLTYPB==1){memfwtdata=3.0*ntr*ns*fac2;}

        /* full seismograms */
	memseismograms1=(nseismograms+1)*ntr_glob*ns*fac2;

	membuffer=2.0*fdo3*(NY+NX)*fac2;
	buffsize=2.0*2.0*fdo3*(NX +NY)*sizeof(MPI_FLOAT);
	memtotal=memdyn+memmodel+memseismograms+memfwt+memfwt1+memfwtdata+membuffer+(buffsize*pow(2.0,-20.0));

	if (MYID==0){
	   printf("\n **Message from main (printed by PE %d):\n",MYID);
	   printf(" Size of local grids: NX=%d \t NY=%d\n",NX,NY);
	   printf(" Each process is now trying to allocate memory for:\n");
	   printf(" Dynamic variables: \t\t %6.2f MB\n", memdyn);
	   printf(" Static variables: \t\t %6.2f MB\n", memmodel);
	   printf(" Seismograms: \t\t\t %6.2f MB\n", memseismograms);
	   printf(" Buffer arrays for grid exchange: %6.2f MB\n", membuffer);
	   printf(" Storage of forward modelled wavefields: %6.2f MB\n", memfwt);
	   printf(" Gradients, updated material parameters: %6.2f MB\n", memfwt1);
	   printf(" Residual seismograms: \t %6.2f MB\n", memfwtdata);
	   printf(" Full seismograms: \t %6.2f MB\n", memseismograms1);
	   printf(" Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
	   printf(" ------------------------------------------------ \n");
	   printf(" Total memory required: \t %6.2f MB.\n\n", memtotal);
	}
        	
}



