/*------------------------------------------------------------------------
 *  DENISE Black Edition: 2D isotropic elastic time domain FWI Code 
 *
 *
 *  Authors:
 *  -----------  
 * 
 *  D. Koehn    (FWI code + updates)
 *  D. De Nil   (FWI code + updates)
 *  L. Rehor    (viscoelastic modelling, Butterworth-filter)
 *  A. Kurzmann (original step length estimation)
 *  M. Schaefer (source wavelet inversion)
 *  S. Heider   (time-windowing)
 *  T. Bohlen   (original FD forward code) 
 *  L. Zhang    (towed streamer, pressure inversion)
 *  D. Wehner   (gravity modelling and inversion)
 *  V. Kazei    (shot parallelization of the PSV module)
 *  
 *  
 *  In case of questions contact the author:
 *	Dr. Daniel Koehn, Kiel University, Institute of Geoscience,
 *	Otto-Hahn-Platz 1, D-24098 Kiel, Germany, ph: +49 431 880 4566,
 *	mailto:dkoehn@geophysik.uni-kiel.de,
 *	Homepage: http://www.geophysik.uni-kiel.de/~dkoehn
 *
 *
 *  DENISE Black Edition is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, version 2.0 of the License only. 
 *  
 *  DENISE Black Edition is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with DENISE Black Edition (see file LICENSE.md) 
 *
 *  If you show modelling/inversion results in a paper or presentation please 
 *  give a reference to the following papers:
 *
 *  Daniel Koehn, Denise De Nil, Andre Kurzmann, Anna Przebindowska and Thomas Bohlen (2012): 
 *  On the influence of model parametrization in elastic full waveform tomography, 
 *  Geophysical Journal International, 191(1), 325-345.
 *
 *  Daniel Koehn (2011): Time Domain 2D Elastic Full Waveform Tomography, PhD-Thesis, Kiel University
 *  Available at: http://nbn-resolving.de/urn:nbn:de:gbv:8-diss-67866 
 * 
 *  
 *  Thank you for your co-operation, 
 *  Daniel Koehn
 * 
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */
#include "cseife.h"

int main(int argc, char **argv){
char * fileinp;
FILE *fpsrc;

/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

setvbuf(stdout, NULL, _IONBF, 0);
		
/* print program name, version etc to stdout*/
if (MYID == 0) info(stdout);

/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
FILEINP1=argv[2];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open Denise input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
                err(" --- ");
	}
}

/* read input file *.inp */
read_par(FP);
 
/* Init shot parallelization*/
COLOR = MYID / (NPROCX * NPROCY);
MPI_Comm_split(MPI_COMM_WORLD, COLOR, MYID, &SHOT_COMM);
MPI_Comm_rank(SHOT_COMM, &MYID_SHOT);

/* Init subdomain communication*/
MPI_Comm_split(MPI_COMM_WORLD, MYID_SHOT, MYID, &DOMAIN_COMM);

NCOLORS = NP / (NPROCX * NPROCY);

/*printf("MYID: %d \t COLOR: %d \n", MYID, COLOR);
printf("NP: %d \t NCOLORS: %d \n", NP, NCOLORS);
printf("NX: %d \t NY: %d \n", NPROCX, NPROCY);*/

MPI_Barrier(MPI_COMM_WORLD);

count_src();

/*printf("Number of shots %d \n", NSHOTS);*/

/* check if parameters for PHYSICS and MODE are correct */
check_mode_phys();

/* ---------------------------------------------------- */
/* Forward, FWI, RTM and RTMOD modules (2D PSV-problem) */
/* ---------------------------------------------------- */
if(PHYSICS==1){
  physics_PSV();
}

/* ---------------------------------------------------- */
/* Forward, FWI, RTM and RTMOD modules (2D AC-problem)  */
/* ---------------------------------------------------- */
if(PHYSICS==2){
  physics_AC();
}

/* ---------------------------------------------------- */
/* Forward, RTM modules (2D PSV VTI-problem) */
/* ---------------------------------------------------- */
if(PHYSICS==3){
  physics_VTI();
}

/* ---------------------------------------------------- */
/* Forward, RTM modules (2D PSV TTI-problem) */
/* ---------------------------------------------------- */
if(PHYSICS==4){
  physics_TTI();
}

/* ---------------------------------------------------- */
/* Forward, FWI, RTM and RTMOD modules (2D SH-problem)  */
/* ---------------------------------------------------- */
if(PHYSICS==5){
  physics_SH();
}

MPI_Comm_free(&SHOT_COMM);

MPI_Finalize();
return 0;	

}/*main*/
