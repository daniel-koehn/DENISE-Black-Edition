/*------------------------------------------------------------------------
 *   Write 2D snapshot for current timestep  to file                                   
 *   last update 24/05/2002
 *
 *  T. Bohlen
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap_rsg(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx, float **syy, float **u, float **pi){

	/* 
		different data formats of output:
		SNAP_FORMAT=1  :  SU (IEEE)
		SNAP_FORMAT=2  :  ASCII
		SNAP_FORMAT=3  :  BINARY (IEEE)
		
		different types:
		SNAP=1 : values in vx and vy
		SNAP=2 : -(vx+vy) (pressure field)
		SNAP=3 : divergence of vx and vy (energy of compressional waves)
		         and curl of vx and vy (energy of shear waves)
		SNAP=4 : both particle velocities (type=1) and energy (type=3)
		*/


	int i,j;
	float amp, dh24, vyx, vxy, vxx, vyy;
	float vxxs, vxys, vyxs, vyys;
	float **divfield, **curlfield;
	char snapfile_x[STRING_SIZE], snapfile_y[STRING_SIZE], snapfile_div[STRING_SIZE];
	char snapfile_rot[STRING_SIZE], snapfile_p[STRING_SIZE], ext[8], wm[1];
	FILE *fpx1, *fpy1, *fpx2, *fpy2;

	extern float DH, DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int NX, NY,  SNAP_FORMAT, SNAP;
	extern int MYID, POS[3], IDX, IDY;



	switch(SNAP_FORMAT){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}
	sprintf(snapfile_x,"%s%s.x.%i%i",SNAP_FILE,ext,POS[1],POS[2]);
	sprintf(snapfile_y,"%s%s.y.%i%i",SNAP_FILE,ext,POS[1],POS[2]);
	sprintf(snapfile_div,"%s%s.div.%i%i",SNAP_FILE,ext,POS[1],POS[2]);
	sprintf(snapfile_rot,"%s%s.rot.%i%i",SNAP_FILE,ext,POS[1],POS[2]);
	sprintf(snapfile_p,"%s%s.p.%i%i",SNAP_FILE,ext,POS[1],POS[2]);

	fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);
	
	
	if (nsnap==1) 
		sprintf(wm,"w");
	else 
		sprintf(wm,"a");
		

	switch(SNAP){
	case 1 :
		fprintf(fp,"%s\n", snapfile_x);
		fprintf(fp,"%s\n\n", snapfile_y);
		
		fpx1=fopen(snapfile_x,wm);
		fpy1=fopen(snapfile_y,wm);
		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
				amp=vx[j][i];
				writedsk(fpx1,amp,SNAP_FORMAT);
				amp=vy[j][i];
				writedsk(fpy1,amp,SNAP_FORMAT);
			}
		fclose(fpx1);
		fclose(fpy1);
		break;


	case 2 :
		fprintf(fp,"%s\n\n",snapfile_p);
		fpx1=fopen(snapfile_p,wm);

		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
				amp=-sxx[j][i]-syy[j][i];
				writedsk(fpx1,amp,SNAP_FORMAT);
			}
		fclose(fpx1);
		break;

	case 4 :
		fprintf(fp,"%s\n", snapfile_x);
		fprintf(fp,"%s\n", snapfile_y);
		fpx1=fopen(snapfile_x,wm);
		fpy1=fopen(snapfile_y,wm);
		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
				amp=vx[j][i];
				writedsk(fpx1,amp,SNAP_FORMAT);
				amp=vy[j][i];
				writedsk(fpy1,amp,SNAP_FORMAT);
			}
		fclose(fpx1);
		fclose(fpy1);
	case 3 :
		/* output of the curl of the velocity field according to Dougherty and
				                  Stephen (PAGEOPH, 1988) */
		/*if (NY1<=2) error("NY1 must be greater than 2.");*/
		fprintf(fp,"%s\n", snapfile_div);
		fprintf(fp,"%s\n\n", snapfile_rot);
		fpx2=fopen(snapfile_div,wm);
		fpy2=fopen(snapfile_rot,wm);
		curlfield  =  matrix(-1,NY+2,-1,NX+2);

		dh24=1.0/(DH*2.0);
		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
		       		/* in the rotated coordinate system: */	
				vxxs=vx[j][i+1]-vx[j+1][i];
				vxys=vx[j+1][i+1]-vx[j][i];
				vyxs=vy[j][i+1]-vy[j+1][i];
				vyys=vy[j+1][i+1]-vy[j][i];
				/* in the cartesian coordinate system: */
             			vyx=(vyys+vyxs)*dh24;
             			vxy=(vxys-vxxs)*dh24;
				curlfield[j][i]=(vxy-vyx)*sqrt(u[j][i]);
			}

		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
				writedsk(fpy2,curlfield[j][i],SNAP_FORMAT);
			}
		free_matrix(curlfield,-1,NY+2,-1,NX+2);

		/* output of the divergence of the velocity field according to Dougherty and
				                  Stephen (PAGEOPH, 1988) */
		divfield  =  matrix(-1,NY+2,-1,NX+2);
		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
		       		/* in the rotated coordinate system: */	
				vxxs=vx[j][i+1]-vx[j+1][i];
				vxys=vx[j+1][i+1]-vx[j][i];
				vyxs=vy[j][i+1]-vy[j+1][i];
				vyys=vy[j+1][i+1]-vy[j][i];
				/* in the cartesian coordinate system: */
          			vxx=(vxys+vxxs)*dh24;
             			vyy=(vyys-vyxs)*dh24;
				divfield[j][i]=(vxx+vyy)*sqrt(pi[j][i]);
			}


		for (i=1;i<=NX;i+=IDX)
			for (j=1;j<=NY;j+=IDY){
				writedsk(fpx2,divfield[j][i],SNAP_FORMAT);
			}

		free_matrix(divfield,-1,NY+2,-1,NX+2);
		fclose(fpx2);
		fclose(fpy2);
		break;
	}


}


