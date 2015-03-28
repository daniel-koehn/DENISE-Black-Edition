/*
 *   Create dissipative boundarie around the model grid
 *   The dissipative coefficients are stored in the 2-D array
 *   absorb_coeff. The interior of the model is weighted by the
 *   coefficient 1. In the absorbing frame the coefficients 
 *   are less than one. Coefficients are computed using 
 *   exponential damping (see Cerjan et al., 1985, 
 *   Geophysics, 50, 705-708)
 *
 *   last update K. Witte 22.12.01 :-)
 */

#include "fd.h"

void absorb(float ** absorb_coeff)
{

	/* extern variables */

	extern float DH, FW, DAMPING;
	extern int FREE_SURF, NX, NY, BOUNDARY;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern FILE *FP;
	
	/* local variables */
	int i, j, ifw, ii, jj, xb, yb, xe, ye;
	float amp, a, *coeff;
	/*char modfile[STRING_SIZE];*/
	
	if (MYID==0)
	{
		fprintf(FP,"\n **Message from absorb (printed by PE %d):\n",MYID);
		fprintf(FP," Coefficcients for absorbing frame are now calculated.\n");
		fprintf(FP," Width of dissipative frame (meter)= %f\n",FW);
		fprintf(FP," Percentage of exponential damping = %5.2f\n",DAMPING);
	}

	amp=1.0-DAMPING/100.0;   /* amplitude at the edge of the numerical grid */
	ifw=iround(FW/DH);  /* frame width in gridpoints */
	coeff=vector(1,ifw);
	a=sqrt(-log(amp)/((ifw-1)*(ifw-1)));
	
	for (i=1;i<=ifw;i++)
		coeff[i]=exp(-(a*a*(ifw-i)*(ifw-i)));
	
	if (MYID==0)
	{
		fprintf(FP," Table of coefficients \n # \t coeff \n");
		/*printf(" ifw=%d \t a=%f amp=%f \n", ifw,a,amp); */
		for (i=1;i<=ifw;i++)
			fprintf(FP," %d \t %5.3f \n", i, coeff[i]);
	}	
	

	/* initialize array of coefficients with one */
	for (j=1;j<=NY;j++)
	for (i=1;i<=NX;i++) absorb_coeff[j][i]=1.0;


	/* compute coefficients for left and right grid boundaries (x-direction) */
	if ((!BOUNDARY) && (POS[1]==0))
	{
		yb=1; ye=NY; 
		for (i=1;i<=ifw;i++)
		{
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][i]=coeff[i];
		}
	}
			
	if ((!BOUNDARY) && (POS[1]==NPROCX-1))
	{
		yb=1; ye=NY;
		for (i=1;i<=ifw;i++){
			ii=NX-i+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][ii]=coeff[i];
		}
	}
	

	/* compute coefficients for top and bottom grid boundaries (y-direction) */

	if ((POS[2]==0) && (!(FREE_SURF)))
	{
		xb=1; xe=NX;
		for (j=1;j<=ifw;j++)
		{
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[j][i]=coeff[j];
		}
	}

	if (POS[2]==NPROCY-1)
	{
		xb=1; xe=NX;
		for (j=1;j<=ifw;j++)
		{
			jj=NY-j+1;
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[jj][i]=coeff[j];
		}
	}


/*	sprintf(modfile,"absorb_coeff.bin");

	writemod(modfile,absorb_coeff,3); 

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); 
*/

	free_vector(coeff,1,ifw);
}



