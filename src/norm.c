/*
 * Normalize Gradient
 *
 * Daniel Koehn
 */

#include "fd.h"

float norm(float ** waveconv, int iter, int sws)
{

	/* extern variables */

    extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2, SPAT_FILT_SIZE, SPAT_FILT_1;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1, filtsize;

	float **waveconvtmp, **waveconvtmps, grad, normgauss;
	int winx, winy;
	float min2, max2, min1, max1, norm;
	float tmp_max, tmp_min;
	FILE *fp_grad;
	
	min2 = waveconv[1][1];
	max2 = waveconv[1][1];
	
	/* estimate min and max of the global gradient */
	 for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){
	            
	       if(waveconv[j][i]<min2){min2=waveconv[j][i];}
	       if(waveconv[j][i]>max2){max2=waveconv[j][i];}
	            
	    }
	 }
	
   /* Compute maximum value on all CPUs */
   tmp_max = 0.0;
   MPI_Allreduce(&max2,&tmp_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
   max2 = tmp_max;

   /* Compute minimum value on all CPUs */
   tmp_min = 0.0;
   MPI_Allreduce(&min2,&tmp_min,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
   min2 = tmp_min;

	/* Normalize gradient */
	norm = max2;
	
	/*if(sqrt(min2*min2)>sqrt(max2*max2)){norm=sqrt(min2*min2);}*/
	
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){	    	    
	       waveconv[j][i] = waveconv[j][i]/norm;	    
	    }
	}                
	
    printf("MYID=%d \t norm_fac=%e \n",MYID,norm);
    return norm;
}



