/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void taper_grad_shot(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int ishot)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2, VTON, SRTON;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1;

	/*extern int GRADT1, GRADT2, GRADT3, GRADT4;*/
	float amp, a, *window, grad_tap, **waveconvtmp;
	char modfile[STRING_SIZE];
	
	extern float SRTRADIUS;
	extern int SRTSHAPE, FILTSIZE;
        float **m, **edgemat, **mm, **msum, minm, maxm, x, y, rad, **taper_coeff_glob;
        float maxrad;

	/*SRTSHAPE=1;
	SRTRADIUS=25.0;
	filtsize=2;*/
	


        /* =================================== */
        /* taper source and receiver positions */
	/* =================================== */
	
        /* Convert from meters to gridpoints -> minimum 5x5 gridpoints */
        srctaper_gridpt = (int)(ceil(2.0*SRTRADIUS/DH));
        if (srctaper_gridpt<5)  srctaper_gridpt = 5;

        m               = matrix(1,srctaper_gridpt,1,srctaper_gridpt);
        edgemat         = matrix(1,4,1,1);
        mm              = matrix(1,NYG,1,NXG);
        msum            = matrix(1,NYG,1,NXG);
        taper_coeff_glob= matrix(1,NYG,1,NXG);
	waveconvtmp     = matrix(0,NY+1,0,NX+1);

        for (iy=1;iy<=NYG;iy++)
                for (ix=1;ix<=NXG;ix++)  msum[iy][ix] = 1.0;

        MPI_Barrier(MPI_COMM_WORLD);

        /*****************************/
        /* Taper at source positions */
        /*****************************/
                
	a = 1.0;
        maxrad = sqrt(2.0*SRTRADIUS*SRTRADIUS);
        for (j=1;j<=srctaper_gridpt;j++) {
                for (i=1;i<=srctaper_gridpt;i++) {
                        x = ((float)i-((float)srctaper_gridpt)/2.0-0.5)*DH;
                        y = ((float)j-((float)srctaper_gridpt)/2.0-0.5)*DH;
                        rad = sqrt(x*x+y*y);

                        switch (SRTSHAPE) {
                        case 1:
                                m[j][i] = erf(a*rad/maxrad);
                                break;
                        case 2:
                                if (rad>0)      m[j][i] = log(rad);
                                else            m[j][i] = 0.0;
                                break;
                        }
                }
        }

        /* generate local taper matrix */
        minm = minimum_m(m,srctaper_gridpt,srctaper_gridpt);
        for (j=1;j<=srctaper_gridpt;j++)
                for (i=1;i<=srctaper_gridpt;i++)  m[j][i] -= minm;

        /* normalize taper matrix to max of values at the centre of all 4 taper area edges,     */
        /* not the global maximum, which is located at the corners                              */
        edgemat[1][1] = m[1][srctaper_gridpt/2];
        edgemat[2][1] = m[srctaper_gridpt/2][1];
        edgemat[3][1] = m[srctaper_gridpt/2][srctaper_gridpt];
        edgemat[4][1] = m[srctaper_gridpt][srctaper_gridpt/2];
        maxm = maximum_m(edgemat,1,4);
        for (j=1;j<=srctaper_gridpt;j++)
                for (i=1;i<=srctaper_gridpt;i++) {
                        m[j][i] /= maxm;
                        if (m[j][i]>1.0)  m[j][i] = 1.0;
                }
        /* get central position within the taper */
        ijc = (int)(ceil((float)srctaper_gridpt/2));

        /*********************/
        /* loop over sources */
        /*for (n=1;n<=nshots;n++) {*/
	n=ishot;
        for (iy=1;iy<=NYG;iy++)
                for (ix=1;ix<=NXG;ix++)  mm[iy][ix] = 1.0;

        i = iround(srcpos[1][n]/DH);
        j = iround(srcpos[2][n]/DH);
        for (iy=1;iy<=srctaper_gridpt;iy++) {
                for (ix=1;ix<=srctaper_gridpt;ix++) {
                        xx = i + ix - ijc;
                        yy = j + iy - ijc;
                        if ((xx<1) || (xx>NXG) || (yy<1) || (yy>NYG))  continue;
                        mm[yy][xx] = m[iy][ix];
                }
        }

/*                      for (iy=1;iy<=NYG;iy++)
                                for (ix=1;ix<=NXG;ix++)  msum[iy][ix] += mm[iy][ix];
*/
        for (iy=1;iy<=NYG;iy++)
               for (ix=1;ix<=NXG;ix++)
                        if (msum[iy][ix] > mm[iy][ix])
                                msum[iy][ix] = mm[iy][ix];

               /* }*/

        /***********************/
        /* loop over receivers */
        /*for (n=1;n<=ntr;n++) {
                for (iy=1;iy<=NYG;iy++)
                        for (ix=1;ix<=NXG;ix++)  mm[iy][ix] = 1.0;
                 i = recpos[1][n];
                 j = recpos[2][n];
                 for (iy=1;iy<=srctaper_gridpt;iy++) {
                        for (ix=1;ix<=srctaper_gridpt;ix++) {
                                xx = i + ix - ijc;
                                yy = j + iy - ijc;
                                if ((xx<1) || (xx>NXG) || (yy<1) || (yy>NYG))  continue;
                                mm[yy][xx] = m[iy][ix];
                        }
                }*/

/*                      for (iy=1;iy<=NYG;iy++)    Die kommenden zwei Zeilen wurden von Daniel auskommentiert.
                                for (ix=1;ix<=NXG;ix++)  msum[iy][ix] += mm[iy][ix];
*/
                       /* for (iy=1;iy<=NYG;iy++)
                                for (ix=1;ix<=NXG;ix++)
                                        if (msum[iy][ix] > mm[iy][ix])
                                                msum[iy][ix] = mm[iy][ix];

                }*/

 
        minm = minimum_m(msum,NXG,NYG);
        for (iy=1;iy<=NYG;iy++)
                for (ix=1;ix<=NXG;ix++)  msum[iy][ix] -= minm;

        maxm = maximum_m(msum,NXG,NYG);
        for (iy=1;iy<=NYG;iy++)
                for (ix=1;ix<=NXG;ix++) {
                        taper_coeff_glob[iy][ix] = msum[iy][ix]/maxm;
                         if ((POS[1]==((ix-1)/NX)) && (POS[2]==((iy-1)/NY))){
                                ii = ix-POS[1]*NX;
                                jj = iy-POS[2]*NY;
                                /*Diese Zeile wurde von Daniel auskommentiert. taper_coeff[jj][ii] = taper_coeff_glob[iy][ix] * ((float)(iy*DH));*/
				/*taper_coeff[jj][ii] = ((float)(iy*DH)) * ((float)(iy*DH)) * ((float)(iy*DH));*/
				taper_coeff[jj][ii]=msum[iy][ix]/maxm;
                        }
                }
	
			
	/* apply taper on local gradient */
        for (j=1;j<=NY;j++){
           for (i=1;i<=NX;i++){
           waveconv[j][i]*=taper_coeff[j][i];
           waveconvtmp[j][i] = waveconv[j][i];
           }
        }	
	
	
	 /* apply filter at shot and receiver points */
	/*for (n=1;n<=nshots;n++)
	{*/
	n=ishot;
		i1 = iround(srcpos[1][n]/DH);
		j1 = iround(srcpos[2][n]/DH);
		fprintf(FP,"\n Shot %d (printed by PE %d):\n",n,MYID);
		fprintf(FP,"\n i1: %d (printed by PE %d):\n",i1,MYID);
		fprintf(FP,"\n j1: %d (printed by PE %d):\n",j1,MYID);
		
		for (i=i1-FILTSIZE;i<=i1+FILTSIZE;i++){
		   for (j=j1-FILTSIZE;j<=j1+FILTSIZE;j++){
		       if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
			     ii = i-POS[1]*NX;
			     jj = j-POS[2]*NY;
			     /*printf("\n ii: %d (printed by PE %d):\n",ii,MYID);
			     printf("\n jj: %d (printed by PE %d):\n",jj,MYID);*/
			     /*waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
				       		+ 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
						+ 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
			      waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;*/
			      if (jj>0){
			         waveconvtmp[jj][ii] = 0.0;
			         taper_coeff[jj][ii] = 0.0;
			      }
		        }
		    }
		}
		
		
	/*}*/

	
	/*for (n=1;n<=ntr;n++)
	{
		i1 = recpos[1][n];
		j1 = recpos[2][n];
		
            for (i=i1-filtsize;i<=i1+filtsize;i++){
		   for (j=j1-filtsize;j<=j1+filtsize;j++){
		
		       if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
			    ii = i-POS[1]*NX;
			    jj = j-POS[2]*NY;   */
			    /* Die kommenden 4 Zeilen wurden von Daniel auskommentiert. waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
						+ 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
						+ 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
			      waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;*/
			
		/*	waveconvtmp[jj][ii] = 0.0;
			
		       }
		   }
	     }	       
	}*/	
	
		
        /* apply taper on local gradient */
        for (j=1;j<=NY;j++){
           for (i=1;i<=NX;i++){
           waveconv[j][i] = waveconvtmp[j][i];
           }
        }
	
	
                free_matrix(m,1,srctaper_gridpt,1,srctaper_gridpt);
                free_matrix(edgemat,1,4,1,1);
                free_matrix(mm,1,NYG,1,NXG);
                free_matrix(msum,1,NYG,1,NXG);
                free_matrix(taper_coeff_glob,1,NYG,1,NXG);
		free_matrix(waveconvtmp,0,NX+1,0,NY+1);
        
	
	

 
		
	MPI_Barrier(MPI_COMM_WORLD);
        sprintf(modfile,"taper_coeff_%i.bin",ishot);
        writemod(modfile,taper_coeff,3); 
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(modfile,3); 
	



	
	
}



