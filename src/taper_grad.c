/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void taper_grad(float **waveconv, float **taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int sws)
{

        /* extern variables */

        extern float DH;
        extern float SRTRADIUS, EXP_TAPER_GRAD_HOR;
        extern int SRTSHAPE, FILTSIZE;
        extern int FREE_SURF, NX, NY, NXG, NYG;
        extern int NPROCX, NPROCY, MYID, POS[3];
        extern int GRADT1, GRADT2, GRADT3, GRADT4;
        extern char TFILE[STRING_SIZE];
        extern FILE *FP;

        /* local variables */
        int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength, taperlength2, VTON, SRTON;
        int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1;


        float amp, a, *window, grad_tap, **waveconvtmp;
        char modfile[STRING_SIZE];

        float **m, **edgemat, **mm, **msum, minm, maxm, x, y, rad, **taper_coeff_glob;
        float maxrad;

        FILE *fp_taper;

        /*SRTSHAPE=2;
	SRTRADIUS=25.0;
	filtsize=2;*/

        /* =============== */
        /* Vertical taper  */
        /* =============== */

        if (sws == 1)
        {

                /* define taper geometry */
                taperlength = GRADT2 - GRADT1;
                taperlength2 = GRADT4 - GRADT3;
                ifw = GRADT4 - GRADT1;

                /*printf("%d \t %d \t %d \t %d \n",GRADT1, GRADT2, GRADT3, GRADT4);
	printf("%d \t %d \t %d \n",taperlength, taperlength2,ifw);*/

                if (MYID == 0)
                {
                        fprintf(FP, "\n **Message from taper_grid (printed by PE %d):\n", MYID);
                        fprintf(FP, " Coefficients for gradient taper are now calculated.\n");
                }

                waveconvtmp = matrix(0, NY + 1, 0, NX + 1);
                window = vector(1, ifw);

                /* Gaussian window */
                a = 3; /* damping coefficient */
                for (i = 1; i <= ifw; i++)
                {
                        window[i] = 1.0;

                        if (i <= taperlength)
                        {
                                window[i] = exp(-(1.0 / 2.0) * (a * (i - taperlength) / (taperlength / 2.0)) * (a * (i - taperlength) / (taperlength / 2.0)));
                        }

                        if (i >= (ifw - taperlength2))
                        {
                                window[i] = exp(-(1.0 / 2.0) * (a * (i - (ifw - taperlength2)) / (taperlength2 / 2.0)) * (a * (i - (ifw - taperlength2)) / (taperlength2 / 2.0)));
                        }
                }

                /* loop over global grid */
                for (j = 1; j <= NYG; j++)
                {

                        h = 1;
                        for (i = 1; i <= NXG; i++)
                        {

                                grad_tap = 0.0;

                                if ((i > GRADT1) && (i < GRADT4))
                                {
                                        grad_tap = window[h];
                                        h++;
                                }

                                if ((POS[1] == ((i - 1) / NX)) &&
                                    (POS[2] == ((j - 1) / NY)))
                                {
                                        ii = i - POS[1] * NX;
                                        jj = j - POS[2] * NY;

                                        taper_coeff[jj][ii] = grad_tap;
                                }
                        }
                }

                /* apply taper on local gradient */
                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] *= taper_coeff[j][i];
                                waveconvtmp[j][i] = waveconv[j][i];
                        }
                }

                /* apply filter at shot and receiver points */
                for (n = 1; n <= nshots; n++)
                {
                        i = iround(srcpos[1][n] / DH);
                        j = iround(srcpos[2][n] / DH);
                        if ((POS[1] == ((i - 1) / NX)) && (POS[2] == ((j - 1) / NY)))
                        {
                                ii = i - POS[1] * NX;
                                jj = j - POS[2] * NY;
                                /*waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
						+ 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
						+ 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
			waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;*/

                                waveconvtmp[jj][ii] = 0.0;
                        }
                }

                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] = waveconvtmp[j][i];
                        }
                }

                for (n = 1; n <= ntr; n++)
                {
                        i = recpos[1][n];
                        j = recpos[2][n];
                        if ((POS[1] == ((i - 1) / NX)) && (POS[2] == ((j - 1) / NY)))
                        {
                                ii = i - POS[1] * NX;
                                jj = j - POS[2] * NY;
                                /*waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
						+ 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
						+ 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
			waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;*/

                                waveconvtmp[jj][ii] = 0.0;
                        }
                }

                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] = waveconvtmp[j][i];
                        }
                }

                sprintf(modfile, "%s_coeff_vert.bin", TFILE);

                writemod(modfile, taper_coeff, 3);

                MPI_Barrier(MPI_COMM_WORLD);

                if (MYID == 0)
                        mergemod(modfile, 3);

                free_vector(window, 1, ifw);
                free_matrix(waveconvtmp, 0, NX + 1, 0, NY + 1);
        }

        /* ======================== */
        /* Horizontal Taper         */
        /* ======================== */

        if (sws == 2)
        {
                /* define taper geometry */
                taperlength = GRADT2 - GRADT1;
                taperlength2 = GRADT4 - GRADT3;
                ifw = GRADT2 - GRADT1 + 1;

                if (MYID == 0)
                {
                        fprintf(FP, "\n **Message from taper_grid (printed by PE %d):\n", MYID);
                        fprintf(FP, " Coefficients for gradient taper are now calculated.\n");
                        printf("%d \t %d \t %d \t %d \n", GRADT1, GRADT2, GRADT3, GRADT4);
                        printf("%d \t %d \t %d \n", taperlength, taperlength2, ifw);
                }

                waveconvtmp = matrix(0, NY + 1, 0, NX + 1);
                window = vector(1, ifw);

                /* Gaussian window */
                a = 3; /* damping coefficient */
                for (i = 1; i <= ifw; i++)
                {
                        window[i] = 1.0;

                        if (i <= taperlength)
                        {
                                window[i] = exp(-(1.0 / 2.0) * (a * (i - taperlength) / (taperlength / 2.0)) * (a * (i - taperlength) / (taperlength / 2.0)));
                        }
                }

                /* loop over global grid */
                h = 1;
                for (j = 1; j <= NYG; j++)
                {
                        for (i = 1; i <= NXG; i++)
                        {

                                grad_tap = 0.0;

                                if ((j >= GRADT1) && (j <= GRADT2))
                                {
                                        grad_tap = window[h];
                                }

                                if (j > GRADT2)
                                {

                                        /*grad_tap=((float)(j*DH))*((float)(j*DH))*((float)(j*DH));*/
                                        grad_tap = ((float)(j * DH)) * ((float)(j * DH));
                                        /*grad_tap=(float)(j*DH);*/

                                        /*grad_tap=((float)((j*j)/(GRADT2*GRADT2)));*/
                                        /*grad_tap=1.0;*/

                                        grad_tap = pow((float)(j * DH), EXP_TAPER_GRAD_HOR);
                                }

                                /*grad_tap=((float)(j*DH));*/

                                if ((POS[1] == ((i - 1) / NX)) &&
                                    (POS[2] == ((j - 1) / NY)))
                                {
                                        ii = i - POS[1] * NX;
                                        jj = j - POS[2] * NY;

                                        taper_coeff[jj][ii] = grad_tap;
                                }
                        }

                        if (j >= GRADT1)
                        {
                                h++;
                        }
                }

                /* apply taper on local gradient */
                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] *= taper_coeff[j][i];
                                waveconvtmp[j][i] = waveconv[j][i];
                        }
                }

                /* apply filter at shot and receiver points */
                /*for (n=1;n<=nshots;n++)
        {
                i = iround(srcpos[1][n]/DH);
                j = iround(srcpos[2][n]/DH);
                if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
                        ii = i-POS[1]*NX;
                        jj = j-POS[2]*NY;
                        waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
                                                + 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
                                                + 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
                        waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;
                }
        }

        for (j=1;j<=NY;j++){
                for (i=1;i<=NX;i++){
                        waveconv[j][i] = waveconvtmp[j][i];
                }
        }

        for (n=1;n<=ntr;n++)
        {
                i = recpos[1][n];
                j = recpos[2][n];
                if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
                        ii = i-POS[1]*NX;
                        jj = j-POS[2]*NY;
                        waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
                                                + 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
                                                + 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
                        waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;
                }
        }

        for (j=1;j<=NY;j++){
                for (i=1;i<=NX;i++){
                        waveconv[j][i] = waveconvtmp[j][i];
                }
        }*/

                sprintf(modfile, "%s_coeff_hor.bin", TFILE);
                writemod(modfile, taper_coeff, 3);
                MPI_Barrier(MPI_COMM_WORLD);
                if (MYID == 0)
                        mergemod(modfile, 3);

                free_vector(window, 1, ifw);
                free_matrix(waveconvtmp, 0, NX + 1, 0, NY + 1);
        } /* end of sws==2 */

        /* =================================== */
        /* taper source and receiver positions */
        /* =================================== */

        if (sws == 3)
        {
                /* Convert from meters to gridpoints -> minimum 5x5 gridpoints */
                srctaper_gridpt = (int)(ceil(2.0 * SRTRADIUS / DH));
                if (srctaper_gridpt < 5)
                        srctaper_gridpt = 5;

                m = matrix(1, srctaper_gridpt, 1, srctaper_gridpt);
                edgemat = matrix(1, 4, 1, 1);
                mm = matrix(1, NYG, 1, NXG);
                msum = matrix(1, NYG, 1, NXG);
                taper_coeff_glob = matrix(1, NYG, 1, NXG);
                waveconvtmp = matrix(0, NY + 1, 0, NX + 1);

                for (iy = 1; iy <= NYG; iy++)
                        for (ix = 1; ix <= NXG; ix++)
                                msum[iy][ix] = 1.0;

                MPI_Barrier(MPI_COMM_WORLD);

                /*****************************/
                /* Taper at source positions */
                /*****************************/

                a = 1.0;
                maxrad = sqrt(2.0 * SRTRADIUS * SRTRADIUS);
                for (j = 1; j <= srctaper_gridpt; j++)
                {
                        for (i = 1; i <= srctaper_gridpt; i++)
                        {
                                x = ((float)i - ((float)srctaper_gridpt) / 2.0 - 0.5) * DH;
                                y = ((float)j - ((float)srctaper_gridpt) / 2.0 - 0.5) * DH;
                                rad = sqrt(x * x + y * y);

                                switch (SRTSHAPE)
                                {
                                case 1:
                                        m[j][i] = erf(a * rad / maxrad);
                                        break;
                                case 2:
                                        if (rad > 0)
                                                m[j][i] = log(rad);
                                        else
                                                m[j][i] = 0.0;
                                        break;
                                }
                        }
                }

                /* generate local taper matrix */
                minm = minimum_m(m, srctaper_gridpt, srctaper_gridpt);
                for (j = 1; j <= srctaper_gridpt; j++)
                        for (i = 1; i <= srctaper_gridpt; i++)
                                m[j][i] -= minm;

                /* normalize taper matrix to max of values at the centre of all 4 taper area edges,     */
                /* not the global maximum, which is located at the corners                              */
                edgemat[1][1] = m[1][srctaper_gridpt / 2];
                edgemat[2][1] = m[srctaper_gridpt / 2][1];
                edgemat[3][1] = m[srctaper_gridpt / 2][srctaper_gridpt];
                edgemat[4][1] = m[srctaper_gridpt][srctaper_gridpt / 2];
                maxm = maximum_m(edgemat, 1, 4);
                for (j = 1; j <= srctaper_gridpt; j++)
                        for (i = 1; i <= srctaper_gridpt; i++)
                        {
                                m[j][i] /= maxm;
                                if (m[j][i] > 1.0)
                                        m[j][i] = 1.0;
                        }

                /* get central position within the taper */
                ijc = (int)(ceil((float)srctaper_gridpt / 2));

                /*********************/
                /* loop over sources */
                for (n = 1; n <= nshots; n++)
                {
                        for (iy = 1; iy <= NYG; iy++)
                                for (ix = 1; ix <= NXG; ix++)
                                        mm[iy][ix] = 1.0;

                        i = iround(srcpos[1][n] / DH);
                        j = iround(srcpos[2][n] / DH);
                        for (iy = 1; iy <= srctaper_gridpt; iy++)
                        {
                                for (ix = 1; ix <= srctaper_gridpt; ix++)
                                {
                                        xx = i + ix - ijc;
                                        yy = j + iy - ijc;
                                        if ((xx < 1) || (xx > NXG) || (yy < 1) || (yy > NYG))
                                                continue;
                                        mm[yy][xx] = m[iy][ix];
                                }
                        }

                        /*                      for (iy=1;iy<=NYG;iy++)
                                for (ix=1;ix<=NXG;ix++)  msum[iy][ix] += mm[iy][ix];
*/
                        for (iy = 1; iy <= NYG; iy++)
                                for (ix = 1; ix <= NXG; ix++)
                                        if (msum[iy][ix] > mm[iy][ix])
                                                msum[iy][ix] = mm[iy][ix];
                }

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

                minm = minimum_m(msum, NXG, NYG);
                for (iy = 1; iy <= NYG; iy++)
                        for (ix = 1; ix <= NXG; ix++)
                                msum[iy][ix] -= minm;

                maxm = maximum_m(msum, NXG, NYG);
                for (iy = 1; iy <= NYG; iy++)
                        for (ix = 1; ix <= NXG; ix++)
                        {
                                taper_coeff_glob[iy][ix] = msum[iy][ix] / maxm;

                                if ((POS[1] == ((ix - 1) / NX)) && (POS[2] == ((iy - 1) / NY)))
                                {
                                        ii = ix - POS[1] * NX;
                                        jj = iy - POS[2] * NY;
                                        /*Diese Zeile wurde von Daniel auskommentiert: taper_coeff[jj][ii] = taper_coeff_glob[iy][ix] * ((float)(iy*DH));*/
                                        /*taper_coeff[jj][ii] = ((float)(iy*DH)) * ((float)(iy*DH)) * ((float)(iy*DH));*/
                                        taper_coeff[jj][ii] = msum[iy][ix] / maxm;
                                }
                        }

                /* apply taper on local gradient */
                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] *= taper_coeff[j][i];
                                waveconvtmp[j][i] = waveconv[j][i];
                        }
                }

                /* apply filter at shot and receiver points */
                for (n = 1; n <= nshots; n++)
                {
                        i1 = iround(srcpos[1][n] / DH);
                        j1 = iround(srcpos[2][n] / DH);

                        for (i = i1 - FILTSIZE; i <= i1 + FILTSIZE; i++)
                        {
                                for (j = j1 - FILTSIZE; j <= j1 + FILTSIZE; j++)
                                {
                                        if ((POS[1] == ((i - 1) / NX)) && (POS[2] == ((j - 1) / NY)))
                                        {
                                                ii = i - POS[1] * NX;
                                                jj = j - POS[2] * NY;
                                                /*waveconvtmp[jj][ii] = 1*waveconv[jj][ii]
				       		+ 8*(waveconv[jj-1][ii] + waveconv[jj+1][ii] + waveconv[jj][ii-1] + waveconv[jj][ii+1])
						+ 4*(waveconv[jj-1][ii-1] + waveconv[jj-1][ii+1] + waveconv[jj+1][ii+1] + waveconv[jj+1][ii-1]);
			      waveconvtmp[jj][ii] = waveconvtmp[jj][ii]/49;*/
                                                if (jj > 0)
                                                {
                                                        waveconvtmp[jj][ii] = 0.0;
                                                        taper_coeff[jj][ii] = 0.0;
                                                }
                                        }
                                }
                        }
                }

                /*for (n=1;n<=ntr;n++)
	{
		i1 = recpos[1][n];
		j1 = recpos[2][n];
		
            for (i=i1-FILTSIZE;i<=i1+FILTSIZE;i++){
		   for (j=j1-FILTSIZE;j<=j1+FILTSIZE;j++){
		
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
                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] = waveconvtmp[j][i];
                        }
                }

                free_matrix(m, 1, srctaper_gridpt, 1, srctaper_gridpt);
                free_matrix(edgemat, 1, 4, 1, 1);
                free_matrix(mm, 1, NYG, 1, NXG);
                free_matrix(msum, 1, NYG, 1, NXG);
                free_matrix(taper_coeff_glob, 1, NYG, 1, NXG);
                free_matrix(waveconvtmp, 0, NX + 1, 0, NY + 1);

                MPI_Barrier(MPI_COMM_WORLD);
                sprintf(modfile, "%s_coeff_sources.bin", TFILE);
                writemod(modfile, taper_coeff, 3);
                MPI_Barrier(MPI_COMM_WORLD);
                if (MYID == 0)
                        mergemod(modfile, 3);
        }

        /* ======================== */
        /* Read Taper from file     */
        /* ======================== */

        if ((sws >= 4) && (sws <= 6))
        {

                if (MYID == 0)
                {
                        fprintf(FP, "\n **Message from taper_grid (printed by PE %d):\n", MYID);
                        fprintf(FP, " Coefficients for gradient taper are now calculated.\n");
                }

                if (sws == 4)
                {
			sprintf(modfile, "%s.bin", TFILE);
                        fp_taper = fopen(modfile, "r");
                }
                if (sws == 5)
                {
			sprintf(modfile, "%s_u.bin", TFILE);
                        fp_taper = fopen(modfile, "r");
                }
                if (sws == 6)
                {
			sprintf(modfile, "%s_rho.bin", TFILE);
                        fp_taper = fopen(modfile, "r");
                }

                /* loop over global grid */
                for (i = 1; i <= NXG; i++)
                {
                        for (j = 1; j <= NYG; j++)
                        {

                                fread(&grad_tap, sizeof(float), 1, fp_taper);

                                if ((POS[1] == ((i - 1) / NX)) && (POS[2] == ((j - 1) / NY)))
                                {
                                        ii = i - POS[1] * NX;
                                        jj = j - POS[2] * NY;

                                        taper_coeff[jj][ii] = grad_tap;
                                }
                        }
                }

                for (j = 1; j <= NY; j++)
                {
                        for (i = 1; i <= NX; i++)
                        {
                                waveconv[j][i] *= taper_coeff[j][i];
                        }
                }

                fclose(fp_taper);

                sprintf(modfile, "%s_coeff_file.bin", TFILE);
                writemod(modfile, taper_coeff, 3);

                MPI_Barrier(MPI_COMM_WORLD);
                if (MYID == 0)
                        mergemod(modfile, 3);
        }
}
