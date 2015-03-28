/*------------------------------------------------------------------------
 *   Read extern source wavelet                                26 Mar 1997
 *
 *  written by T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float *rd_sour(int *nts,FILE* fp_source){

	/* local variables */
	float *psource;
	int i, c;

	if (fp_source==NULL) err(" Source file could no be opened !");
	/* fscanf(fp_source,"%i", nts); */
        *nts=0;
        while ((c=fgetc(fp_source)) != EOF)
         if (c=='\n') ++(*nts);
        rewind(fp_source);
			/*printf(" Number of samples (nts) in source file: %i\n",*nts);*/
		  
	psource=vector(1,*nts);
	for (i=1;i<=*nts;i++) fscanf(fp_source,"%e\n",&psource[i]);
	fclose(fp_source);
	return psource;
}
