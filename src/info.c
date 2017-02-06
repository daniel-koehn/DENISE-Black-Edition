/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," *******************************************************************************\n");
	fprintf(fp," This is program DENISE Black-Edition                                           \n");
	fprintf(fp," Parallel 2-D elastic Finite Difference FWI code                                \n");
	fprintf(fp,"                                                                                \n");
	fprintf(fp," Forward/FWI/RTM codes written by D. Koehn and D. De Nil                        \n");
	fprintf(fp," 2D isotropic PSV forward code partly based on FDVEPS written by T. Bohlen      \n");
	fprintf(fp," Institute of Geosciences, Kiel University, Germany                           \n\n");
	fprintf(fp," See README.md file and LICENSE.md for redistribution conditions.               \n");
	fprintf(fp," *******************************************************************************\n");
	fprintf(fp,"\n");

}
