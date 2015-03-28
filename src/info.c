/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ***********************************************************\n");
	fprintf(fp," This is program DENISE. Version 1.0                        \n");
	fprintf(fp," Parallel 2-D elastic Finite Difference FWT code            \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," FWT code written by D. Koehn                               \n");
	fprintf(fp," forward code written by  T. Bohlen                         \n");
	fprintf(fp," Institute of Geosciences, Kiel University, Germany         \n\n");
	fprintf(fp," See COPYING file for copying and redistribution conditions.\n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");

}
