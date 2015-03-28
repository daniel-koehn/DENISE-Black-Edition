/*------------------------------------------------------------------------
 *   Write note to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void note(FILE *fp){

extern char LOG_FILE[STRING_SIZE];
extern int MYID, LOG;

	fprintf(fp," Please note: \n");
	fprintf(fp," Each processing element (PE) is writing log information during program \n");
	fprintf(fp," execution to %s.PE .\n",LOG_FILE);
	fprintf(fp," See corresponding log-files for further information on program status.\n");
	fprintf(fp," Information about overall program execution \n");
	fprintf(fp," (numerical artefacts, accuracy, computing times etc) \n");
	fprintf(fp," will be written by PE 0 to ");
	if (LOG==1) fprintf(fp," standard output. \n");
	    else    fprintf(fp," %s.%i .\n",LOG_FILE,MYID);
}
