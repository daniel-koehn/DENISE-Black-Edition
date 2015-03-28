/*------------------------------------------------------------------------
 *   Write seismograms to disk                                  
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  inseis(FILE *fp, int comp, float **section, int ntr, int ns, int sws){

	/* declaration of extern variables */
	extern int NDT, MYID;
	extern char DATA_DIR[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern float  TIME, DH, DT, REFREC[4];
        char data[STRING_SIZE];
	const float xshift=800.0, yshift=800.0;
        FILE *fpdata;
	
	if(sws==1){  /* open seismic data */
	   sprintf(data,"%s.shot%d",DATA_DIR,comp);
	}
	
	if(sws==2){ /* open forward modelled data */
	   sprintf(data,"%s.shot%d.%d",SEIS_FILE_P,comp,MYID);
	}
	
	printf("%s\n",data);
	
	fpdata = fopen(data,"r");

	/* declaration of local variables */
	int i,j;
	segy tr;
	int tracl1 ;
	float xr, yr, x, y, dump;

		for(tracl1=1;tracl1<=ntr;tracl1++){        /* SEGY (without file-header) */
                        
			fread(&tr,240,1,fpdata);
			fread(&tr.data,4,ns,fpdata);
			
			section[tracl1][1]=0.0;
			
			for(j=1;j<=ns;j++){
			  dump=tr.data[j];
			  section[tracl1][j+1]=dump;
			  /*printf("%i \t %i \t %e \n",tracl1,j,section[tracl1][j]);*/
			}
		}



	fclose(fpdata);
}
