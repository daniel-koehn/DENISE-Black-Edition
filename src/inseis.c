/*------------------------------------------------------------------------
 *   Write seismograms to disk                                  
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  inseis(FILE *fp, int comp, float **section, int ntr, int ns, int sws, int iter){

	/* declaration of extern variables */
	extern int NDT, MYID;
	extern char DATA_DIR[STRING_SIZE], SEIS_FILE_P[STRING_SIZE], SEIS_FILE_VX[STRING_SIZE];
	extern char DATA_DIR_T0[STRING_SIZE];
	extern char SEIS_FILE_VY[STRING_SIZE];
	extern float  TIME, DH, DT, REFREC[4];
        char data[STRING_SIZE];
	const float xshift=800.0, yshift=800.0;
        FILE *fpdata;
	
	if(sws==1){  /* open seismic data vx */
	   sprintf(data,"%s_x.su.shot%d",DATA_DIR,comp);
	}
	
	if(sws==2){  /* open seismic data vy */
	   sprintf(data,"%s_y.su.shot%d",DATA_DIR,comp);
	}
	
	if(sws==3){ /* open forward modelled data vx*/
	   sprintf(data,"%s.shot%d.%d",SEIS_FILE_VX,comp,MYID);
	}
	
	if(sws==4){ /* open forward modelled data vy*/
	   sprintf(data,"%s.shot%d.%d",SEIS_FILE_VY,comp,MYID);
	}
	
	if(sws==5){ /* open old data residuals vx*/
	   sprintf(data,"%s.shot%d_it-1.%d",SEIS_FILE_VX,comp,MYID);
	}
	
	if(sws==6){ /* open old data residuals vy*/
	   sprintf(data,"%s.shot%d_it-1.%d",SEIS_FILE_VY,comp,MYID);
	}
	
	if(sws==7){  /* open convolved seismic data vx */
	   sprintf(data,"%s_x.su.conv.shot%d",DATA_DIR,comp);
	}
	
	if(sws==8){  /* open convolved seismic data vy */
	   sprintf(data,"%s_y.su.conv.shot%d",DATA_DIR,comp);
	}
	
	if(sws==9){  /* if TIMELAPSE == 1, open synthetic seismic data vx at t0 */
	   sprintf(data,"%s_x.su.shot%d",DATA_DIR_T0,comp);
	}
	
	if(sws==10){  /* if TIMELAPSE == 1, open synthetic seismic data vy at t0 */
	   sprintf(data,"%s_y.su.shot%d",DATA_DIR_T0,comp);
	}
	
	/*printf("%s\n",data);*/
	
	fpdata = fopen(data,"r");
        if (fpdata==NULL) err(" Field data file could not be opened !");
        
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
			  }
			  
			
		}



	fclose(fpdata);
}
