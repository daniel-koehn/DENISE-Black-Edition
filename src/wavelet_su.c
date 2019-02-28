/*------------------------------------------------------------------------
 *   Read multiple source wavelets from SU file                                  
 *   
 *   Kiel, 17.02.2018
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  wavelet_su(int comp, float **section, int nsrc, int ns, int nsrc_loc, float ** srcpos_loc){

	/* declaration of extern variables */	
	extern char SIGNAL_FILE[STRING_SIZE];
        char data[STRING_SIZE];	
	float ** signals;
        FILE *fpdata;
		
	sprintf(data,"%s_shot%d.su",SIGNAL_FILE,comp);		
	
	signals = fmatrix(1,nsrc,1,ns);
	
	/*printf("%s\n",data);*/
	
	fpdata = fopen(data,"r");
        if (fpdata==NULL) err(" Source wavelet file could not be opened !");
        
	/* declaration of local variables */
	int i,j,l,pick_src;
	segy tr;
	float xr, yr, x, y, dump;

	/* read all source wavelets */
	for(l=1;l<=nsrc;l++){        /* SEGY (without file-header) */
                        
		fread(&tr,240,1,fpdata);
		fread(&tr.data,4,ns,fpdata);
			
		signals[l][1]=0.0;
			
			  
		for(j=1;j<=ns;j++){
		    dump=tr.data[j];			    			    
	            signals[l][j+1]=dump;
		}
			  
			
	}

	fclose(fpdata);
	
	for (l=1;l<=nsrc_loc;l++) {	
	
		pick_src = (int)srcpos_loc[3][l];
		/* printf("pick_src = %d of %d shots \n",pick_src,nsrc); */
	
	        for(j=1;j<=ns;j++){			  
		    section[l][j] = signals[pick_src][j];
	        }
		
	}
			
	free_matrix(signals,1,nsrc,1,ns);
	
}
