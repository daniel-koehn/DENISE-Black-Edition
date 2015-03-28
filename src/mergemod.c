/*------------------------------------------------------------------------
 *   merge model files written by the different processes to 
 *   a single file                                 
 *   last update 08/10/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void mergemod(char modfile[STRING_SIZE], int format){



	extern int  NYG, MYID, NPROCX, NPROCY;
	extern int NX, NY, NPROC, IDX, IDY;
	extern FILE *FP;


	char file[STRING_SIZE];
	FILE *fp[NPROCY][NPROCX], *fpout;
	int i, j, ip, jp;
	float a;


	fprintf(FP,"\n **Message from mergemod (printed by PE %d):\n",MYID);
	fprintf(FP," PE %d starts merge of %d model files \n",MYID,NPROC);	

	fprintf(FP,"\n writing merged model file to  %s \n",modfile);
	fpout=fopen(modfile,"w");

	
	
	
	fprintf(FP," Opening model files: %s.??? ",modfile);
	for (ip=0;ip<=NPROCX-1; ip++)
  	for (jp=0;jp<=NPROCY-1; jp++){
      		sprintf(file,"%s.%i%i",modfile,ip,jp);
      		fp[jp][ip]=fopen(file,"r");
      		if (fp[jp][ip]==NULL) err("merge: can't read model file !"); 
      	}

	fprintf(FP," ... finished. \n");



	fprintf(FP," Copying...");

  	 for (ip=0;ip<=NPROCX-1; ip++){
      		for (i=1;i<=NX;i+=IDX){
			for (jp=0;jp<=NPROCY-1; jp++){
	    			for (j=1;j<=NY;j+=IDY){
	      			a=readdsk(fp[jp][ip],format);
	      			writedsk(fpout,a,format);
	       			}
	   		}
	 	}
      	}
	fprintf(FP," ... finished. \n");

	for (ip=0;ip<=NPROCX-1; ip++)
   	for (jp=0;jp<=NPROCY-1; jp++){
      		fclose(fp[jp][ip]);
      	}
	fclose(fpout);
	
	fprintf(FP," Use \n");
	fprintf(FP," ximage n1=%d < %s  label1=Y label2=X title=%s \n",
      			((NYG-1)/IDY)+1,modfile,modfile);
	fprintf(FP," to visualize model. \n");



}


