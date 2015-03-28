/*------------------------------------------------------------------------
 *   merge snapshots files written by the different processes to 
 *   a single file                                 
 *   last update 25/05/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void merge(int nsnap, int type){


extern char SNAP_FILE[STRING_SIZE];
extern int NXG, NYG, SNAP_FORMAT, NPROCX, NPROCY;
extern int NX, NY, IDX, IDY;
extern FILE *FP;
extern float DH;

char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10];
FILE *fp[NPROCY][NPROCX], *fpout;
int i, j, ip, jp, n;
float a, max;





switch(SNAP_FORMAT){
case 1: 
   sprintf(ext,".su");
   break;
case 2: 
   sprintf(ext,".asc");
   break;
case 3: 
   sprintf(ext,".bin");
   break;
}

switch(type){
case 1: 
   fprintf(FP," x-component of particle velocity");
   strcat(ext,".x");
   break;
case 2: 
   fprintf(FP," y-component of particle velocity");
   strcat(ext,".y");
   break;
case 4: 
   fprintf(FP," P-wave energyfield");
   strcat(ext,".div");
   break;
case 5: 
   fprintf(FP," S-wave energyfield");
   strcat(ext,".rot");
   break;
case 6: 
   fprintf(FP," pressure");
   strcat(ext,".p");
   break;
default: 
   err(" merge: cannot find snapfiles! ");
   break;
}

sprintf(mfile,"%s%s",SNAP_FILE,ext);
fprintf(FP," (files: %s.??? ).\n",mfile);

sprintf(outfile,"%s%s",SNAP_FILE,ext);
fprintf(FP,"\n writing merged snapshot file to  %s \n",outfile);
fpout=fopen(outfile,"w");



fprintf(FP," Opening snapshot files: %s.??? ",mfile);


for (ip=0;ip<=NPROCX-1; ip++)
   for (jp=0;jp<=NPROCY-1; jp++){
      sprintf(file,"%s.%i%i",mfile,ip,jp);
      fp[jp][ip]=fopen(file,"r");
      if (fp[jp][ip]==NULL) err("merge: can't read snapfile !"); 
      }

fprintf(FP," ... finished. \n");


fprintf(FP," Copying...");

max=0.0;
for (n=0;n<=nsnap-2; n++){
   for (ip=0;ip<=NPROCX-1; ip++){
      for (i=1;i<=NX;i+=IDX){
	 for (jp=0;jp<=NPROCY-1; jp++){
	    for (j=1;j<=NY;j+=IDY){
	       a=readdsk(fp[jp][ip],SNAP_FORMAT);
	       if (a>max) max=a;
	       writedsk(fpout,a,SNAP_FORMAT);
	       }
	    }
	 }
      }
   }
fprintf(FP," ... finished. \n");

for (ip=0;ip<=NPROCX-1; ip++)
   for (jp=0;jp<=NPROCY-1; jp++){
      fclose(fp[jp][ip]);
      }


if (SNAP_FORMAT==3){
   fprintf(FP," Use \n");
   fprintf(FP," xmovie n1=%d n2=%d  < %s loop=1 label1=Y label2=X title=%%g d1=%f d2=%f f1=%f f2=%f clip=%e \n",
      		((NYG-1)/IDY)+1,((NXG-1)/IDY)+1,outfile,DH,DH,DH,DH,max/10.0);
   fprintf(FP," to play the movie. \n");
   }


}


