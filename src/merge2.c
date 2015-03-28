/*------------------------------------------------------------------------
 *   merge snapshots files written at different timesteps 
 *   to a single file                                 
 *   last update 05/04/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void merge2(int nsnap, int type){



	char file[STRING_SIZE], mfile[STRING_SIZE], ext[10];
	FILE *fpw, *fpr;
	int i, j, n, n1, n2;
	float amp;


	extern char SNAP_FILE[STRING_SIZE];
	extern int NXG, NYG, SNAP_FORMAT;
	extern FILE *FP;



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


	fprintf(FP,"\n **Message from merge:\n");
	fprintf(FP," Start merge of %d snapfiles of \n",nsnap);	

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

	fprintf(FP," Reading: ");

	fpw=fopen(mfile,"w");
	for (n=1;n<=nsnap; n++){
		sprintf(file,"%s.%i",mfile,n);
		fprintf(FP,"%s  ",file);
		fpr=fopen(file,"r");
		if (fpr==NULL) err("merge: can't read snapfile !"); 
		for (i=1;i<=NXG;i++)
		for (j=1;j<=NYG;j++){
			amp=readdsk(fpr,SNAP_FORMAT);
			writedsk(fpw,amp,SNAP_FORMAT);
		}
		fclose(fpr);
	}
	fclose(fpw);
	
	fprintf(FP,"\n merged snapshot file written to  %s \n",mfile);


	n1=NYG;
	n2=NXG;

	if (SNAP_FORMAT==3){
		fprintf(FP," Use \n");
		fprintf(FP," xmovie n1=%d n2=%d < %s loop=1 label1=Y label2=X title=%%g \n",
			n1,n2,mfile);
		fprintf(FP," to visualize snapshot. \n");
	}


}


