
#include "fd.h"

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy){

	int i,j;
	char myid[5];
	char checkptfile[STRING_SIZE];
	FILE *fp;

	extern int MYID;
	extern char  CHECKPTFILE[STRING_SIZE];



	sprintf(checkptfile,"%s",CHECKPTFILE);
	sprintf(myid,".%d",MYID);
	strcat(checkptfile,myid);



	fp=fopen(checkptfile,"wb");
	if (fp==NULL) err("CHECKPTFILE can't be opened !");
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){

		fwrite( &vx[j][i],sizeof(float),1,fp);
		fwrite( &vy[j][i],sizeof(float),1,fp);
		fwrite(&sxx[j][i],sizeof(float),1,fp);
		fwrite(&syy[j][i],sizeof(float),1,fp);
		fwrite(&sxy[j][i],sizeof(float),1,fp);
		}
	}


	fclose(fp);

}
