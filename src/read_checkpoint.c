
#include "fd.h"

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy){

	int i,j;
	char myid[5];
	FILE *fp;
	char checkptfile[STRING_SIZE];
	extern int MYID;
	extern char  CHECKPTFILE[STRING_SIZE];



	sprintf(checkptfile,"%s",CHECKPTFILE);
	sprintf(myid,".%d",MYID);
	strcat(checkptfile,myid);



	fp=fopen(checkptfile,"rb");
	if (fp==NULL) err("CHECKPTFILE can't be opened !");
	
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){

		fread( &vx[j][i],sizeof(float),1,fp);
		fread( &vy[j][i],sizeof(float),1,fp);
		fread(&sxx[j][i],sizeof(float),1,fp);
		fread(&syy[j][i],sizeof(float),1,fp);
		fread(&sxy[j][i],sizeof(float),1,fp);
		}
	}

	fclose(fp);

}
