/*------------------------------------------------------------------------
 *  compute receiver positions or read them from external file        
 *  last update: 11.02.02
 *
 *  T. Bohlen
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int **receiver(FILE *fp, int *ntr, int ishot){

	/* declaration of extern variables */
	extern  char REC_FILE[STRING_SIZE];
	extern float DH, REFREC[4], FW;
	extern int READREC, NX;
	extern int MYID_SHOT, N_STREAMER;
	extern float REC_INCR_X, REC_INCR_Y;

	int **recpos1, **recpos, nxrec=0, nyrec=0, nzrec=0;
	int   itr=1, itr1=0, itr2=0, recflag=0, c, ifw, n, i, j;
	int nxrec1, nxrec2, nyrec1, nyrec2, nzrec1=1, nzrec2=1;
	float xrec, yrec, aa, bb;
	char REC_FILE1[STRING_SIZE];
	
	FILE *fpr;

	if (MYID_SHOT==0)
	{

     	     if (READREC){ /* read receiver positions from file */

		if(READREC==1){
		    sprintf(REC_FILE1,"%s.dat",REC_FILE);
     		    fprintf(fp,"\n Reading receiver positions from single file: \n\t%s\n",REC_FILE1);
		}

		if(READREC==2){
		    sprintf(REC_FILE1,"%s_shot_%i.dat",REC_FILE,ishot);
     		    fprintf(fp,"\n Reading receiver positions from multiple files: \n\t%s\n",REC_FILE1);
		}

		fpr=fopen(REC_FILE1,"r");
     		if (fpr==NULL) err(" Receiver file could not be opened !");
     		*ntr=0;
     		while ((c=fgetc(fpr)) != EOF)
     			if (c=='\n') ++(*ntr);
     		rewind(fpr);
     
     		recpos1=imatrix(1,3,1,*ntr);
		
		aa = (float)(ishot-1)*REC_INCR_X;
                bb = (float)(ishot-1)*REC_INCR_Y;
		
     		for (itr=1;itr<=*ntr;itr++){
     			fscanf(fpr,"%f%f\n",&xrec, &yrec);
     			
			if(itr<=N_STREAMER){
			  recpos1[1][itr]=iround((xrec+REFREC[1]+aa)/DH);
     			  recpos1[2][itr]=iround((yrec+REFREC[2]+bb)/DH);
			}else{
			  recpos1[1][itr]=iround((xrec+REFREC[1])/DH);
     			  recpos1[2][itr]=iround((yrec+REFREC[2])/DH);
			}
     			recpos1[3][itr]=iround((0.0+REFREC[3])/DH);
     		}
     		fclose(fpr);
     		fprintf(fp," Message from function receiver (written by PE %d):\n",MYID_SHOT);/***/
     		fprintf(fp," Number of receiver positions found: %i\n",*ntr);
     
     		/* check if more than one receiver is located
     				         at the same gridpoint */
     		for (itr=1;itr<=(*ntr-1);itr++)
     			for (itr1=itr+1;itr1<=*ntr;itr1++)
     				if ((recpos1[1][itr]==recpos1[1][itr1])
     				    && (recpos1[2][itr]==recpos1[2][itr1])
     				    && (recpos1[3][itr]==recpos1[3][itr1]))
     					recpos1[1][itr1]=-(++recflag);
     
     		recpos=imatrix(1,3,1,*ntr-recflag);
     		for (itr=1;itr<=*ntr;itr++)
     			if (recpos1[1][itr]>0){
     				recpos[1][++itr2]=recpos1[1][itr];
     				recpos[2][itr2]=recpos1[2][itr];
     				recpos[3][itr2]=recpos1[3][itr];
     			}
     
     		*ntr=itr2;
     		if (recflag>0){
     			fprintf(fp,"\n\n");
     			fprintf(fp," Warning:\n");
     			fprintf(fp," Several receivers located at the same gridpoint !\n");
     			fprintf(fp," Number of receivers reduced to %i\n", *ntr);
     			fprintf(fp,"\n\n");
     		}
     
     		free_imatrix(recpos1,1,3,1,*ntr);
  
     	}
     
     	/*   fprintf(fp,"Gridpoints of receiver positions (x,y,z):\n");
     		for (itr=1;itr<=*ntr;itr++)
     				fprintf(fp,"%i\t%i\t%i\n",recpos[1][itr],recpos[2][itr],recpos[3][itr]);*/
     
    

	}


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);
	if (MYID_SHOT!=0) recpos=imatrix(1,3,1,*ntr);
	MPI_Bcast(&recpos[1][1],(*ntr)*3,MPI_INT,0,MPI_COMM_WORLD);

/*	if (MYID_SHOT==0)
	{
		fprintf(fp,"\n **Message from function receiver (written by PE %d):\n",MYID_SHOT);
		fprintf(fp," Number of receiver positions found: %i\n",*ntr);
		fprintf(fp," Receiver positions (in gridpoints) in the global model-system:\n");
		fprintf(fp," x  \ty \n");
		fprintf(fp," -  \t- \n");
		for (l=1;l<=*ntr;l++)
			fprintf(fp," %i\t%i\n",recpos[1][l],recpos[2][l]);
		fprintf(fp,"\n\n");
	}
*/

	return recpos;
}
