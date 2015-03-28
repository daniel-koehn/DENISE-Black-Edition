/*------------------------------------------------------------------------
 *  Read 2D gravity station coordinates  
 *
 *  D. Koehn
 *  Kiel, the 7th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float **read_grav_pos(int *ngrav){

	/* declaration of extern variables */
	extern char GRAV_STAT_POS[STRING_SIZE];
        extern int MYID, NGRAVB;
        extern float DH;

        /* local variables */
	float **recpos;
	int   itr=1, itr1=0, c;
	float xrec, yrec;
        FILE *fpr;

     	/* read receiver positions from file */
     	
        if(MYID==0){
           printf("\n Reading gravity station positions from file: \n\t%s\n",GRAV_STAT_POS);
        }	
 
        fpr=fopen(GRAV_STAT_POS,"r");
     	if ((fpr==NULL)&&(MYID==0)) err(" Gravity modelling station file could not be opened !");

	*ngrav=0;
	while ((c=fgetc(fpr)) != EOF)
		if (c=='\n') ++(*ngrav);
	rewind(fpr);

	recpos=matrix(1,2,1,*ngrav);
	for (itr=1;itr<=*ngrav;itr++){
		fscanf(fpr,"%f%f\n",&xrec, &yrec);
		recpos[1][itr]=xrec + (DH*NGRAVB);
		recpos[2][itr]=yrec;
	}
	fclose(fpr);

        if(MYID==0){
	  printf(" Message from function read_grav_pos (written by PE %d):\n",MYID);/***/
	  printf(" Number of gravity stations found: %i\n",*ngrav);
        }
 
       if(MYID==0){
     	   printf("Positions of gravity stations (x,y):\n");
     	   for (itr=1;itr<=*ngrav;itr++){
     	      printf("%e\t%e\n",recpos[1][itr],recpos[2][itr]);
           }
       }

       return recpos;
}
