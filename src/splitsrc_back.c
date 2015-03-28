/*  Computation of local source coordinates
 * 
 *	 
 *   last update 12.02.02, T. Bohlen
 */
 
 #include "fd.h"
float **splitsrc_back(int **recpos,int *nsrc_loc, int nsrc)
{

	extern int IENDX, IENDY, POS[3], MYID, REC1, REC2;
	extern float DH;

	int a,b,i=0,j,k,found=0;
	float ** srcpos_dummy, **srcpos_local=NULL;
	srcpos_dummy = matrix(1,6,1,nsrc);

	for (j=1;j<=nsrc;j++) {
		a=(recpos[1][j]-1)/IENDX;
		b=(recpos[2][j]-1)/IENDY;


		if ((POS[1]==a)&&(POS[2]==b)) {
			found = 1;
			i++;
			srcpos_dummy[1][i] = ((recpos[1][j]-1)%IENDX)+1;
			srcpos_dummy[2][i] = ((recpos[2][j]-1)%IENDY)+1;
			srcpos_dummy[3][i] = 0.0;
			srcpos_dummy[4][i] = 0.0;
			srcpos_dummy[5][i] = 125.0;
			srcpos_dummy[6][i] = 1.0;
		}
	}
   
	if (i>0) srcpos_local = matrix(1,6,1,i);
	REC2=REC1+(i-1);
	for (k=1;k<=i;k++){
		srcpos_local[1][k] = srcpos_dummy[1][k];
		srcpos_local[2][k] = srcpos_dummy[2][k];
		srcpos_local[3][k] = srcpos_dummy[3][k];
		srcpos_local[4][k] = srcpos_dummy[4][k];
		srcpos_local[5][k] = srcpos_dummy[5][k];
		srcpos_local[6][k] = srcpos_dummy[6][k];
	}
	free_matrix(srcpos_dummy,1,6,1,nsrc);

	printf("\n **Message from splitsrc:\n");
	printf(" Splitting of receiver positions from global to local grids finished.\n");
	printf(" MYID= %d \t \t no. of sources= %d\n",MYID,i);
	

	printf("\n **Message from splitsrc:\n");
	printf(" Table of local source positions (in gridpoints), time-shift, centre frequency and amplitude:\n");
	printf(" MYID\t  x\t  y\t  z\t  tshift  fc\t  amp\n");

	for (j=1;j<=i;j++)
	    printf(" %3d\t%4.0f\t%4.0f\t%4.0f\t%6.2f\t%6.2f\t%6.2f\n",MYID,srcpos_local[1][j],srcpos_local[2][j],srcpos_local[3][j],srcpos_local[4][j],srcpos_local[5][j],srcpos_local[6][j]);

        
        *nsrc_loc=i;
	return srcpos_local;

}
