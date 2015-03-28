/*------------------------------------------------------------------------
 *   calculate step length for material parameter update
 *   
 *   waveconv = conjugated gradient direction
 *   gradp    = preconditionated gradient
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
float calc_opt_step(float *  L2t, float ** waveconv, float ** gradg, float * epst, int sws, float C_vp){

extern int NX, NY, IDX, IDY, MYID;
extern float EPSILON, EPSILON_u, EPSILON_rho;
int i, j, n;
float opteps, H1, H2, critmult;
float *x, *b, **A;

critmult = 5.0;
n = 3; /* number of test calculations */

A =  matrix(1,n,1,n);
x = vector(1,n);
b = vector(1,n);

/* calculate optimal step size after Tarantola (1986)*/
/*H1=0.0;
H2=0.0;
for (i=1;i<=NX;i=i+IDX){
     for (j=1;j<=NY;j=j+IDY){
         H1 += waveconv[j][i]*(1.0/C_vp)*waveconv[j][i];
	 H2 += waveconv[j][i]*gradg[j][i];
     }
  }
    
H1=exchange_L2(H1,1,1); 
H2=exchange_L2(H2,1,1);*/

/* calculate optimal step length for Vp update */
/*if(sws==1){
opteps = (H2/((L2t[sws]*L2t[sws]/(EPSILON*EPSILON))+H1));
if(fabs(opteps) > (10.0 * fabs(EPSILON)) ){opteps=EPSILON;}
/*opteps = EPSILON;*/
/*}*/ 

/* calculate optimal step length for Vs update */
/*if(sws==2){
opteps = EPSILON_u * ((L2t[sws]*L2t[4])/(L2t[sws]*L2t[sws]));
if(fabs(opteps) > (10.0 * fabs(EPSILON_u)) ){opteps=EPSILON_u;}
} */

/* calculate optimal step size by line search */

/* fit parabola function to L2 norm */

/* define coefficient matrix A */
for (i=1;i<=n;i++){
   A[i][3]=(epst[i]*epst[i]);
   A[i][2]=(epst[i]);
   A[i][1]=(1.0);
}

/* define RHS vector b */
for (i=1;i<=n;i++){
   b[i]=(L2t[i]);
}

/* solve matrix equation using LU decomposition */
/*LU_decomp(A,x,b,n);*/
solvelin(A,b,x,n,1);

/* calculate optimal step length -> extremum of the parabola */
opteps = -x[2]/(2.0*x[3]);

/* if L2[1] < L2[2] < L2[3]*/
if (((2.0*x[3])< 0.0)&&(L2t[2] > L2t[3])){
opteps = epst[3];
}

/* if L2[1] > L2[2] > L2[3] */
/*if (((2.0*x[3])< 0.0)&&(L2t[2] > L2t[1])){
opteps = epst[3];
}*/

/* if opteps < 50.0 set opteps=50.0 */
if (opteps > epst[3]){
opteps = epst[3];
}

if (opteps < 0.0){
opteps = epst[1];
}

/*if (epst[3]==0){
opteps = epst[1]/2.0;
}*/

/*opteps = epst[1];*/

/*if (isnan(opteps))
{opteps = epst[1];}*/


free_matrix(A,1,n,1,n);
free_vector(x,1,n);
free_vector(b,1,n);
return opteps;		
}

