/*------------------------------------------------------------------------
 *   stress free surface condition
 *   last update 03/01/04, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_elastic(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy,
float ** sxy, float  **  pi, float  **  u, float ** rho, float * hc){


	int i,j,k;
	float fjm, g;
	float  vxx, vyy;
	float  dh24, dthalbe;	
	extern float DT, DH;
	extern int NX, INVMAT1;

	dthalbe=DT/2.0;
	dh24=1.0/DH;

	j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
	for (i=1;i<=NX;i++){

		/*Mirroring the components of the stress tensor to make
			a stress free surface (method of imaging)*/
		syy[j][i]=0.0;
		
		/*syy[j-1][i]=-syy[j+1][i];
		sxy[j-1][i]=-sxy[j][i];*/
		
		for (k=1;k<=4;k++){syy[j-k][i]=-syy[j+k][i];}
		for (k=1;k<=4;k++){sxy[j-k][i]=-sxy[j+k-1][i];} 
                
		/*vxx=(vx[j][i]-vx[j][i-1])*(dh24);*/
		vxx =  dh24*(hc[1]*(vx[j][i]-vx[j][i-1])
		     + hc[2]*(vx[j][i+1]-vx[j][i-2])
		     + hc[3]*(vx[j][i+2]-vx[j][i-3])
		     + hc[4]*(vx[j][i+3]-vx[j][i-4]));
		     
		/*vyy=(vy[j][i]-vy[j-1][i])*(dh24);*/
		vyy =  dh24*(hc[1]*(vy[j][i]-vy[j-1][i])
			   + hc[2]*(vy[j+1][i]-vy[j-2][i])  
			   + hc[3]*(vy[j+2][i]-vy[j-3][i])  
			   + hc[4]*(vy[j+3][i]-vy[j-4][i]));
		      
                if(INVMAT1==3){
		fjm=u[j][i]*2.0;
		g=pi[j][i];}
		
		if(INVMAT1==1){
		fjm=rho[j][i] * u[j][i] * u[j][i] * 2.0;
		g=rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}

		
		/*sxx[j][i]+= DT*((4.0*((g*fjm)+(fjm*fjm))/(g+2*fjm))*vxx);*/
		sxx[j][i]+= -DT*((g*g)/(g+fjm)*vxx+g*vyy);
		

	}
}
