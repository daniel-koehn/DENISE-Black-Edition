/*------------------------------------------------------------------------
 *   stress free surface condition
 *   T. Bohlen
 *   last update 2011/10/06, L. Groos
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_PML(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy,
float ** sxy, float ***p, float ***q, float  **  ppi, float  **  pu, float **prho, float **ptaup, float **ptaus, float *etajm, float *peta, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx){


	int i,j,m,h,h1,l;
	int fdoh;
	float bjm, djm, e, fjm, g;
	float  vxx, vyy, sump=0.0;
	float  dh24, dthalbe;	
	float *pts, ws, sumu, sumpi, mu, pi;
	extern float DT, DH, *FL;
	extern int NX, INVMAT1, L;
        extern int FW, BOUNDARY;
        extern int NPROCX, NPROCY, POS[3], MYID; 
	extern int FDORDER;
	extern float TS;
	
	
	fdoh = FDORDER/2;
	dthalbe=DT/2.0;
	dh24=1.0/DH;
	
	
	
	/* vector for maxwellbodies */
		pts=vector(1,L);
		for (l=1;l<=L;l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
		}
	
	
		/*ws=2.0*PI*FL[1];*/
		ws=2.0*PI*(1.0/TS);
	
		sumu=0.0;
		sumpi=0.0;
		for (l=1;l<=L;l++){
			sumu=sumu+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
			sumpi=sumpi+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
		}		
		
	

	j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
	for (i=1;i<=NX;i++){
		
		
		for (l=1;l<=L;l++){
			etajm[l]=peta[l];
		}
		
		
		
		/*Mirroring the components of the stress tensor to make
			a stress free surface (method of imaging)*/
		syy[j][i]=0.0;
		
		/* since syy is zero on the free surface also the
		corresponding memory-variables must set to zero */
		for (l=1;l<=L;l++) q[j][i][l]=0.0;
		
		
		
		/* now updating the stress component sxx and the memory-
		variables p[j][i][l] at the free surface */

		/* first calculate spatial derivatives of components
			of particle velocities */
		
		vxx = 0.0;
		vyy = 0.0;
		for (m=1; m<=fdoh; m++) {
		
			/*Mirroring the components of the stress tensor to make
			a stress free surface (method of imaging)*/
			syy[j-m][i]=-syy[j+m][i];
			sxy[j-m][i]=-sxy[j+m-1][i];
		
			vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]);
			vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]);
		}
		vxx *= dh24;
		vyy *= dh24;
		
		
		/*for (k=1;k<=4;k++){syy[j-k][i]=-syy[j+k][i];}
		for (k=1;k<=4;k++){sxy[j-k][i]=-sxy[j+k-1][i];}*/ 
                
		/*vxx=(vx[j][i]-vx[j][i-1])*(dh24);*/
		/*vxx =  dh24*(hc[1]*(vx[j][i]-vx[j][i-1])
		     + hc[2]*(vx[j][i+1]-vx[j][i-2])
		     + hc[3]*(vx[j][i+2]-vx[j][i-3])
		     + hc[4]*(vx[j][i+3]-vx[j][i-4]));*/
		     
		/*vyy=(vy[j][i]-vy[j-1][i])*(dh24);*/
		/*vyy =  dh24*(hc[1]*(vy[j][i]-vy[j-1][i])
			   + hc[2]*(vy[j+1][i]-vy[j-2][i])  
			   + hc[3]*(vy[j+2][i]-vy[j-3][i])  
			   + hc[4]*(vy[j+3][i]-vy[j-4][i]));*/

             
	     
             /* apply PML boundary */    
             /* left boundary */
             if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];                 
             }

             /* right boundary */
             if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=NX-FW+1)){
                
                        h1 = (i-NX+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1];                                            
             } 
	     
	     
	     
		
		
		
			
		if (INVMAT1==1){
			mu=(pu[j][i]*pu[j][i]*prho[j][i])/(1.0+sumu*ptaus[j][i]);
			pi=(ppi[j][i]*ppi[j][i]*prho[j][i])/(1.0+sumpi*ptaup[j][i]);
		}
		if (INVMAT1==3){
			mu=pu[j][i]/(1.0+sumu*ptaus[j][i]);
			pi=(ppi[j][i]+2*pu[j][i])/(1.0+sumpi*ptaup[j][i]);
		}
			     
	     
	     
	     
	     /* sums used in updating sxx */
		sump=0.0;
		for (l=1;l<=L;l++) sump+=p[j][i][l];


		fjm=mu*2.0*(1.0+L*ptaus[j][i]);
		g=pi*(1.0+L*ptaup[j][i]);

		/* partially updating sxx */
		sxx[j][i]+= -(DT*(g-fjm)*(g-fjm)*vxx/g)-(DT*(g-fjm)*vyy)-(dthalbe*sump);

		/* updating the memory-variable p[j][i][l] at the free surface */
		sump=0.0;
		for (l=1;l<=L;l++){
			bjm=etajm[l]/(1.0+(etajm[l]*0.5));
			djm=2.0*mu*ptaus[j][i];
			e=pi*ptaup[j][i];
			p[j][i][l]+=bjm*(((djm-e)*((fjm/g)-1.0)*vxx)-((djm-e)*vyy));
			sump+=p[j][i][l];
		}
		/*completely updating the stress sxx */
		sxx[j][i]+=(dthalbe*sump);   
		
	}
	free_vector(pts,1,L);
}
