/* $Id: update_s_elastic_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_hc(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float **rho, float *hc, int infoout){


	int i,j, m, fdoh;
	float fipjp, f, g;
	float  vxx, vyy, vxy, vyx;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, INVMAT1;
	extern FILE *FP;
	double time1, time2;
	
	

	dhi = DT/DH;
	fdoh = FDORDER/2;

	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	


	switch (FDORDER){

	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				/* Compute values for shearmodulus u[j][i], 
				    P-wave modulus pi[j][i],
				    tau for S-waves and P-waves taus[j][i], 
				    taup[j][i] at staggered grid points: */
		       
				/* spatial derivatives of the components of the velocities */
				/* using Holberg coefficients */
				vxx = hc[1]*(vx[j][i]  -vx[j][i-1])*dhi;
				vyy = hc[1]*(vy[j][i]  -vy[j-1][i])*dhi;
				vyx = hc[1]*(vy[j][i+1]-vy[j][i])*dhi;
				vxy = hc[1]*(vx[j+1][i]-vx[j][i])*dhi;
		       
		                ux[j][i] = vxx;
				uy[j][i] = vyy;
				
				uxy[j][i] = (vxy + vyx);
				
				/* updating components of the stress tensor, partially */
				fipjp=uipjp[j][i];
				f=u[j][i];
				g=pi[j][i];

				sxy[j][i]=(fipjp*(vxy+vyx));
				sxx[j][i]=(g*(vxx+vyy))+(2.0*f*vxx);
				syy[j][i]=(g*(vxx+vyy))+(2.0*f*vyy);
			}
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				      )*dhi;
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				      )*dhi;
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				      )*dhi;
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				      )*dhi;
		                
				ux[j][i] += DT*vxx;
				uy[j][i] += DT*vyy;
		       
				fipjp=uipjp[j][i]*DT;
				f=u[j][i]*DT;
				g=pi[j][i]*DT;

				sxy[j][i]+=(fipjp*(vxy+vyx));
				sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
			}
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				      )*dhi;
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				      )*dhi;
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				      )*dhi;
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				      )*dhi;
				
				ux[j][i] += DT*vxx;
				uy[j][i] += DT*vyy;      
		       
				fipjp=uipjp[j][i]*DT;
				f=u[j][i]*DT;
				g=pi[j][i]*DT;

				sxy[j][i]+=(fipjp*(vxy+vyx));
				sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
			}
		}
		break;

	case 8:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4])
				      )*dhi;
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i])
				      )*dhi;
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3])
				      )*dhi;
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i])
				      )*dhi;
                                
				/* save sum of strain tensor exx + eyy */
		                ux[j][i] += vxx;
				uy[j][i] += vyy;
				
				uxy[j][i] += (vxy + vyx);
		       
				fipjp=uipjp[j][i];
				
				/* lambda - mu relationship*/
				if (INVMAT1==3){
				f=u[j][i];
                                g = pi[j][i];}  
				
				if (INVMAT1==1){
				f = rho[j][i] * u[j][i] * u[j][i];
                                g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
				sxy[j][i] += fipjp*(vxy+vyx);
				sxx[j][i] += (g*(vxx+vyy))+(2.0*f*vxx);
				syy[j][i] += (g*(vxx+vyy))+(2.0*f*vyy);
			}
		}
		break;

	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4])
				       + hc[5]*(vx[j][i+4]-vx[j][i-5])
				      )*dhi;
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i])
				       + hc[5]*(vy[j+4][i]-vy[j-5][i])
				      )*dhi;
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3])
				       + hc[5]*(vy[j][i+5]-vy[j][i-4])
				      )*dhi;
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i])
				       + hc[5]*(vx[j+5][i]-vx[j-4][i])
				      )*dhi;
				      
				ux[j][i] += DT*vxx;
				uy[j][i] += DT*vyy;      
		       		
				fipjp=uipjp[j][i]*DT;
				f=u[j][i]*DT;
				g=pi[j][i]*DT;
				
				sxy[j][i]+=(fipjp*(vxy+vyx));
				sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				
				/*
				piv = pi[j][i]*(vxx+vyy);
				sxy[j][i]+=(uipjp[j][i]*(vxy+vyx));
				sxx[j][i]+=piv-(2.0*u[j][i]*vyy);
				syy[j][i]+=piv-(2.0*u[j][i]*vxx);
				*/
			}
		}
		break;
		
	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4])
				       + hc[5]*(vx[j][i+4]-vx[j][i-5])
				       + hc[6]*(vx[j][i+5]-vx[j][i-6])
				      )*dhi;
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i])
				       + hc[5]*(vy[j+4][i]-vy[j-5][i])
				       + hc[6]*(vy[j+5][i]-vy[j-6][i])
				      )*dhi;
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3])
				       + hc[5]*(vy[j][i+5]-vy[j][i-4])
				       + hc[6]*(vy[j][i+6]-vy[j][i-5])
				      )*dhi;
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i])
				       + hc[5]*(vx[j+5][i]-vx[j-4][i])
				       + hc[6]*(vx[j+6][i]-vx[j-5][i])
				      )*dhi;
				      
				ux[j][i] += DT*vxx;
				uy[j][i] += DT*vyy;
		       
				fipjp=uipjp[j][i]*DT;
				f=u[j][i]*DT;
				g=pi[j][i]*DT;

				sxy[j][i]+=(fipjp*(vxy+vyx));
				sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
			}
		}
		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = 0.0;
				vyy = 0.0;
				vyx = 0.0;
				vxy = 0.0;
				for (m=1; m<=fdoh; m++) {
					vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
					vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
					vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
					vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
				}	
	
				fipjp=uipjp[j][i]*DT;
				f=u[j][i]*DT;
				g=pi[j][i]*DT;	

				sxy[j][i]+=(fipjp*(vxy+vyx))*dhi;
				sxx[j][i]+=((g*(vxx+vyy))-(2.0*f*vyy))*dhi;
				syy[j][i]+=((g*(vxx+vyy))-(2.0*f*vxx))*dhi;
			}
		}
		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
