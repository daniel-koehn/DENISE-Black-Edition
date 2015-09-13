/* $Id: update_s_elastic_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_PML(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float **rho, float *hc, int infoout,
        float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int sws){


	int i,j, m, fdoh, h, h1;
	float fipjp, f, g;
	float  vxx, vyy, vxy, vyx;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, INVMAT1, INVMAT, FW;
        extern int FREE_SURF, BOUNDARY, GRAD_FORM;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
		
	dhi = DT/DH;
	fdoh = FDORDER/2;
	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}

        for (j=ny1;j<=ny2;j++){
	    for (i=nx1;i<=nx2;i++){

               vxx = 0.0;
	       vyx = 0.0;
	       vxy = 0.0;
	       vyy = 0.0;

               /* calculate spatial derivatives */
               for (m=1;m<=fdoh;m++){
	     
	          vxx += hc[m] * (vx[j][i+m-1] - vx[j][i-m]  ) * dhi;
		  vyy += hc[m] * (vy[j+m-1][i] - vy[j-m][i]  ) * dhi;
		  vyx += hc[m] * (vy[j][i+m]   - vy[j][i-m+1]) * dhi;
		  vxy += hc[m] * (vx[j+m][i]   - vx[j-m+1][i]) * dhi;    
	     
	       }

               /* left boundary */                                         
               if ((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
                }

                /* right boundary */                                         
                if ((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
                 }

	         /* top boundary */                                         
                 if ((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

                  }
	
	          /* bottom boundary */                                         
                  if ((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
                   }

                   fipjp=uipjp[j][i];			      
			      
                   /* lambda - mu relationship*/
		   if (INVMAT1==3){
		      f = u[j][i];
                      g = pi[j][i];
		   }  
				
		   if (INVMAT1==1){
			
	              f = rho[j][i] * u[j][i] * u[j][i];
                      g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
		      /* save time derivative of forward wavefield for gradient calculation */
		      if ((INVMAT<=1)&&(GRAD_FORM==2)){
	                  ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
		          uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
			  uxy[j][i] = fipjp*(vyx+vxy)/DT;
		       } 
			      
		       if ((INVMAT<=1)&&(GRAD_FORM==3)){
			      
			   if (sws==0){
	                      ux[j][i] = vxx/DT;
			      uy[j][i] = vyy/DT;   
			      uxy[j][i] = (vyx+vxy)/DT;
		           }
				
			   if (sws==1){
	                      ux[j][i] += vxx;
			      uy[j][i] += vyy;
			      uxy[j][i] += 0.5*(vyx+vxy);
			   }
				  
		        }
			      
		        if ((INVMAT<=1)&&(GRAD_FORM==4)){
			      
			   if (sws==0){
	                      ux[j][i] = vxx;
			      uy[j][i] = vyy;
			      uxy[j][i] = vyx+vxy;
			   }
				  
			}
				
			sxy[j][i] += fipjp*(vyx+vxy);
			sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
			syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);

            }
	}
						

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
