/*  Update_s_elastic_PML_PSV
 *
 *  updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *  by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *  and second order accuracy in time for the elastic PSV problem
 *   
 *  Daniel Koehn
 *  Kiel, 24.07.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_PML_PSV(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float **rho, float *hc, int infoout,
        float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int mode){


	int i,j, m, fdoh, h, h1;
	float fipjp, f, g;
	float  vxx, vyy, vxy, vyx;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, INVMAT1, FW;
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
	


	switch (FDORDER){

	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i]))*dhi; 

        /* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 			      			      
			      
			      			      				
			      sxy[j][i] += fipjp*(vyx+vxy);
			      sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
			      syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);
			}
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i]))*dhi; 

        /* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 
                              			      			      
				
				sxy[j][i] += fipjp*(vyx+vxy);
				sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
				syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);
			}
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i]))*dhi; 

        /* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
			
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 
			      			      				
				sxy[j][i] += fipjp*(vyx+vxy);
				sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
				syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);
			}
		}
		break;

	case 8:
/*#pragma acc enter data \
copyin(rho[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(uipjp[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(uxy[ny1:ny2-ny1+1][nx1:nx2-nx1+1],ux[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(vy[ny1-4:ny2-ny1+8][nx1-3:nx2-nx1+8],vx[ny1-3:ny2-ny1+8][nx1-4:nx2-nx1+8])\
copyin(uy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(u[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(hc[1:4],a_x_half[nx1:nx2-nx1+1])\
copyin(pi[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(psi_vyx[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(b_x_half[nx1:nx2-nx1+1],a_y[ny1:ny2-ny1+1])\
copyin(psi_vxx[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(b_y[ny1:ny2-ny1+1],K_x_half[nx1:nx2-nx1+1],K_y[ny1:ny2-ny1+1],POS[1:2],a_y_half[ny1:ny2-ny1+1],a_x[nx1:nx2-nx1+1])\
copyin(psi_vxy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])\
copyin(b_x[nx1:nx2-nx1+1],K_x[nx1:nx2-nx1+1],K_y_half[ny1:ny2-ny1+1],b_y_half[ny1:ny2-ny1+1])\
copyin(psi_vyy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])*/


#pragma acc parallel private(i,j)
#pragma acc loop independent gang 
    for (j=ny1;j<=ny2;j++){
#pragma acc loop independent vector 
	for (i=nx1;i<=nx2;i++){

			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i]))*dhi; 

        /* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 
			      
			      				
				sxy[j][i] += fipjp*(vyx+vxy);
				sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
				syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);

   }}

/*
#pragma exit data                                   \
copyout(sxy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])          \
copyout(sxx[ny1:ny2-ny1+1][nx1:nx2-nx1+1])          \
copyout(syy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])          \
copyout(psi_vyx[ny1:ny2-ny1+1][nx1:nx2-nx1+1])      \
copyout(psi_vxx[ny1:ny2-ny1+1][nx1:nx2-nx1+1])      \
copyout(psi_vxy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])      \
*/
/*copyout(K_y_half[ny1:ny2-ny1+1])
copyout(b_y_half[ny1:ny2-ny1+1]) \ 
copyout(psi_vyy[ny1:ny2-ny1+1][nx1:nx2-nx1+1])
copyout(b_y[ny1:ny2-ny1+1],K_x_half[nx1:nx2-nx1+1],K_y[ny1:ny2-ny1+1],a_y_half[ny1:ny2-ny1+1],a_x[nx1:nx2-nx1+1]) 
copyout(b_x_half[nx1:nx2-nx1+1],a_y[ny1:ny2-ny1+1]) \
copyout(b_x[nx1:nx2-nx1+1],K_x[nx1:nx2-nx1+1]) */

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
				      
	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 
			      			      				
				sxy[j][i] += fipjp*(vyx+vxy);
				sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
				syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);
			
				      
			
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
				
	 /* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
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
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

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
                                  g = pi[j][i];}  
				
			      if (INVMAT1==1){
	                          f = rho[j][i] * u[j][i] * u[j][i];
                                  g = rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (g*(vxx+vyy)+(2.0*f*vxx))/DT;
				uy[j][i] = (g*(vxx+vyy)+(2.0*f*vyy))/DT;
				uxy[j][i] = fipjp*(vyx+vxy)/DT;
			      } 			      			      			      
				
				sxy[j][i] += fipjp*(vyx+vxy);
				sxx[j][i] += g*(vxx+vyy)+(2.0*f*vxx);
				syy[j][i] += g*(vxx+vyy)+(2.0*f*vyy);      
				
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
