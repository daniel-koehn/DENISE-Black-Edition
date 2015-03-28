/* $Id: update_s_visc_PML.c,v 1.1.1.1 2011/10/06 22:44:52 groos Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_PML(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** pi, float ** u, float ** uipjp, float **rho, float *hc, int infoout,
	float ***r, float ***p, float ***q, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, 
      float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
      float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
      float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx){



	int i,j, m, fdoh, h, h1, l;
	float  vxx, vyy, vxy, vyx;
	float  dhi, dthalbe;	
	extern float DT, DH;
	extern int MYID, FDORDER, FW, L;
        extern int FREE_SURF, BOUNDARY, INVMAT;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
	
	float sumr=0.0, sump=0.0, sumq=0.0;
	
	

	/*dhi = DT/DH;*/
	dhi=1.0/DH;
	fdoh = FDORDER/2;
	dthalbe = DT/2.0;

	
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
	
	/* computing sums of the old memory variables */
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				sumr+=r[j][i][l];
				sump+=p[j][i][l];
				sumq+=q[j][i][l];
			}
			
			
                        /* updating components of the stress tensor, partially */
			sxy[j][i] += (fipjp[j][i]*(vxy+vyx))+(dthalbe*sumr);
			sxx[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vyy)+(dthalbe*sump);
			syy[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vxx)+(dthalbe*sumq);
				
			
			/* now updating the memory-variables and sum them up*/
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*(vxy+vyx)));
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vyy));
				q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vxx));
				sumr += r[j][i][l];
				sump += p[j][i][l];
				sumq += q[j][i][l];
			}
			
			
			/* and now the components of the stress tensor are
			   completely updated */
			sxy[j][i]+=(dthalbe*sumr);
			sxx[j][i]+=(dthalbe*sump);
			syy[j][i]+=(dthalbe*sumq);

			/* save forward wavefield for gradient calculation */
			if(INVMAT==0){
			  ux[j][i] = vxx;
		          uy[j][i] = vyy;
			  uxy[j][i] = vxy + vyx;
			}

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

	/* computing sums of the old memory variables */
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				sumr+=r[j][i][l];
				sump+=p[j][i][l];
				sumq+=q[j][i][l];
			}


                        /* updating components of the stress tensor, partially */
			sxy[j][i] += (fipjp[j][i]*(vxy+vyx))+(dthalbe*sumr);
			sxx[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vyy)+(dthalbe*sump);
			syy[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vxx)+(dthalbe*sumq);
			
			
			/* now updating the memory-variables and sum them up*/
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*(vxy+vyx)));
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vyy));
				q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vxx));
				sumr += r[j][i][l];
				sump += p[j][i][l];
				sumq += q[j][i][l];
			}

			/* and now the components of the stress tensor are
			   completely updated */
			sxy[j][i]+=(dthalbe*sumr);
			sxx[j][i]+=(dthalbe*sump);
			syy[j][i]+=(dthalbe*sumq);

			/* save forward wavefield for gradient calculation */
			if(INVMAT==0){
			  ux[j][i] = vxx;
		          uy[j][i] = vyy;
			  uxy[j][i] = vxy + vyx;
			}

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

	/* computing sums of the old memory variables */
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				sumr+=r[j][i][l];
				sump+=p[j][i][l];
				sumq+=q[j][i][l];
			}


                        /* updating components of the stress tensor, partially */
			sxy[j][i] += (fipjp[j][i]*(vxy+vyx))+(dthalbe*sumr);
			sxx[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vyy)+(dthalbe*sump);
			syy[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vxx)+(dthalbe*sumq);
			
			
			/* now updating the memory-variables and sum them up*/
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*(vxy+vyx)));
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vyy));
				q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vxx));
				sumr += r[j][i][l];
				sump += p[j][i][l];
				sumq += q[j][i][l];
			}

			/* and now the components of the stress tensor are
			   completely updated */
			sxy[j][i]+=(dthalbe*sumr);
			sxx[j][i]+=(dthalbe*sump);
			syy[j][i]+=(dthalbe*sumq);

			/* save forward wavefield for gradient calculation */
			if(INVMAT==0){
			  ux[j][i] = vxx;
		          uy[j][i] = vyy;
			  uxy[j][i] = vxy + vyx;
			}

			}
		}
		break;

	case 8:

    for (j=ny1;j<=ny2;j++){
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

	/* computing sums of the old memory variables */
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				sumr+=r[j][i][l];
				sump+=p[j][i][l];
				sumq+=q[j][i][l];
			}


                        /* updating components of the stress tensor, partially */
			sxy[j][i] += (fipjp[j][i]*(vxy+vyx))+(dthalbe*sumr);
			sxx[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vyy)+(dthalbe*sump);
			syy[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vxx)+(dthalbe*sumq);
			
			
			/* now updating the memory-variables and sum them up*/
			sumr=sump=sumq=0.0;
			for (l=1;l<=L;l++){
				r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*(vxy+vyx)));
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vyy));
				q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vxx));
				sumr += r[j][i][l];
				sump += p[j][i][l];
				sumq += q[j][i][l];
			}

			/* and now the components of the stress tensor are
			   completely updated */
			sxy[j][i]+=(dthalbe*sumr);
			sxx[j][i]+=(dthalbe*sump);
			syy[j][i]+=(dthalbe*sumq);
			
			/* save forward wavefield for gradient calculation */
			if(INVMAT==0){
			  ux[j][i] = vxx;
		          uy[j][i] = vyy;
			  uxy[j][i] = vxy + vyx;
			}

   }}
		break;

/*	case 10:
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
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);*/
				
				/*
				piv = pi[j][i]*(vxx+vyy);
				sxy[j][i]+=(uipjp[j][i]*(vxy+vyx));
				sxx[j][i]+=piv-(2.0*u[j][i]*vyy);
				syy[j][i]+=piv-(2.0*u[j][i]*vxx);
				*/
/*			}
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
*/
		
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
				vxx *= dhi;
				vyy *= dhi;
				vyx *= dhi;
				vxy *= dhi;

				sumr=sump=sumq=0.0;
				for (l=1;l<=L;l++){
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += (fipjp[j][i]*(vxy+vyx))+(dthalbe*sumr);
				sxx[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vyy)+(dthalbe*sump);
				syy[j][i] += (g[j][i]*(vxx+vyy))-(2.0*f[j][i]*vxx)+(dthalbe*sumq);

				sumr=sump=sumq=0.0;
				for (l=1;l<=L;l++){
					r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*(vxy+vyx)));
					p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vyy));
					q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy))+(2.0*d[j][i][l]*vxx));
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+=(dthalbe*sumr);
				sxx[j][i]+=(dthalbe*sump);
				syy[j][i]+=(dthalbe*sumq);
			}
		}
		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
