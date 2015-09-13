/* $Id: update_v_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen 
 *
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void update_v_PML(int nx1, int nx2, int ny1, int ny2, int nt,
	float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float **  uttx, float **  utty,float ** sxx, float ** syy,
	float ** sxy, float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
	float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_sxx_x, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_sxy_x){

	int i, j,l,fdoh,m, h, h1;
	float amp, dtdh, azi_rad;
	float vxtmp, vytmp;
        float sxx_x, syy_y, sxy_y, sxy_x;
	
	extern float DT, DH, ANGLE;
	double time1, time2;
	extern int MYID, QUELLTYP, QUELLTYPB, CHECKPTREAD, FDORDER;
        extern int FDORDER, INVMAT1, GRAD_FORM;
        extern int FREE_SURF, BOUNDARY, FW;
        extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;

	fdoh = FDORDER/2;
         
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}

	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */

        for (j=ny1;j<=ny2;j++){
	   for (i=nx1;i<=nx2;i++){

              sxx_x = 0.0;
	      sxy_x = 0.0;
	      sxy_y = 0.0;
	      syy_y = 0.0;
	      
	      /* calculate spatial derivatives*/
	      for (m=1; m<=fdoh; m++) {
	      
		 sxx_x +=  hc[m]*( sxx[j][i+m]  - sxx[j][i-m+1] );
		 sxy_x +=  hc[m]*( sxy[j][i+m-1]  - sxy[j][i-m] );
		 sxy_y +=  hc[m]*( sxy[j+m-1][i]  - sxy[j-m][i] );
		 syy_y +=  hc[m]*( syy[j+m][i]  - syy[j-m+1][i] );
				
	      }
				                                                     
	      /* left boundary */                                         
              if ((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
              }

              /* right boundary */                                         
              if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                           h = i;                                 
                           
                           psi_sxx_x[j][h1] = b_x_half[h1] * psi_sxx_x[j][h1] + a_x_half[h1] * sxx_x;                                                
			   sxx_x = sxx_x / K_x_half[h1] + psi_sxx_x[j][h1];

                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
              }

         
	      /* top boundary */                                         
              if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
                           syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i]; 
                        
                           psi_sxy_y[j][i] = b_y[j] * psi_sxy_y[j][i] + a_y[j] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[j] + psi_sxy_y[j][i]; 
              }
		

	      /* bottom boundary */                                         
              if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		           h = j;		            			          		       
                           
                           psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * syy_y;                                                
			   syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i];
                           
                           psi_sxy_y[h1][i] = b_y[h1] * psi_sxy_y[h1][i] + a_y[h1] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i]; 

              }
	                             
              if (GRAD_FORM==1){

                           if(sw==0){
                              vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;                 
                              vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}                 

                           if(sw==1){ 
                              vxp1[j][i] += DT*vx[j][i];
                              vyp1[j][i] += DT*vy[j][i];}

              }
                           
              if((GRAD_FORM==2)||(GRAD_FORM==3)||(GRAD_FORM==4)){

                           if(sw==0){
                              vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;
                              vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}
                           
                           if(sw==1){
			      vxp1[j][i] = vx[j][i];
			      vyp1[j][i] = vy[j][i];}

              }

              vx[j][i] += DT*rip[j][i]*(sxx_x+sxy_y)/DH;
	      vy[j][i] += DT*rjp[j][i]*(sxy_x+syy_y)/DH; 		         
      
            }
        }

 	/* Forward Modelling (sw==0) */
	if(sw==0){
	
	        for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    azi_rad=srcpos_loc[7][l]*PI/180;
		    QUELLTYP=(int)srcpos_loc[8][l];
		
		    if(QUELLTYP==2){vx[j][i] +=  signals1[l][nt];}  /* single force in x */
		    if(QUELLTYP==3){vy[j][i] +=  signals1[l][nt];}  /* single force in y */
		    if(QUELLTYP==4){vx[j][i] +=  sin(azi_rad) * signals1[l][nt];    /* rotated force in x */
		                    vy[j][i] +=  cos(azi_rad) * signals1[l][nt];}  /* rotated force in y */          
		              
		}
	 }
		
	/* Backpropagation (sw==1) */
	if(sw==1){
	
		for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    
                    if((GRAD_FORM==1)||(GRAD_FORM==3)||(GRAD_FORM==4)){
		       if(QUELLTYPB==1){vx[j][i] += signals[l][nt];    /* single force in x */
		                        vy[j][i] += signals1[l][nt];}  /* + single force in y */

		       if(QUELLTYPB==2){vy[j][i] += signals1[l][nt];}  /* single force in y */
		       if(QUELLTYPB==3){vx[j][i] += signals[l][nt];}   /* single force in x */
                    }

                    if(GRAD_FORM==2){
                       if(QUELLTYPB==1){vx[j][i] += rip[j][i]*signals[l][nt];    /* single force in x */
                                        vy[j][i] += rjp[j][i]*signals1[l][nt];}  /* + single force in y */

                       if(QUELLTYPB==2){vy[j][i] += rjp[j][i]*signals1[l][nt];}  /* single force in y */
                       if(QUELLTYPB==3){vx[j][i] += rip[j][i]*signals[l][nt];}   /* single force in x */
                    }


		}
	}                         			

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
