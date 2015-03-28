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
	dtdh = DT*DT/DH;
        
     /* drad = PI/180.0; 
        angle = 135.0; */
         
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}


	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */

	switch (FDORDER){
	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
                           sxx_x =  hc[1]*(sxx[j][i+1]-sxx[j][i]);

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1]);

		           sxy_y =  hc[1]*(sxy[j][i]-sxy[j-1][i]);

                           syy_y = hc[1]*(syy[j+1][i]-syy[j][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
                           
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
		            			   
       		           /*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
                           syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
                           
                           psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * syy_y;                                                
			   syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i];
                           
                           psi_sxy_y[h1][i] = b_y[h1] * psi_sxy_y[h1][i] + a_y[h1] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i]; 

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;                  
                                 vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

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
		break;
		
	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
                           sxx_x =  hc[1]*(sxx[j][i+1]-sxx[j][i])
					    + hc[2]*(sxx[j][i+2]-sxx[j][i-1]);

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2]);

		           sxy_y =  hc[1]*(sxy[j][i]-sxy[j-1][i])
					    + hc[2]*(sxy[j+1][i]-sxy[j-2][i]);

                           syy_y = hc[1]*(syy[j+1][i]-syy[j][i])
					    + hc[2]*(syy[j+2][i]-syy[j-1][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
                           
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
		            			   
       		           /*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
                           syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
                           
                           psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * syy_y;                                                
			   syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i];
                           
                           psi_sxy_y[h1][i] = b_y[h1] * psi_sxy_y[h1][i] + a_y[h1] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i]; 

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;                  
                                 vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

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
	     
		break;
		
	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
                           sxx_x =  hc[1]*(sxx[j][i+1]-sxx[j][i])
					    + hc[2]*(sxx[j][i+2]-sxx[j][i-1])
					    + hc[3]*(sxx[j][i+3]-sxx[j][i-2]);

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3]);

		           sxy_y =  hc[1]*(sxy[j][i]-sxy[j-1][i])
					    + hc[2]*(sxy[j+1][i]-sxy[j-2][i])
					    + hc[3]*(sxy[j+2][i]-sxy[j-3][i]);

                           syy_y = hc[1]*(syy[j+1][i]-syy[j][i])
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
                           
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
		            			   
       		           /*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
                           syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
                           
                           psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * syy_y;                                                
			   syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i];
                           
                           psi_sxy_y[h1][i] = b_y[h1] * psi_sxy_y[h1][i] + a_y[h1] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i]; 

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;                  
                                 vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

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

	     
		break;
		
	case 8:

for (j=ny1;j<=ny2;j++){
	for (i=nx1;i<=nx2;i++){

                           sxx_x =  hc[1]*(sxx[j][i+1]-sxx[j][i])
					    + hc[2]*(sxx[j][i+2]-sxx[j][i-1])
					    + hc[3]*(sxx[j][i+3]-sxx[j][i-2])
					    + hc[4]*(sxx[j][i+4]-sxx[j][i-3]);

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4]);

		           sxy_y =  hc[1]*(sxy[j][i]-sxy[j-1][i])
					    + hc[2]*(sxy[j+1][i]-sxy[j-2][i])
					    + hc[3]*(sxy[j+2][i]-sxy[j-3][i])
					    + hc[4]*(sxy[j+3][i]-sxy[j-4][i]);

                           syy_y = hc[1]*(syy[j+1][i]-syy[j][i])
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
                           sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
                           
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
		            			   
       		           /*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
                           syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
                           
                           psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * syy_y;                                                
			   syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i];
                           
                           psi_sxy_y[h1][i] = b_y[h1] * psi_sxy_y[h1][i] + a_y[h1] * sxy_y;                                                
                           sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i]; 

        }                       
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*(sxx_x+sxy_y)/DH;                 
                                 vyp1[j][i] = rjp[j][i]*(sxy_x+syy_y)/DH;}                 

                              if(sw==1){ 
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }
                           
                           if(GRAD_FORM==2){

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

	      
		break;
		
	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
					    + hc[2]*(sxx[j][i+2]-sxx[j][i-1])
					    + hc[3]*(sxx[j][i+3]-sxx[j][i-2])
					    + hc[4]*(sxx[j][i+4]-sxx[j][i-3])
					    + hc[5]*(sxx[j][i+5]-sxx[j][i-4])
					    + hc[1]*(sxy[j][i]  -sxy[j-1][i])
					    + hc[2]*(sxy[j+1][i]-sxy[j-2][i])
					    + hc[3]*(sxy[j+2][i]-sxy[j-3][i])
					    + hc[4]*(sxy[j+3][i]-sxy[j-4][i])
					    + hc[5]*(sxy[j+4][i]-sxy[j-5][i])
					   )*dtdh*rip[j][i];
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i])
					    + hc[5]*(syy[j+5][i]-syy[j-4][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4])
					    + hc[5]*(sxy[j][i+4]-sxy[j][i-5])
					   )*dtdh*rjp[j][i];
			}
		}
		break;

	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
					    + hc[2]*(sxx[j][i+2]-sxx[j][i-1])
					    + hc[3]*(sxx[j][i+3]-sxx[j][i-2])
					    + hc[4]*(sxx[j][i+4]-sxx[j][i-3])
					    + hc[5]*(sxx[j][i+5]-sxx[j][i-4])
					    + hc[6]*(sxx[j][i+6]-sxx[j][i-5])
					    + hc[1]*(sxy[j][i]  -sxy[j-1][i])
					    + hc[2]*(sxy[j+1][i]-sxy[j-2][i])
					    + hc[3]*(sxy[j+2][i]-sxy[j-3][i])
					    + hc[4]*(sxy[j+3][i]-sxy[j-4][i])
					    + hc[5]*(sxy[j+4][i]-sxy[j-5][i])
					    + hc[6]*(sxy[j+5][i]-sxy[j-6][i])
					   )*dtdh*rip[j][i];
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i])
					    + hc[5]*(syy[j+5][i]-syy[j-4][i])
					    + hc[6]*(syy[j+6][i]-syy[j-5][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4])
					    + hc[5]*(sxy[j][i+4]-sxy[j][i-5])
					    + hc[6]*(sxy[j][i+5]-sxy[j][i-6])
					   )*dtdh*rjp[j][i];
			}
		}
		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
					vxtmp +=   hc[m]*( sxx[j][i+m]   - sxx[j][i-m+1] )
						 + hc[m]*( sxy[j+m-1][i] - sxy[j-m][i]   );
							
					vytmp +=   hc[m]*( syy[j+m][i]   - syy[j-m+1][i] )
						 + hc[m]*( sxy[j][i+m-1] - sxy[j][i-m]   );
				}
					
				vx[j][i] += vxtmp*dtdh*rip[j][i];
				vy[j][i] += vytmp*dtdh*rjp[j][i];
			}
		}
		break;

	} /* end of switch(FDORDER) */


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
		              
		}}
		
		/* Backpropagation (sw==1) */
		if(sw==1){
		for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    
		    if(QUELLTYPB==1){vx[j][i] += signals[l][nt];    /* single force in x */
		                     vy[j][i] += signals1[l][nt];}  /* + single force in y */

		    if(QUELLTYPB==2){vy[j][i] += signals1[l][nt];}  /* single force in y */
		    if(QUELLTYPB==3){vx[j][i] += signals[l][nt];}   /* single force in x */

		}}                         
			
	
	/*if (FW>0.0)
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			vx[j][i]*=absorb_coeff[j][i];
			vy[j][i]*=absorb_coeff[j][i];
			sxy[j][i]*=absorb_coeff[j][i];
			sxx[j][i]*=absorb_coeff[j][i];
			syy[j][i]*=absorb_coeff[j][i];
		
	}}*/

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
