/*  Update_v_PML_AC
 *
 *  updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *  by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *  and second order accuracy in time for the acoustic problem
 *   
 *  Daniel Koehn
 *  Kiel, 10.06.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void update_v_PML_AC(int nx1, int nx2, int ny1, int ny2, int nt,
	float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float **  uttx, float **  utty,float ** p,
	float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
	float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_p_x, float ** psi_p_y){

	int i, j,l,fdoh,m, h, h1;
	float amp, dtdh, azi_rad;
	float vxtmp, vytmp;
        float p_x, p_y;
	
	extern float DT, DH, ANGLE;
	double time1, time2;
	extern int MYID, QUELLTYP, QUELLTYPB, FDORDER;
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

                           p_x =  hc[1]*(p[j][i+1]-p[j][i]);
                           p_y = hc[1]*(p[j+1][i]-p[j][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];

              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                         
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];                           

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                  
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }
                           
                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 		         
      
     }
}
		break;
		
	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

                           p_x =  hc[1]*(p[j][i+1]-p[j][i])
					    + hc[2]*(p[j][i+2]-p[j][i-1]);  

                           p_y = hc[1]*(p[j+1][i]-p[j][i])
					    + hc[2]*(p[j+2][i]-p[j-1][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                         
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];
                            

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                  
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }
                           
                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 		         
      
     }
}
	     
		break;
		
	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

                           p_x =  hc[1]*(p[j][i+1]-p[j][i])
					    + hc[2]*(p[j][i+2]-p[j][i-1])
					    + hc[3]*(p[j][i+3]-p[j][i-2]);

                           p_y = hc[1]*(p[j+1][i]-p[j][i])
					    + hc[2]*(p[j+2][i]-p[j-1][i])
					    + hc[3]*(p[j+3][i]-p[j-2][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                         
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];                           

        }                       
                           
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                  
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}   

                              if(sw==1){
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }

                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }
                           
                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 		         
      
     }
}

	     
		break;
		
	case 8:

for (j=ny1;j<=ny2;j++){
	for (i=nx1;i<=nx2;i++){

                           p_x =  hc[1]*(p[j][i+1]-p[j][i])
					    + hc[2]*(p[j][i+2]-p[j][i-1])
					    + hc[3]*(p[j][i+3]-p[j][i-2])
					    + hc[4]*(p[j][i+4]-p[j][i-3]);

                           p_y = hc[1]*(p[j+1][i]-p[j][i])
					    + hc[2]*(p[j+2][i]-p[j-1][i])
					    + hc[3]*(p[j+3][i]-p[j-2][i])
					    + hc[4]*(p[j+4][i]-p[j-3][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                        
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];
                           

        }                       
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                 
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}                 

                              if(sw==1){ 
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }
                           
                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }   

                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 		         
      
     }
}

	      
		break;
		
	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				

			   p_x =  hc[1]*(p[j][i+1]-p[j][i])
					    + hc[2]*(p[j][i+2]-p[j][i-1])
					    + hc[3]*(p[j][i+3]-p[j][i-2])
					    + hc[4]*(p[j][i+4]-p[j][i-3])
					    + hc[5]*(p[j][i+5]-p[j][i-4]);

                           p_y = hc[1]*(p[j+1][i]-p[j][i])
					    + hc[2]*(p[j+2][i]-p[j-1][i])
					    + hc[3]*(p[j+3][i]-p[j-2][i])
					    + hc[4]*(p[j+4][i]-p[j-3][i])
					    + hc[5]*(p[j+5][i]-p[j-4][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];

              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                         
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];                           

        }                       

                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                 
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}                 

                              if(sw==1){ 
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }
                           
                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }

                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 
				
			}
		}
		break;

	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
			   p_x =  hc[1]*(p[j][i+1]-p[j][i])
					    + hc[2]*(p[j][i+2]-p[j][i-1])
					    + hc[3]*(p[j][i+3]-p[j][i-2])
					    + hc[4]*(p[j][i+4]-p[j][i-3])
					    + hc[5]*(p[j][i+5]-p[j][i-4])
					    + hc[6]*(p[j][i+6]-p[j][i-5]);

                           p_y = hc[1]*(p[j+1][i]-p[j][i])
					    + hc[2]*(p[j+2][i]-p[j-1][i])
					    + hc[3]*(p[j+3][i]-p[j-2][i])
					    + hc[4]*(p[j+4][i]-p[j-3][i])
					    + hc[5]*(p[j+5][i]-p[j-4][i])
					    + hc[6]*(p[j+6][i]-p[j-5][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                           psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           /*psi_p_x[j][i] = b_x_half[i] * psi_p_x[j][i] + a_x_half[i] * p_x;                                                
                           p_x = p_x / K_x_half[i] + psi_p_x[j][i];*/
                           
                           psi_p_x[j][h1] = b_x_half[h1] * psi_p_x[j][h1] + a_x_half[h1] * p_x;                                                
			   p_x = p_x / K_x_half[h1] + psi_p_x[j][h1];

         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		           psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i]; 
                         
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       		           /*psi_p_y[j][i] = b_y_half[j] * psi_p_y[j][i] + a_y_half[j] * p_y;                                                
                           p_y = p_y / K_y_half[j] + psi_p_y[j][i];*/
                           
                           psi_p_y[h1][i] = b_y_half[h1] * psi_p_y[h1][i] + a_y_half[h1] * p_y;                                                
			   p_y = p_y / K_y_half[h1] + psi_p_y[h1][i];                           

        }                       
                           if(GRAD_FORM==1){

                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;                 
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;}                 

                              if(sw==1){ 
                                 vxp1[j][i] += DT*vx[j][i];
                                 vyp1[j][i] += DT*vy[j][i];}

                           }
                           
                           if(GRAD_FORM==2){

                              if(sw==1){
				 vxp1[j][i] = vx[j][i];
			         vyp1[j][i] = vy[j][i];			      
 			      }
                           
                              if(sw==0){
                                 vxp1[j][i] = rip[j][i]*p_x/DH;
                                 vyp1[j][i] = rjp[j][i]*p_y/DH;			      			         

                              }
			   }

                           vx[j][i] += DT*rip[j][i]*p_x/DH;
		           vy[j][i] += DT*rjp[j][i]*p_y/DH; 
			
			
			}
		}
		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
					vxtmp +=   hc[m]*( p[j][i+m]   - p[j][i-m+1] );							
					vytmp +=   hc[m]*( p[j+m][i]   - p[j-m+1][i] );
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
		    if(QUELLTYP==4){vx[j][i] +=  -sin(azi_rad) * signals1[l][nt];    /* rotated force in x */
		                    vy[j][i] +=   cos(azi_rad) * signals1[l][nt];}   /* rotated force in y */          
		              
		}}
		
		/* Backpropagation (sw==1) */
		if(sw==1){
		for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    
                    if((GRAD_FORM==1)||(GRAD_FORM==2)){
		       if(QUELLTYPB==1){vx[j][i] += signals[l][nt];    /* single force in x */
		                        vy[j][i] += signals1[l][nt];}  /* + single force in y */

		       if(QUELLTYPB==2||QUELLTYPB==6||QUELLTYPB==7){vy[j][i] += signals1[l][nt];}  /* single force in y */
		       if(QUELLTYPB==3||QUELLTYPB==5||QUELLTYPB==7){vx[j][i] += signals[l][nt];}   /* single force in x */
                    }

		}}                         
				

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
