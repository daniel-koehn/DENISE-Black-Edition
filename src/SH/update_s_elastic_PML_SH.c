/*  Update_s_elastic_PML_SH
 *
 *  updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *  by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *  and second order accuracy in time for the elastic SH problem
 *   
 *  Daniel Koehn
 *  Kiel, 09.12.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_PML_SH(int nx1, int nx2, int ny1, int ny2,
	float ** vz, float **  uz, float **  uzx, float **   syz, float **   sxz,
	float ** ujp, float ** uip, float **rho, float *hc, int infoout,
        float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_vzx, float ** psi_vzy, int mode){

	int i,j, m, h, h1, l;
	float  vzx, vzy;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, FW, L, GRAD_FORM;
        extern int FREE_SURF, BOUNDARY, ADJ_SIGN;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;	
	
	dhi=DT/DH;
	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	
	switch (FDORDER){

	case 2:
		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]))*dhi;
                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i]))*dhi;


        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

	  		/* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
						
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}
			
		    }
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]) + hc[2]*(vz[j][i+2]-vz[j][i-1]))*dhi;
                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i]) + hc[2]*(vz[j+2][i]-vz[j-1][i]))*dhi;


        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

		        /* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
	
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}

		    }
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]) + hc[2]*(vz[j][i+2]-vz[j][i-1]) + hc[3]*(vz[j][i+3]-vz[j][i-2]))*dhi;
                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i]) + hc[2]*(vz[j+2][i]-vz[j-1][i]) + hc[3]*(vz[j+3][i]-vz[j-2][i]))*dhi;


        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

			/* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
	
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}

		    }
		}
		break;

	case 8:
		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]) 
			     + hc[2]*(vz[j][i+2]-vz[j][i-1])
			     + hc[3]*(vz[j][i+3]-vz[j][i-2])
			     + hc[4]*(vz[j][i+4]-vz[j][i-3]))*dhi;

                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i])
		             + hc[2]*(vz[j+2][i]-vz[j-1][i])
		             + hc[3]*(vz[j+3][i]-vz[j-2][i])
			     + hc[4]*(vz[j+4][i]-vz[j-3][i]))*dhi; 

        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

			/* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
	
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}

		    }
		}
		break;

	case 10:

		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]) 
			     + hc[2]*(vz[j][i+2]-vz[j][i-1])
			     + hc[3]*(vz[j][i+3]-vz[j][i-2])
			     + hc[4]*(vz[j][i+4]-vz[j][i-3])
			     + hc[5]*(vz[j][i+5]-vz[j][i-4]))*dhi;

                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i])
		             + hc[2]*(vz[j+2][i]-vz[j-1][i])
		             + hc[3]*(vz[j+3][i]-vz[j-2][i])
			     + hc[4]*(vz[j+4][i]-vz[j-3][i])
			     + hc[5]*(vz[j+5][i]-vz[j-4][i]))*dhi; 

        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

			/* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
	
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}

		    }
		}
		break;
		
	case 12:

		for (j=ny1;j<=ny2;j++){
		    for (i=nx1;i<=nx2;i++){

                        vzx = (hc[1]*(vz[j][i+1]-vz[j][i]) 
			     + hc[2]*(vz[j][i+2]-vz[j][i-1])
			     + hc[3]*(vz[j][i+3]-vz[j][i-2])
			     + hc[4]*(vz[j][i+4]-vz[j][i-3])
			     + hc[5]*(vz[j][i+5]-vz[j][i-4])
			     + hc[6]*(vz[j][i+6]-vz[j][i-5]))*dhi;

                        vzy = (hc[1]*(vz[j+1][i]-vz[j][i])
		             + hc[2]*(vz[j+2][i]-vz[j-1][i])
		             + hc[3]*(vz[j+3][i]-vz[j-2][i])
			     + hc[4]*(vz[j+4][i]-vz[j-3][i])
			     + hc[5]*(vz[j+5][i]-vz[j-4][i])
			     + hc[6]*(vz[j+6][i]-vz[j-5][i]))*dhi; 

        		/* left boundary */                                         
        		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vzx[j][i] = b_x_half[i] * psi_vzx[j][i] + a_x_half[i] * vzx;
                        	vzx = vzx / K_x_half[i] + psi_vzx[j][i];                 
         		}

        		/* right boundary */                                         
        		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        	h1 = (i-nx2+2*FW);
                        	h = i;                        
                        
                        	psi_vzx[j][h1] = b_x_half[h1] * psi_vzx[j][h1] + a_x_half[h1] * vzx;
				vzx = vzx / K_x_half[h1] + psi_vzx[j][h1];
                                           
         		}

			/* top boundary */                                         
        		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        	psi_vzy[j][i] = b_y_half[j] * psi_vzy[j][i] + a_y_half[j] * vzy;                                                                 
                        	vzy = vzy / K_y_half[j] + psi_vzy[j][i];
                        	
        		}
	
	  		/* bottom boundary */                                         
        		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        	h1 = (j-ny2+2*FW);                                        
                        	h = j;                                                
                        
                        	psi_vzy[h1][i] = b_y_half[h1] * psi_vzy[h1][i] + a_y_half[h1] * vzy;
				vzy = vzy / K_y_half[h1] + psi_vzy[h1][i];
        
        		}
	
			/* update stress components */
			sxz[j][i] += ADJ_SIGN*uip[j][i]*vzx;
			syz[j][i] += ADJ_SIGN*ujp[j][i]*vzy;

			/* save forward wavefield for gradient calculation */
			if(mode==0){
			  uz[j][i] = vzx/DT;
			  uzx[j][i] = vzy/DT;
			}

		    }
		}
		break;

		
	default:

		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
