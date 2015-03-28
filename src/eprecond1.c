#include "fd.h"
void eprecond1(float ** We, float ** Ws, float ** Wr){
  
	extern int NX, NY, IDX, IDY, DTINV, EPRECOND;
        extern int POS[3], NXG;
        extern float DH;
	int i, j, k, l, ii, jj;
	float maxWetmp, maxWe, EPSILON_WE, x, y, xmin, xmax;
        	
        xmin = DH;
        xmax = NXG*DH;

	EPSILON_WE = 0.005;
	/* EPSILON_WE = 0.05; */
	
	maxWetmp=0.0;
	
	if(EPRECOND==1){
	  /* calculate energy weighting */
	  for (i=1;i<=NX;i=i+IDX){
	 	  for (j=1;j<=NY;j=j+IDY){
			     
		      We[j][i]=sqrt(Ws[j][i]*Wr[j][i]); /* energy weighted source and receiver contribution */
		    
		      /* estimate maximum We on this CPU*/
		      if(We[j][i]>maxWetmp){
		        maxWetmp = We[j][i]; 
		      }
		        
		  }
	  }
	}
	
	if(EPRECOND==3){
	/* Forward wavefield + approximation of the receiver Greens function (Plessix and Mulder, 2004) */
	  for (i=1;i<=NX;i=i+IDX){
	 	  for (j=1;j<=NY;j=j+IDY){

                      /* calculate global coordinates */		      
                      ii=i+POS[1]*NX;
                      jj=j+POS[2]*NY;

                      x = ii*DH;
                      y = jj*DH;
  
		      We[j][i]= sqrt(Ws[j][i]) * (asinh((xmax-x)/y)-asinh((xmin-x)/y));
		    
		      /* estimate maximum We on this CPU*/
		      if(We[j][i]>maxWetmp){
		        maxWetmp = We[j][i]; 
		      }
		        
		  }
	  }
	}
                                                                                                                                                                                                             
	/* estimate maximum of We */
	MPI_Allreduce(&maxWetmp,&maxWe,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	/* regularize energy weighting to avoid divison by zero */
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
	                        
			We[j][i] = We[j][i] + (EPSILON_WE*maxWe);
	                                                
		}
	}
		
						
}
