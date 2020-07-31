/*
 * Assemble gradients for each shot (AC problem) 
 *
 * Daniel Koehn
 * Kiel, 12/06/2017
 */

#include "fd.h"

void ass_gradAC(struct fwiPSV *fwiPSV, struct matAC *matAC, int iter){

        /* global variables */
	extern int NX, NY, IDX, IDY, INVMAT1;
        extern int GRAD_FORM;
        extern int INV_VP_ITER, INV_RHO_ITER;
	extern float DT;

	/* local variables */
	int i, j;
        float lamss;	

	/* calculate gradient for lambda, Vp or Zp */
	/* --------------------------------------- */	

  	/* norm(fwiPSV.waveconv_shot);
  	norm(fwiPSV.waveconv_u_shot);  
  	norm(fwiPSV.waveconv_rho_shot); */

  	for (i=1;i<=NX;i=i+IDX){
     		for (j=1;j<=NY;j=j+IDY){

       		/* calculate lambda gradient */           
       		(*fwiPSV).waveconv_lam[j][i] = - DT * (*fwiPSV).waveconv_shot[j][i];
	 
       		if(INVMAT1==3){

         		if(GRAD_FORM==1){
	 
        	    		lamss = (*matAC).ppi[j][i];
	    
		    		if(lamss>0.0){
        	       		  (*fwiPSV).waveconv_lam[j][i] = (1.0/(4.0 * lamss * lamss)) * (*fwiPSV).waveconv_lam[j][i];
		    		}
	    
          		}
			
			if(GRAD_FORM==2){
	 
        	    		lamss = (*matAC).ppi[j][i];
	    
		    		if(lamss>0.0){
        	       		  (*fwiPSV).waveconv_lam[j][i] = - (1.0/(lamss * lamss)) * (*fwiPSV).waveconv_lam[j][i];
		    		}
	    
          		}
			

          		(*fwiPSV).waveconv_shot[j][i] = (*fwiPSV).waveconv_lam[j][i];
		}
			 
        	if(INVMAT1==1){

	  		/* calculate Vp gradient */ 
	  		if(GRAD_FORM==1){
	   
              			lamss = (*matAC).prho[j][i] * (*matAC).ppi[j][i] * (*matAC).ppi[j][i];
	      
	      			if(lamss>0.0){
	          			(*fwiPSV).waveconv_lam[j][i] = (1.0/(4.0 * lamss * lamss)) * (*fwiPSV).waveconv_lam[j][i];
	      			}
	      
	      			(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matAC).ppi[j][i] * (*matAC).prho[j][i] * (*fwiPSV).waveconv_lam[j][i];

          		}

	          	if(GRAD_FORM==2){                             	
	   
	      			lamss = (*matAC).prho[j][i] * (*matAC).ppi[j][i] * (*matAC).ppi[j][i];
	      
	      			if(lamss>0.0){
		  			(*fwiPSV).waveconv_lam[j][i] = -(1.0/(lamss * lamss)) * (*fwiPSV).waveconv_lam[j][i];
	      			}
	      
	      			(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matAC).ppi[j][i] * (*matAC).prho[j][i] * (*fwiPSV).waveconv_lam[j][i]; 
	      
          		}
		    		   
		 		
		}
		 
        	if(INVMAT1==2){	
	   		/* calculate Zp gradient */
           		(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matAC).ppi[j][i] * (*fwiPSV).waveconv_lam[j][i];	   
		}
	
        	if(iter<INV_VP_ITER){
           		(*fwiPSV).waveconv_shot[j][i] = 0.0;
        	}
	                                                                       
      		}
   	}


	/* calculate gradient for density */
	/* ------------------------------ */
	for (i=1;i<=NX;i=i+IDX){
    		for (j=1;j<=NY;j=j+IDY){

       		/* calculate density gradient rho' */
       		(*fwiPSV).waveconv_rho_s[j][i]= - DT * (*fwiPSV).waveconv_rho_shot[j][i];
				 
       		if(INVMAT1==1){
          		/* calculate density gradient */
          		(*fwiPSV).waveconv_rho_shot[j][i] = ((((*matAC).ppi[j][i] * (*matAC).ppi[j][i]) * (*fwiPSV).waveconv_lam[j][i]) 
                         		                     + (*fwiPSV).waveconv_rho_s[j][i]);
       		}
		 
       		if(INVMAT1==3){
          		/* calculate density gradient */
          		(*fwiPSV).waveconv_rho_shot[j][i] = (*fwiPSV).waveconv_rho_s[j][i];
       		}

       		if(iter<INV_RHO_ITER){
          		(*fwiPSV).waveconv_rho_shot[j][i] = 0.0;
       		}
	 
    		}
	}
	
}



