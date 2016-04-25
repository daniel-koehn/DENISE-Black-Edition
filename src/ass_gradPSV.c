/*
 * Assemble gradients for each shot (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 24/04/2016
 */

#include "fd.h"

void ass_gradPSV(struct fwiPSV *fwiPSV, struct matPSV *matPSV, int iter){

        /* global variables */
	extern int NX, NY, IDX, IDY, INVMAT1;
        extern int GRAD_FORM;
        extern int INV_VP_ITER, INV_VS_ITER, INV_RHO_ITER;
	extern float DT;

	/* local variables */
	int i, j;
        float muss, lamss;	

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
	 
		     		muss = (*matPSV).pu[j][i];
        	    		lamss = (*matPSV).ppi[j][i];
	    
		    		if((muss>0.0)||(lamss>0.0)){
        	       		  (*fwiPSV).waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * (*fwiPSV).waveconv_lam[j][i];
		    		}
	    
          		}

          		(*fwiPSV).waveconv_shot[j][i] = (*fwiPSV).waveconv_lam[j][i];
		}
			 
        	if(INVMAT1==1){

	  		/* calculate Vp gradient */ 
	  		if(GRAD_FORM==1){
	   
              			muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
              			lamss = (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] - 2.0 *  muss;
	      
	      			if((muss>0.0)||(lamss>0.0)){
	          			(*fwiPSV).waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * (*fwiPSV).waveconv_lam[j][i];
	      			}
	      
	      			(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matPSV).ppi[j][i] * (*matPSV).prho[j][i] * (*fwiPSV).waveconv_lam[j][i];

          		}

	          	if(GRAD_FORM==2){                             	
	   
	      			muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
	      			lamss = (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] - 2.0 *  muss;
	      
	      			if((muss>0.0)||(lamss>0.0)){
		  			(*fwiPSV).waveconv_lam[j][i] = (1.0/(4.0*(lamss+muss) * (lamss+muss))) * (*fwiPSV).waveconv_lam[j][i];
	      			}
	      
	      			(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matPSV).ppi[j][i] * (*matPSV).prho[j][i] * (*fwiPSV).waveconv_lam[j][i]; 
	      
          		}
		    		   
		 		
		}
		 
        	if(INVMAT1==2){	
	   		/* calculate Zp gradient */
           		(*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matPSV).ppi[j][i] * (*fwiPSV).waveconv_lam[j][i];	   
		}
	
        	if(iter<INV_VP_ITER){
           		(*fwiPSV).waveconv_shot[j][i] = 0.0;
        	}
	                                                                       
      		}
   	}

	/* calculate gradient for mu, Vs or Zs */
	/* ----------------------------------- */

	for (i=1;i<=NX;i=i+IDX){
   		for (j=1;j<=NY;j=j+IDY){
		 
      		/* calculate mu gradient */ 
      		(*fwiPSV).waveconv_mu[j][i] = -DT * (*fwiPSV).waveconv_u_shot[j][i];
		 
      		if(INVMAT1==1){		
         		/* calculate Vs gradient */		 
         		(*fwiPSV).waveconv_u_shot[j][i] = (- 4.0 * (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*fwiPSV).waveconv_lam[j][i]) + 2.0 * (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*fwiPSV).waveconv_mu[j][i];         	 
      		}
		 
      		if(INVMAT1==2){
        		/* calculate Zs gradient */
        		(*fwiPSV).waveconv_u_shot[j][i] = (- 4.0 * (*matPSV).pu[j][i] * (*fwiPSV).waveconv_lam[j][i]) + (2.0 * (*matPSV).pu[j][i] * (*fwiPSV).waveconv_mu[j][i]);
      		}
		 
      		if(INVMAT1==3){
        		/* calculate u gradient */
        		(*fwiPSV).waveconv_u_shot[j][i] = (*fwiPSV).waveconv_mu[j][i];
      		}

      		if(iter<INV_VS_ITER){
         		(*fwiPSV).waveconv_u_shot[j][i] = 0.0;
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
          		(*fwiPSV).waveconv_rho_shot[j][i] = (((((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i])-(2.0 * (*matPSV).pu[j][i] * (*matPSV).pu[j][i])) * (*fwiPSV).waveconv_lam[j][i]) 
                         		                     + ((*matPSV).pu[j][i] * (*matPSV).pu[j][i] * (*fwiPSV).waveconv_mu[j][i]) + (*fwiPSV).waveconv_rho_s[j][i]);
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



