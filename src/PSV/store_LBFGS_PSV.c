/*------------------------------------------------------------------------
 * Store gradients and model parameters for 
 * Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS)
 * optimization (Nocedal and Wright, 2006) in case of the PSV problem
 * 
 * Daniel Koehn
 * Kiel, 25.07.2016
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"

void store_LBFGS_PSV(float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float ** waveconv, float ** gradp, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float **ppi, float ** pu, float ** prho, int nxnyi, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec){

	extern int NX, NY, IDX, IDY, SPATFILTER;
	extern int HESSIAN, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	extern int POS[3], MYID;
        extern float C_vp, C_vs, C_rho;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, k, h, h1, h2;
	float betaz, betan, gradplastiter, gradclastiter, betar, beta;
	float gamma_LBFGS, sum_nom, sum_denom;
        float LBFGSTMP, LBFGSTMP1, LBFGSTMP2, LBFGSTMP3, modellastiter, norm_fac, norm_fac_u, norm_fac_rho;
        float beta_LBFGS_1;
        int ki, itershift, iter1;
	extern FILE *FP;
	FILE *FP3, *FP4, *FP6, *FP5, *FP7;
	
        itershift = 1;

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Vp/Zp/lambda ========================================================================= */
/* ===================================================================================================================================================== */
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv[j][i] = C_vp * waveconv[j][i];
   }
}

/* apply median filter at source positions */
/*median_src(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv);
}

if(SPATFILTER==2){
  smooth2(waveconv);
}
  
/* Normalize gradient to maximum value */
/*norm_fac_u=norm(waveconv_u,iter,2);
if(MYID==0){printf("norm_fac_u=%e \n",norm_fac_u);}*/
  
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp[j][i] = waveconv[j][i];
   }
}

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Vs/Zs/mu ============================================================================= */
/* ===================================================================================================================================================== */
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_u[j][i] = C_vs * waveconv_u[j][i];
   }
}

/* apply median filter at source positions */
/*median_src(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv_u);
}

if(SPATFILTER==2){
  smooth2(waveconv_u);
}
  
/* Normalize gradient to maximum value */
/*norm_fac_u=norm(waveconv_u,iter,2);
if(MYID==0){printf("norm_fac_u=%e \n",norm_fac_u);}*/
  
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_u[j][i] = waveconv_u[j][i];
   }
}

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT rho ================================================================================== */
/* ===================================================================================================================================================== */

/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_rho[j][i] = C_rho * waveconv_rho[j][i];
   }
}

/* apply median filter at source positions */
/*median_src(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv_rho);
}

if(SPATFILTER==2){
  smooth2(waveconv_rho);
}
   
/* Normalize gradient to maximum value */
/*norm_fac_rho=norm(waveconv_rho,iter,3);
if(MYID==0){printf("norm_fac_rho=%e \n",norm_fac_rho);}*/

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_rho[j][i] = waveconv_rho[j][i];
   }
} 

/* apply spatial wavelength filter */
/*if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_rho,iter,3);}*/

/* -------------------- */
/* store l-BFGS vectors */
/* -------------------- */

if(iter>1){

   /* load old models and gradients - rho and store them in the LBFGS vectors */
   /* ------------------------------------------------------------------------ */

   sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");
   
   sprintf(jac1,"%s_p_mrho.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   /*iter1 = iter-itershift;*/ /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */
   
   h = NLBFGS_vec*(LBFGS_pointer-1) + 1; /* locate current initial position in LBFGS-vector */
   
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS[h] = waveconv_rho[j][i]-gradplastiter;

	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS[h] = prho[j][i]-modellastiter;
          
          h++;
 
       }
     }
     
     fclose(FP6);
     fclose(FP7);
   
   /* load old models and gradients - Vs and store them in the LBFGS vectors */
   /* ----------------------------------------------------------------------- */
   sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vs.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS[h] = waveconv_u[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS[h] = pu[j][i]-modellastiter;  
          
          h++;
          
       }
     }
     
     fclose(FP6);
     fclose(FP7);

   /* load old models and gradients - Vp and store them in the LBFGS vectors */
   /* ----------------------------------------------------------------------- */
   sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vp.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS[h] = waveconv[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS[h] = ppi[j][i]-modellastiter;  
          
          h++;
          
       }
     }
     
     fclose(FP6);
     fclose(FP7);

     /* save q_LBFGS for all material parameters */    
     h=1;
 
     for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
                             
	     q_LBFGS[h] = waveconv_rho[j][i];
	     h++;
                                                                 
         }
     }                                                                     
                                                                                   
     for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
          
	     q_LBFGS[h] = waveconv_u[j][i];
	     h++;	   
	      
         }
     }

     for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
          
	     q_LBFGS[h] = waveconv[j][i];
	     h++;	   
	      
         }
     }
 
}
	
}
