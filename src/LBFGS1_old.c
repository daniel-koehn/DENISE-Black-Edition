/*------------------------------------------------------------------------
 * Module for the Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
 * method for the elastic multiparameter inversion of
 * vp, vs, rho and lambda, mu, rho respectively
 * 
 * Daniel Koehn
 * ----------------------------------------------------------------------*/

#include "fd.h"

void LBFGS1(float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, int nfstart_jac, float ** waveconv, float C_vp, float ** gradp, float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, 
            float * alpha_LBFGS, float **ppi, float ** pu, float ** prho, int nxnyi, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec){

	extern int NX, NY, IDX, IDY, SPATFILTER;
	extern int HESSIAN, INVMAT, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, k, h, h1, h2;
	float betaz, betan, gradplastiter, gradclastiter, betar, beta;
	float gamma_LBFGS, sum_nom, sum_denom;
        float LBFGSTMP, LBFGSTMP1, LBFGSTMP2, LBFGSTMP3, modellastiter, norm_fac, norm_fac_u, norm_fac_rho;
        int ki, itershift, iter1;
	extern FILE *FP;
	FILE *FP3, *FP4, *FP6, *FP5, *FP7;
	
        itershift = 1;

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Vp/Zp/lambda ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==1)||(INVMAT==0)){
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv[j][i] = C_vp * waveconv[j][i];
   }
}

if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}
 
/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);
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

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

}

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Vs/Zs/mu ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==3)||(INVMAT==0)){
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_u[j][i] = C_vs * waveconv_u[j][i];
   }
}

if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}
 
/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);
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

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_u_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

}

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT rho ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==2)||(INVMAT==0)){

/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_rho[j][i] = C_rho * waveconv_rho[j][i];
   }
}

if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,6);
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

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_rho_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}
}

/* calculate H^-1 * waveconv, using the L-BFGS method, if iter > 1 */
/* --------------------------------------------------------------------- */

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
     
     /* calculate improved first guess Hessian gamma_LBFGS */
     h1 = NLBFGS_vec*(LBFGS_pointer-1) + 1;
     h2 = NLBFGS_vec*LBFGS_pointer; 
     
     sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0);
     sum_denom = dotp(y_LBFGS,y_LBFGS,h1,h2,0);
     gamma_LBFGS = sum_nom/sum_denom;
     
     /*printf("gamma_LBFGS = %e \n",gamma_LBFGS);*/
         
     /* update variable rho for all LBFGS-iterations and all parameter classes*/
     for(k=1;k<=NLBFGS;k++){
          
        h1 = NLBFGS_vec*(k-1) + 1;
        h2 = NLBFGS_vec*k;
        sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0); 
	
	if(sum_nom>0.0){rho_LBFGS[k] = 1.0/sum_nom;}
	else{rho_LBFGS[k] = 0.0;}
	   
	/*if(MYID==0){                                                
	printf("rho_LBFGS = %e of k = %d \n",rho_LBFGS[k],k);}*/
	                                                       
     }
     
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

     /* update alpha_LBFGS and q_LBFGS */
     for(k=1;k<=NLBFGS;k++){
		
       h1 = NLBFGS_vec*(k-1) + 1;
       h2 = NLBFGS_vec*k;
       sum_nom = dotp(s_LBFGS,q_LBFGS,h1,h2,1);
       alpha_LBFGS[k] = rho_LBFGS[k] * sum_nom;
       
       /* update q for all material parameters */
       h = NLBFGS_vec*(k-1) + 1;
       for (i=1;i<=NLBFGS_vec;i++){
           q_LBFGS[i] = q_LBFGS[i] - alpha_LBFGS[k] * y_LBFGS[h];
           h++;
       }
     }
	 
       /* Multiply gradient with approximated Hessian */
       for (i=1;i<=NLBFGS_vec;i++){
           r_LBFGS[i] = gamma_LBFGS * q_LBFGS[i];
       }

     /* calculate H^-1 * waveconv[j][i] */
     for(k=1;k<=NLBFGS;k++){
        
        h1 = NLBFGS_vec*(k-1) + 1;
        h2 = NLBFGS_vec*k;
        /* calculate beta_LBFGS*/   
        sum_nom = dotp(y_LBFGS,r_LBFGS,h1,h2,1);
        beta_LBFGS[k] = rho_LBFGS[k] * sum_nom;

        h = NLBFGS_vec*(k-1) + 1;
        for (i=1;i<=NLBFGS_vec;i++){
	   r_LBFGS[i] = r_LBFGS[i] + s_LBFGS[h]*(alpha_LBFGS[k]-beta_LBFGS[k]);
	   h++;
        }
         
     }

     /* update gradients */
     h=1;
     
     /* density */
     for (i=1;i<=NX;i=i+IDX){   
        for (j=1;j<=NY;j=j+IDY){
                                 
	    waveconv_rho[j][i] = r_LBFGS[h];
	    h++;  
                                                                  
	}
     }
                                                                               
     /* Vs */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            waveconv_u[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }

     /* Vp */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            waveconv[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }


     /* Denormalize Gradients */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
            
           waveconv[j][i] = waveconv[j][i] * C_vp;
	   waveconv_u[j][i] = waveconv_u[j][i] * C_vs;
	   waveconv_rho[j][i] = waveconv_rho[j][i] * C_rho;

        }
     }

}

/* save old models Vs */
/* ------------------ */

    /* save old model */
	sprintf(jac,"%s_p_vp.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&ppi[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vp.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	

/* save old models Vs */
/* ------------------ */

    /* save old model */
	sprintf(jac,"%s_p_vs.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&pu[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vs.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_u[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);


/* save old models Rho */
/* ------------------ */

	sprintf(jac,"%s_p_mrho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&prho[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_mrho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_rho[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g_rho */
        sprintf(jac,"%s_c_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
}
