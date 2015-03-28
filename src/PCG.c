/*------------------------------------------------------------------------
 * Module for the Preconditioned Conjugate Gradient Method (PCG)
 * for the material parameters vp, vs, rho and lambda, mu, rho respectively
 * 
 * Daniel Koehn
 * ----------------------------------------------------------------------*/

#include "fd.h"

void PCG(float ** waveconv, float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float C_vp, float ** gradp, int nfstart_jac,
	     float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho){

	extern int NX, NY, IDX, IDY, SPATFILTER, GRAD_FILTER;
	extern int HESSIAN, INVMAT, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	extern int POS[3], MYID, GRAD_METHOD;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac2[225];
	int i, j;
	float betaz, betan, gradplastiter, gradclastiter, betar, beta;
	extern FILE *FP;
	FILE *FP3, *FP4, *FP6, *FP5;
	
/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT ZP ================================================================================== */
/* ===================================================================================================================================================== */


if((HESSIAN!=1)&&(INVMAT==0)){
	
/* Preconditioning of the gradient */
/* ------------------------------- */

/* apply taper on the gradient */
/* --------------------------- */

if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}

if (SWS_TAPER_FILE){   /* read taper from BIN-File*/
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,4);}   

/* apply median filter at source positions */
/*median_src(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv);
}

if(SPATFILTER==2){
  smooth2(waveconv);
}
  
/* save gradient */
/*sprintf(jac,"%s_g.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv[j][i],sizeof(float),1,FP3);
           }
        }
	
fclose(FP3);

MPI_Barrier(MPI_COMM_WORLD);*/
          
/* merge gradient file */ 
/*sprintf(jac,"%s_g.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);*/
 
/* Normalize gradient to maximum value */
/*norm(waveconv,iter,1);*/

/* apply spatial wavelength filter */
/*if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
        spat_filt(waveconv,iter,1);}*/

/* apply 2D-Gaussian filter*/
if(GRAD_FILTER==1){smooth_grad(waveconv,1);}

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv[j][i] = C_vp * waveconv[j][i];
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
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_p_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	remove(jac);
}



/* calculate conjugate gradient direction, if iter > 1 (after Mora 1987) */
/* --------------------------------------------------------------------- */
if(GRAD_METHOD!=3){

if(iter>1){
   
   sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   if(iter>2){
      sprintf(jac2,"%s_c.old.%i%i",JACOBIAN,POS[1],POS[2]);
      FP5=fopen(jac2,"rb");
   }
   
        /* apply scalar product to obtain the coefficient beta */
     betaz = 0.0;
     betan = 0.0;
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          fread(&gradplastiter,sizeof(float),1,FP6);
	  
	  /*if(gradglastiter==gradg[j][i]) err("TEST1");*/
	  /*if (MYID==10)  printf("TEST beta (MYID=%d) bei (j,i)=(%i,%i): gradg(k-1) = %e, gradg(k) = %e\n",MYID,j,i,gradglastiter,gradg[j][i]);*/
	  
	  /*
	  betaz += (1e5*gradp[j][i]) * ( (1e5*gradg[j][i]) - (1e5*gradglastiter) );
	  betan += (1e5*gradplastiter) * (1e5*gradglastiter);
	  */
	  
	  /* Polak and Ribiere */
	  /*betaz += (gradp[j][i]) * ( (gradg[j][i]) - (gradglastiter) );
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  /* Polak and Ribiere */
	  betaz += (gradp[j][i]) * ( (gradp[j][i]) - (gradplastiter) );
	  betan += (gradplastiter) * (gradplastiter);
	  
	  /* Fletcher and Reeves */
	  /*betaz += (gradp[j][i]) * (gradg[j][i]);
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  
       }
     }
     
     /*printf("TEST: vor exchange (MYID=%d): beta = betaz/betan = %e/%e = %e\n",MYID,betaz,betan,betaz/betan);*/

     /*betaz = exchange_L2(betaz,1,1);
     betan = exchange_L2(betan,1,1);*/
     
     betar = 0.0;
     MPI_Allreduce(&betaz,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betaz = betar;
     
     betar = 0.0;
     MPI_Allreduce(&betan,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betan = betar;
     
     beta = betaz/betan;
     
     /* direction reset */
     if(beta<0.0){beta = 0.0;}

     /*betaVp = beta;*/
     
     printf("\n\nTEST: nach exchange (MYID=%d): beta = %e / %e = %e\n",MYID,betaz,betan,beta);
     
     fseek(FP6,0,SEEK_SET);
     
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	
	  if(iter==2){
             fread(&gradplastiter,sizeof(float),1,FP6);
             waveconv[j][i] = gradp[j][i] + gradplastiter * beta;
          }
   
          if(iter>2){
	     fread(&gradclastiter,sizeof(float),1,FP5);
             waveconv[j][i] = gradp[j][i] + gradclastiter * beta;
          }
	  
       }
     }
     
   fclose(FP6);
   
   if(iter>2){fclose(FP5);}

}

/* output of the conjugate gradient */
if(iter>1){
  sprintf(jac2,"%s_c.old.%i%i",JACOBIAN,POS[1],POS[2]);
  FP5=fopen(jac2,"wb");


  for (i=1;i<=NX;i=i+IDX){
     for (j=1;j<=NY;j=j+IDY){
         fwrite(&waveconv[j][i],sizeof(float),1,FP5);
     }
  }

  fclose(FP5);
  
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac2,"%s_c.old",JACOBIAN);
if (MYID==0) mergemod(jac2,3);  

}

} /* end of if GRAD_METHOD!=3*/

/* output of preconditioned gradient */
sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        /*fwrite(&waveconv[j][i],sizeof(float),1,FP4);*/
	fwrite(&gradp[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_p.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);

}

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Zs ================================================================================== */
/* ===================================================================================================================================================== */

if((HESSIAN!=1)&&(INVMAT==0)){
	
/* Preconditioning of the gradient */
/* ------------------------------- */

/* apply taper on the gradient */
/* --------------------------- */
if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}

if (SWS_TAPER_FILE){   /* read taper from BIN-File*/                          
   taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);}

/* apply median filter at source positions */
/*median_src(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv_u);
}

if(SPATFILTER==2){
  smooth2(waveconv_u);
}
  
/* save gradient */
/*sprintf(jac,"%s_g_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
           }
        }
	
fclose(FP3);

MPI_Barrier(MPI_COMM_WORLD);*/
          
/* merge gradient file */ 
/*sprintf(jac,"%s_g_u.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);*/
 

/* Normalize gradient to maximum value */
/*norm(waveconv_u,iter,2);*/

/* apply spatial wavelength filter */
/*if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_u,iter,2);}*/

/* apply 2D-Gaussian filter*/
if(GRAD_FILTER==1){smooth_grad(waveconv_u,2);}

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_u[j][i] = C_vs * waveconv_u[j][i];
      gradp_u[j][i]=waveconv_u[j][i];
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
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_p_u_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	remove(jac);
}


/* calculate conjugate gradient direction, if iter > 1 (after Mora 1987) */
/* --------------------------------------------------------------------- */
if(GRAD_METHOD!=3){

if(iter>1){

   
   sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   if(iter>2){
      sprintf(jac2,"%s_c_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
      FP5=fopen(jac2,"rb");
   }
   
        /* apply scalar product to obtain the coefficient beta */
     betaz = 0.0;
     betan = 0.0;
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          fread(&gradplastiter,sizeof(float),1,FP6);
	  
	  /*if(gradglastiter==gradg[j][i]) err("TEST1");*/
	  /*if (MYID==10)  printf("TEST beta (MYID=%d) bei (j,i)=(%i,%i): gradg(k-1) = %e, gradg(k) = %e\n",MYID,j,i,gradglastiter,gradg[j][i]);*/
	  
	  /*
	  betaz += (1e5*gradp[j][i]) * ( (1e5*gradg[j][i]) - (1e5*gradglastiter) );
	  betan += (1e5*gradplastiter) * (1e5*gradglastiter);
	  */
	  
	  /* Polak and Ribiere */
	  /*betaz += (gradp_u[j][i]) * ( (gradg_u[j][i]) - (gradglastiter) );
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  /* Polak and Ribiere */
	  betaz += (gradp_u[j][i]) * ( (gradp_u[j][i]) - (gradplastiter) );
	  betan += (gradplastiter) * (gradplastiter);
	  
	  /* Fletcher and Reeves */
	  /*betaz += (gradp[j][i]) * (gradg[j][i]);
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  
       }
     }
     
     /*printf("TEST: vor exchange (MYID=%d): beta = betaz/betan = %e/%e = %e\n",MYID,betaz,betan,betaz/betan);*/

     /*betaz = exchange_L2(betaz,1,1);
     betan = exchange_L2(betan,1,1);*/
     
     betar = 0.0;
     MPI_Allreduce(&betaz,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betaz = betar;
     
     betar = 0.0;
     MPI_Allreduce(&betan,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betan = betar;
     
     beta = betaz/betan;
     
     /* direction reset */
     if(beta<0.0){beta = 0.0;}

     /*betaVs = beta;*/
     printf("\n\nTEST: nach exchange (MYID=%d): beta = %e / %e = %e\n",MYID,betaz,betan,beta);
     
     fseek(FP6,0,SEEK_SET);
     
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	
	  if(iter==2){
             fread(&gradplastiter,sizeof(float),1,FP6);
             waveconv_u[j][i] = gradp_u[j][i] + gradplastiter * beta;
          }
   
          if(iter>2){
	     fread(&gradclastiter,sizeof(float),1,FP5);
             waveconv_u[j][i] = gradp_u[j][i] + gradclastiter * beta;
          }

       }
     }
     
   fclose(FP6);
   
   if(iter>2){fclose(FP5);}

}

/* output of the conjugate gradient */
if(iter>1){
  sprintf(jac2,"%s_c_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
  FP5=fopen(jac2,"wb");


  for (i=1;i<=NX;i=i+IDX){
     for (j=1;j<=NY;j=j+IDY){
         fwrite(&waveconv_u[j][i],sizeof(float),1,FP5);
     }
  }

  fclose(FP5);
  
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac2,"%s_c_u.old",JACOBIAN);
if (MYID==0) mergemod(jac2,3);  

}

} /* end of GRAD_METHOD!=3*/

sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        /*fwrite(&waveconv_u[j][i],sizeof(float),1,FP4);*/
	fwrite(&gradp_u[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_p_u.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);

}

/* =================================================================================================================================================== */

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT rho ================================================================================== */
/* ===================================================================================================================================================== */

if((HESSIAN!=1)&&(INVMAT==0)){
	
/* Preconditioning of the gradient */
/* ------------------------------- */
if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,1);}

if (SWS_TAPER_GRAD_HOR){   /*horizontal gradient taper is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,2);}

if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,3);}

if (SWS_TAPER_FILE){   /* read taper from BIN-File*/                          
   taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,6);}

/* apply median filter at source positions */
/*median_src(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,iter,0);*/

/* apply wavenumber damping */
if(SPATFILTER==1){
  wavenumber(waveconv_rho);
}

if(SPATFILTER==2){
  smooth2(waveconv_rho);
} 

/* save gradient */
/*sprintf(jac,"%s_g_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
           }
        }
	
fclose(FP3);

MPI_Barrier(MPI_COMM_WORLD);*/
          
/* merge gradient file */ 
/*sprintf(jac,"%s_g_rho.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);*/
 

/* Normalize gradient to maximum value */
/*norm(waveconv_rho,iter,3);*/


/* apply spatial wavelength filter */
/*if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_rho,iter,3);}*/

/* apply 2D-Gaussian filter*/
if(GRAD_FILTER==1){smooth_grad(waveconv_rho,3);}

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_rho[j][i] = C_rho * waveconv_rho[j][i];
	gradp_rho[j][i]=waveconv_rho[j][i];
   }
}



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
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_p_rho_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	remove(jac);
}



/* calculate conjugate gradient direction, if iter > 1 (after Mora 1987) */
/* --------------------------------------------------------------------- */

if(GRAD_METHOD!=3){

if(iter>1){

   
   sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   if(iter>2){
      sprintf(jac2,"%s_c_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
      FP5=fopen(jac2,"rb");
   }
   
        /* apply scalar product to obtain the coefficient beta */
     betaz = 0.0;
     betan = 0.0;
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          fread(&gradplastiter,sizeof(float),1,FP6);
	  
	  /*if(gradglastiter==gradg[j][i]) err("TEST1");*/
	  /*if (MYID==10)  printf("TEST beta (MYID=%d) bei (j,i)=(%i,%i): gradg(k-1) = %e, gradg(k) = %e\n",MYID,j,i,gradglastiter,gradg[j][i]);*/
	  
	  /*
	  betaz += (1e5*gradp[j][i]) * ( (1e5*gradg[j][i]) - (1e5*gradglastiter) );
	  betan += (1e5*gradplastiter) * (1e5*gradglastiter);
	  */
	  
	  /* Polak and Ribiere */
	  /*betaz += (gradp_rho[j][i]) * ( (gradg_rho[j][i]) - (gradglastiter) );
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  /* Polak and Ribiere */
	  betaz += (gradp_rho[j][i]) * ( (gradp_rho[j][i]) - (gradplastiter) );
	  betan += (gradplastiter) * (gradplastiter);
	  
	  /* Fletcher and Reeves */
	  /*betaz += (gradp[j][i]) * (gradg[j][i]);
	  betan += (gradplastiter) * (gradglastiter);*/
	  
	  
       }
     }
     
     /*printf("TEST: vor exchange (MYID=%d): beta = betaz/betan = %e/%e = %e\n",MYID,betaz,betan,betaz/betan);*/

     /*betaz = exchange_L2(betaz,1,1);
     betan = exchange_L2(betan,1,1);*/
     
     betar = 0.0;
     MPI_Allreduce(&betaz,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betaz = betar;
     
     betar = 0.0;
     MPI_Allreduce(&betan,&betar,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     betan = betar;
     
     beta = betaz/betan;
     
     /* direction reset */
     if(beta<0.0){beta = 0.0;}

     /*betarho = beta;*/
     printf("\n\nTEST: nach exchange (MYID=%d): beta = %e / %e = %e\n",MYID,betaz,betan,beta);
     
     fseek(FP6,0,SEEK_SET);
     
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	
	  if(iter==2){
             fread(&gradplastiter,sizeof(float),1,FP6);
             waveconv_rho[j][i] = gradp_rho[j][i] + gradplastiter * beta;
          }
   
          if(iter>2){
	     fread(&gradclastiter,sizeof(float),1,FP5);
             waveconv_rho[j][i] = gradp_rho[j][i] + gradclastiter * beta;
          }
	  
	  
	  /*if (iter >= 2)
	  {
	     if (isnan(waveconv_u[j][i]) || isinf(waveconv_u[j][i]))
	     {
		     sum = 0.0;
		     h = 0;
		     for (ii=-1;ii<=1;ii++){
	       		for (jj=-1;jj<=1;jj++){
				if (isnan(waveconv_rho[j+jj][i+ii]) || isinf(waveconv_rho[j+jj][i+ii])) continue;
				sum += waveconv_rho[j+jj][i+ii];
				h++;
			}
		     }
		     if (h>0) waveconv_rho[j][i] = sum / h;
		     else waveconv_rho[j][i] = 0.0;
	     }
	     
	  }*/

       }
     }
     
   fclose(FP6);
   
   if(iter>2){fclose(FP5);}

}

/* output of the conjugate gradient */
if(iter>1){
  sprintf(jac2,"%s_c_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
  FP5=fopen(jac2,"wb");


  for (i=1;i<=NX;i=i+IDX){
     for (j=1;j<=NY;j=j+IDY){
         fwrite(&waveconv_rho[j][i],sizeof(float),1,FP5);
     }
  }

  fclose(FP5);
  
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac2,"%s_c_rho.old",JACOBIAN);
if (MYID==0) mergemod(jac2,3);  

}

} /* end of GRAD_METHOD!=3 */

sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the preconditioned gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        /*fwrite(&waveconv_rho[j][i],sizeof(float),1,FP4);*/
	fwrite(&gradp_rho[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_p_rho.old",JACOBIAN);
if (MYID==0) mergemod(jac,3);

}

	
}
