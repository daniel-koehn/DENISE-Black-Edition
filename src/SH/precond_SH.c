/*------------------------------------------------------------------------
 *  Preconditioning of gradients after shot summation (SH problem)
 * 
 *  Daniel Koehn
 *  Kiel, 13/12/2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void precond_SH(struct fwiSH *fwiSH, struct acq *acq, int nsrc, int ntr_glob, float ** taper_coeff, FILE *FP_GRAV){ 
		
        /* global variables */
	extern int SEISMO, MYID, NX, NY, IDX, IDY, POS[3], L;
        extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
        extern char JACOBIAN[STRING_SIZE];
	
        /* local variables */
        int i, j;
        char jac[STRING_SIZE];
        float tmp;


	/*================== TAPER Vs/Zs/mu ===========================*/
	if (SWS_TAPER_GRAD_VERT){    /* vertical gradient taper is applied */
	   taper_grad((*fwiSH).waveconv_u,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,1);}

	if (SWS_TAPER_GRAD_HOR){    /* horizontal gradient taper is applied */
	   taper_grad((*fwiSH).waveconv_u,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,2);}

	if(SWS_TAPER_GRAD_SOURCES){    /* cylindrical taper around sources is applied */
	   taper_grad((*fwiSH).waveconv_u,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,3);}

	if(SWS_TAPER_FILE){ /* read taper from file */
	   taper_grad((*fwiSH).waveconv_u,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,5);}


	/*================== TAPER Rho ===========================*/
	if (SWS_TAPER_GRAD_VERT){    /* vertical gradient taper is applied */
	   taper_grad((*fwiSH).waveconv_rho,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,1);}

	if (SWS_TAPER_GRAD_HOR){     /* horizontal gradient taper is applied */
	   taper_grad((*fwiSH).waveconv_rho,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,2);}

	if (SWS_TAPER_GRAD_SOURCES){    /* cylindrical taper around sources is applied */
	   taper_grad((*fwiSH).waveconv_rho,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,3);}

	if (SWS_TAPER_FILE){ /* read taper from file */
	   taper_grad((*fwiSH).waveconv_rho,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,6);}


	/*================== TAPER Taus ===========================*/

	if(L){

	   if (SWS_TAPER_GRAD_VERT){    /* vertical gradient taper is applied */
	      taper_grad((*fwiSH).waveconv_ts,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,1);}

	   if (SWS_TAPER_GRAD_HOR){     /* horizontal gradient taper is applied */
	      taper_grad((*fwiSH).waveconv_ts,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,2);}

	   if (SWS_TAPER_GRAD_SOURCES){    /* cylindrical taper around sources is applied */
	      taper_grad((*fwiSH).waveconv_ts,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,3);}

	   if (SWS_TAPER_FILE){ /* read taper from file */
	      taper_grad((*fwiSH).waveconv_ts,taper_coeff,(*acq).srcpos,nsrc,(*acq).recpos,ntr_glob,6);}

	}

	/* output of the seismic gradient for rho after taper  */
	/*sprintf(jac,"%s_seis.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP_GRAV=fopen(jac,"wb");       

	for (i=1;i<=NX;i=i+IDX){
	    for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiPSV).waveconv_rho[j][i];
		fwrite(&tmp,sizeof(float),1,FP_GRAV);
	    }
	}

	fclose(FP_GRAV);

	MPI_Barrier(MPI_COMM_WORLD);*/

	/* merge model file */
        /*sprintf(jac,"%s_seis",JACOBIAN);          
	if (MYID==0) mergemod(jac,3); */
		
}
