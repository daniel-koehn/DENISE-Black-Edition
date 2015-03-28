/*
 * Smoothing model with a median frame to damp inversion artefacts near the sources and receivers
 * M. Schaefer - more or less copied from spat_filt.c - May 2011
 */

#include "fd.h"

void smooth_grad(float ** waveconv, int sws)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	extern int FILT_SIZE_GRAD;
	
	/* local variables */
	int i, j, ii, jj;
	int i1, j1, filtsize, hfs;

	float **model_tmp, **model_med, **filterpart, grad, normgauss, smooth_meter;
	
	char jac_tmp[STRING_SIZE];
	
	FILE *model;
	
	char modfile[STRING_SIZE];
	
	if(MYID==0){
		      	if (FILT_SIZE_GRAD==0)	return;
		      	if (!(FILT_SIZE_GRAD % 2)) {
		      	if (FILT_SIZE_GRAD > 0)	FILT_SIZE_GRAD += 1;
		      	else			FILT_SIZE_GRAD -= 1;
		      	}
	  
		      	hfs = abs(FILT_SIZE_GRAD)/2;
		      	printf("\n hfs: %d \n",hfs);
		      	model_med = matrix(1,NYG,1,NXG);
		      	model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);

		      	filterpart=matrix(1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
		      
		      
		    	if(sws==1){
			sprintf(jac_tmp,"%s_g.old",JACOBIAN);}   

			if(sws==2){
			sprintf(jac_tmp,"%s_g_u.old",JACOBIAN);}    
	
			if(sws==3){
			sprintf(jac_tmp,"%s_g_rho.old",JACOBIAN);}    
	
		      	model=fopen(jac_tmp,"rb");
		      	if (model==NULL) err(" Could not open jacobian file !");
		
		      	/* load merged model */
		      	for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					      fread(&grad, sizeof(float), 1, model);
				      	model_tmp[j][i]=grad;
			      	}	
		      	}
		
		      fclose(model);
		      
		      
		      
		      /* apply 2D-Gaussian filter on vp and vs model */
			      /* extrapolate array */
			      /* left/right boundary */
			      for (j=1;j<=NYG;j++){
				      for (i=-hfs+1;i<=0;i++){
					model_tmp[j][i] = model_tmp[j][1];}
				      for (i=NXG+1;i<=NXG+hfs;i++){
					model_tmp[j][i] = model_tmp[j][NXG];}
			      }
			      /* top/bottom boundary incl. corners */
			      for (j=-hfs+1;j<=0;j++){
				      for (i=-hfs+1;i<=NXG+hfs;i++){
					model_tmp[j][i] = model_tmp[1][i];}
			      }
			      for (j=NYG+1;j<=NYG+hfs;j++){
				      for (i=-hfs+1;i<=NXG+hfs;i++){
					model_tmp[j][i] = model_tmp[NYG][i];}
			      }
		
			/* filter */
			      for (j=1;j<=NYG;j++){
				for (i=1;i<=NXG;i++){
					      /* create a filtersize x filtersize matrix */
					      for (ii=-hfs;ii<=hfs;ii++){
						      for (jj=-hfs;jj<=hfs;jj++){
							      /*if ((ii+hfs+1)<(-hfs+1)) err(" (ii+hfs+1)<(-hfs+1) ! ");
							      if ((ii+hfs+1)>(hfs+NXG)) err(" (ii+hfs+1)>(hfs+NXG) ! ");*/
							      filterpart[jj+hfs+1][ii+hfs+1] = model_tmp[j+jj][i+ii];
						      }
						}
					      /* filter */
						model_med[j][i] = median2d(filterpart,abs(FILT_SIZE_GRAD),abs(FILT_SIZE_GRAD));				
				      }
			      }
			      
			/* output of the preconditioned jacobian */      
			if(sws==1){
			sprintf(jac_tmp,"%s_tmp_g.old",JACOBIAN);}   

			if(sws==2){
			sprintf(jac_tmp,"%s_tmp_g_u.old",JACOBIAN);}    
	
			if(sws==3){
			sprintf(jac_tmp,"%s_tmp_g_rho.old",JACOBIAN);}    
			      
			model=fopen(jac_tmp,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			  fwrite(&model_med[j][i],sizeof(float),1,model);

			  }
			}
			fclose(model);
			
			free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(model_med,1,NXG,1,NYG);
			free_matrix(filterpart,1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
			
		}/* end of if(MYID==0)*/
	MPI_Barrier(MPI_COMM_WORLD);
	smooth_meter=FILT_SIZE_GRAD*DH;
	fprintf(FP,"\n \t ---- Gradient is smoothed with filter length of %d gridpoints which is equivalent to %4.2f meter---- \n",FILT_SIZE_GRAD,smooth_meter);
	
	/* distribute smoothed jacobian on computational nodes */
	if(sws==1){
	sprintf(jac_tmp,"%s_tmp_g.old",JACOBIAN);}   

	if(sws==2){
	sprintf(jac_tmp,"%s_tmp_g_u.old",JACOBIAN);}    
	
	if(sws==3){
	sprintf(jac_tmp,"%s_tmp_g_rho.old",JACOBIAN);}
		
	model=fopen(jac_tmp,"rb");
	if (model==NULL) err(" Could not open jacobian file ! (distribute smoothed jacobian)");
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model);
			   			
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				waveconv[jj][ii]=grad;

			}
		}
	}
		
	fclose(model);
	
	fprintf(FP,"\n \t ---- Smoothed jacobian is distributed on computational nodes ... ---- \n");

}/* end of application condition for the smoothing */
