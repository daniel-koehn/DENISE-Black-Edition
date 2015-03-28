/*
 * Smoothing model with a median frame to damp inversion artefacts near the sources and receivers
 * M. Schaefer - more or less copied from spat_filt.c - May 2011
 */

#include "fd.h"

void smooth_model(float ** ppi, float ** pu, float ** prho, int iter)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char INV_MODELFILE[STRING_SIZE];
	extern FILE *FP;
	extern int MODEL_FILTER, FILT_SIZE;
	
	/* local variables */
	int i, j, ii, jj;
	int i1, j1, filtsize, hfs;

	float **model_tmp_vp, **model_med_vp, **model_tmp_vs, **model_med_vs, **model_tmp_rho, **model_med_rho, **filterpart_vp, **filterpart_rho, **filterpart_vs, grad, normgauss, smooth_meter;
	
	char jac_vp[STRING_SIZE];
	char jac_vs[STRING_SIZE];
	char jac_rho[STRING_SIZE];
	
	FILE *model_vp;
	FILE *model_vs;
	FILE *model_rho;
	
	char modfile_vp[STRING_SIZE];
	char modfile_vs[STRING_SIZE];
	char modfile_rho[STRING_SIZE];
	
	if(MODEL_FILTER==1){
	      if(MYID==0){
		      if (FILT_SIZE==0)	return;
		      if (!(FILT_SIZE % 2)) {
		      if (FILT_SIZE > 0)	FILT_SIZE += 1;
		      else			FILT_SIZE -= 1;
		      }
	  
		      hfs = abs(FILT_SIZE)/2;
		      printf("\n hfs: %d \n",hfs);
		      model_med_vp = matrix(1,NYG,1,NXG);
		      model_tmp_vp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);

		      model_med_vs = matrix(1,NYG,1,NXG);
		      model_tmp_vs = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
		      
		      model_med_rho = matrix(1,NYG,1,NXG);
		      model_tmp_rho = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
		
		      filterpart_vp=matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
		      filterpart_vs=matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
		      filterpart_rho=matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
		      
		      sprintf(jac_vp,"%s_vp.bin",INV_MODELFILE);
		      model_vp=fopen(jac_vp,"rb");
		      if (model_vp==NULL) err(" Could not open vp model file ! ");
		
		      /* load merged model */
		      for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				      fread(&grad, sizeof(float), 1, model_vp);
				      model_tmp_vp[j][i]=grad;
			      }
		      }
		
		      fclose(model_vp);
		      
		      sprintf(jac_vs,"%s_vs.bin",INV_MODELFILE);
		      model_vs=fopen(jac_vs,"rb");
		      if (model_vs==NULL) err(" Could not open vs model file ! ");
		
		      /* load merged model */
		      for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				      fread(&grad, sizeof(float), 1, model_vs);
				      model_tmp_vs[j][i]=grad;
			      }
		      }
		
		      fclose(model_vs);
		      
		      sprintf(jac_rho,"%s_rho.bin",INV_MODELFILE);
		      model_rho=fopen(jac_rho,"rb");
		      if (model_rho==NULL) err(" Could not open rho model file ! ");
		
		      /* load merged model */
		      for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				      fread(&grad, sizeof(float), 1, model_rho);
				      model_tmp_rho[j][i]=grad;
			      }
		      }
		
		      fclose(model_rho);
		      
		      /* apply 2D-Gaussian filter on vp and vs model */
			      /* extrapolate array */
			      /* left/right boundary */
			      for (j=1;j<=NYG;j++){
				      for (i=-hfs+1;i<=0;i++){
					model_tmp_vp[j][i] = model_tmp_vp[j][1];
					model_tmp_vs[j][i] = model_tmp_vs[j][1];
					model_tmp_rho[j][i] = model_tmp_rho[j][1];}
				      for (i=NXG+1;i<=NXG+hfs;i++){
					model_tmp_vp[j][i] = model_tmp_vp[j][NXG];
					model_tmp_vs[j][i] = model_tmp_vs[j][NXG];
					model_tmp_rho[j][i] = model_tmp_rho[j][NXG];}
			      }
			      /* top/bottom boundary incl. corners */
			      for (j=-hfs+1;j<=0;j++){
				      for (i=-hfs+1;i<=NXG+hfs;i++){
					model_tmp_vp[j][i] = model_tmp_vp[1][i];
					model_tmp_vs[j][i] = model_tmp_vs[1][i];
					model_tmp_rho[j][i] = model_tmp_rho[1][i];}
			      }
			      for (j=NYG+1;j<=NYG+hfs;j++){
				      for (i=-hfs+1;i<=NXG+hfs;i++){
					model_tmp_vp[j][i] = model_tmp_vp[NYG][i];
					model_tmp_vs[j][i] = model_tmp_vs[NYG][i];
					model_tmp_rho[j][i] = model_tmp_rho[NYG][i];}
			      }
		
			/* filter */
			      for (j=1;j<=NYG;j++){
				for (i=1;i<=NXG;i++){
					      /* create a filtersize x filtersize matrix */
					      for (ii=-hfs;ii<=hfs;ii++){
						      for (jj=-hfs;jj<=hfs;jj++){
							      /*if ((ii+hfs+1)<(-hfs+1)) err(" (ii+hfs+1)<(-hfs+1) ! ");
							      if ((ii+hfs+1)>(hfs+NXG)) err(" (ii+hfs+1)>(hfs+NXG) ! ");*/
							      filterpart_vp[jj+hfs+1][ii+hfs+1] = model_tmp_vp[j+jj][i+ii];
							      filterpart_vs[jj+hfs+1][ii+hfs+1] = model_tmp_vs[j+jj][i+ii];
							      filterpart_rho[jj+hfs+1][ii+hfs+1] = model_tmp_rho[j+jj][i+ii];
						      }
						}
					      /* filter */
						model_med_vp[j][i] = median2d(filterpart_vp,abs(FILT_SIZE),abs(FILT_SIZE));
						model_med_vs[j][i] = median2d(filterpart_vs,abs(FILT_SIZE),abs(FILT_SIZE));
						model_med_rho[j][i] = median2d(filterpart_rho,abs(FILT_SIZE),abs(FILT_SIZE));				
				      }
			      }
			      
			/* output of the preconditioned model */      
			sprintf(jac_vp,"%s_vp_tmp.bin",INV_MODELFILE);      
			model_vp=fopen(jac_vp,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			  fwrite(&model_med_vp[j][i],sizeof(float),1,model_vp);

			  }
			}
			fclose(model_vp);
			
			/* output of the preconditioned model */
			sprintf(jac_vs,"%s_vs_tmp.bin",INV_MODELFILE);
			model_vs=fopen(jac_vs,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			  fwrite(&model_med_vs[j][i],sizeof(float),1,model_vs);

			  }
			}
			fclose(model_vs);
			
			/* output of the preconditioned model */
			sprintf(jac_rho,"%s_rho_tmp.bin",INV_MODELFILE);
			model_rho=fopen(jac_rho,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			  fwrite(&model_med_rho[j][i],sizeof(float),1,model_rho);

			  }
			}
			fclose(model_rho);
			
			free_matrix(model_tmp_vs,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(model_med_vs,1,NXG,1,NYG);
	
			free_matrix(model_tmp_vp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(model_med_vp,1,NXG,1,NYG);
			
			free_matrix(model_tmp_rho,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(model_med_rho,1,NXG,1,NYG);
			
			free_matrix(filterpart_vp,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
			free_matrix(filterpart_vs,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
			free_matrix(filterpart_rho,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
			
		}/* end of if(MYID==0)*/
	
	
	smooth_meter=FILT_SIZE*DH;
	fprintf(FP,"\n \t ---- V_p, V_s and rho models are smoothed with filter length of %d gridpoints which is equivalent to %4.2f meter---- \n",FILT_SIZE,smooth_meter);
	/*printf(" MYID = %d , Models are smoothed ... \n",MYID);*/
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* distribute smoothed models on computational nodes */
	sprintf(jac_vp,"%s_vp_tmp.bin",INV_MODELFILE);
	model_vp=fopen(jac_vp,"rb");
	if (model_vp==NULL) err(" Could not open vp model file ! ");
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model_vp);
			   			
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				ppi[jj][ii]=grad;

			}
		}
	}
		
	fclose(model_vp);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* distribute smoothed models on computational nodes */
	sprintf(jac_vs,"%s_vs_tmp.bin",INV_MODELFILE);
	model_vs=fopen(jac_vs,"rb");
	if (model_vs==NULL) err(" Could not open vs model file ! ");
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model_vs);
			   			
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				pu[jj][ii]=grad;

			}
		}
	}
		
	fclose(model_vs);
	
	MPI_Barrier(MPI_COMM_WORLD); 
	
	/* distribute smoothed models on computational nodes */
	sprintf(jac_rho,"%s_rho_tmp.bin",INV_MODELFILE);
	model_rho=fopen(jac_rho,"rb");
	if (model_rho==NULL) err(" Could not open rho model file ! ");
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model_rho);
			   			
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				prho[jj][ii]=grad;

			}
		}
	}
		
	fclose(model_rho);
	
	MPI_Barrier(MPI_COMM_WORLD); 

	fprintf(FP,"\n \t ---- Smoothed models are distributed on computational nodes ... ---- \n");

        /*sprintf(modfile_vp,"%s_vp_smoothed_it_%d.bin",INV_MODELFILE,iter);
	writemod(modfile_vp,ppi,3);
                                        
	MPI_Barrier(MPI_COMM_WORLD);
                                                
	if (MYID==0) mergemod(modfile_vp,3);
                                                        
	sprintf(modfile_vs,"%s_vs_smoothed_it_%d.bin",INV_MODELFILE,iter);
	writemod(modfile_vs,pu,3);
	
	MPI_Barrier(MPI_COMM_WORLD);
                                                                                                        
	if (MYID==0) mergemod(modfile_vs,3); 
	
	sprintf(modfile_rho,"%s_rho_smoothed_it_%d.bin",INV_MODELFILE,iter);
	writemod(modfile_rho,prho,3);
                                        
	MPI_Barrier(MPI_COMM_WORLD);
                                                
	if (MYID==0) mergemod(modfile_rho,3);                                                                                                              

	MPI_Barrier(MPI_COMM_WORLD); */
	
	
	}/* end of if(MODEL_FILTER==1)*/
}/* end of application condition for the smoothing */

