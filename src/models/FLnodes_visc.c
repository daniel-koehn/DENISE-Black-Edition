/*
 *   Model homogeneous half space
 *   last update 11.04.02, T. Bohlen
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH, *FL, TAU, DT;
		/* local variables */
	float vp, vs, rhov, ts, tp, muv, piv, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE]; 
	
	FILE *flfile;
	int nodes;
	char cline[256];
	
	float *fldepth, *flrho, *flvp, *flvs;
	
	
	/*-----------------------------------------------------------------------*/
	
	nodes=7;

	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);
	flvs=vector(1,nodes);
	
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	        eta[l]=DT/pts[l];
	}
	
	/*read FL nodes from File*/
	
	flfile=fopen("model_true/model4.fl.dat","r");
	if (flfile==NULL) err(" FL-file could not be opened !");
	
	
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &flvs[l]);
		}
		else l=l-1;
	
	}
	
	if(MYID==0){
	printf(" ------------------------------------------------------------------ \n\n");
	printf(" Information of FL nodes: \n\n");
	printf(" \t depth \t rho \t vp \t vs \n\n");
	
	for (l=1;l<=nodes;l++){
	printf(" \t %f \t %f \t %f \t %f \n\n",fldepth[l],flrho[l],flvp[l],flvs[l]);
	}
	printf(" ------------------------------------------------------------------ \n\n");
	}
	/*-----------------------------------------------------------------------*/
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (l=1;l<nodes;l++){
				for(j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){	
					
				  
					vp=0.0;
					vs=0.0; 
					rhov=0.0;
				  
					vp=(DH*(j-1)-fldepth[l])*(flvp[l+1]-flvp[l])/(fldepth[l+1]-fldepth[l])+flvp[l];
					vp=vp*1000.0;
					vs=(DH*(j-1)-fldepth[l])*(flvs[l+1]-flvs[l])/(fldepth[l+1]-fldepth[l])+flvs[l];
					vs=vs*1000.0;
					rhov=(DH*(j-1)-fldepth[l])*(flrho[l+1]-flrho[l])/(fldepth[l+1]-fldepth[l])+flrho[l];
					rhov=rhov*1000.0; 				

					muv=vs;
					piv=vp;
					ts=TAU;
					tp=TAU;
					
					/* only the PE which belongs to the current global gridpoint 
					  is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) && 
					    (POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						

						u[jj][ii]=muv;
						rho[jj][ii]=rhov;
						pi[jj][ii]=piv;
						taus[jj][ii]=ts;
						taup[jj][ii]=tp;
					}
			     	}
			
			
				for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			  
				vp=0.0; vs=0.0; rhov=0.0;
				vp=flvp[nodes]*1000.0; vs=flvs[nodes]*1000.0; rhov=flrho[nodes]*1000.0;
				
				muv=vs;
				piv=vp;
				ts=TAU;
				tp=TAU;

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
						
					u[jj][ii]=muv;
					rho[jj][ii]=rhov;
					pi[jj][ii]=piv;
					taus[jj][ii]=ts;
					taup[jj][ii]=tp;
					
				}
				}

			}
		}	

		
sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
writemod(modfile,rho,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_0.bin",INV_MODELFILE);
writemod(modfile,u,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vp_it_0.bin",INV_MODELFILE);
writemod(modfile,pi,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

free_vector(fldepth,1,nodes);
free_vector(flrho,1,nodes);
free_vector(flvp,1,nodes);
free_vector(flvs,1,nodes);
free_vector(pts,1,L);
}



