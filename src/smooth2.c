/*------------------------------------------------------------------------
 *   Apply smooth2 filter to gradient
 *   
 *   Daniel Koehn
 *   Kiel, the 4th of January 2014 
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void  smooth2(float ** grad){

	/* declaration of extern variables */
        extern int NX, NY, NXG, NYG, IDX, IDY, GRADT2;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern int SPAT_FILT_SIZE;
	extern float WD_DAMP, WD_DAMP1;
	extern char JACOBIAN[STRING_SIZE];
	
	/* declaration of local variables */
	int i,j, h, fa, fb, itr, ix, iz, n1,n2;
	int tracl1, jj, ii, nmax, * win;
	float gradtmp, ** v, r1, r2, **w, *d, *e, *f, rw;
        
	char jac[STRING_SIZE];
	FILE *fp_grad, *FP3;
	
        /* define parameters */
        r1=WD_DAMP;
        r2=WD_DAMP1;

        n1=NYG;
        n2=NXG;

        /* scale the smoothing parameter */
        r1 = r1*r1*0.25;
        r2 = r2*r2*0.25;
        
        /* allocate space */
        nmax = (n1<n2)?n2:n1;
        win = ivector(0,4);
        w = matrix(0,n2,0,n1);
        d = vector(0,nmax);
        e = vector(0,nmax);
        f = vector(0,nmax);

        /* define windows function */
        win[0] = 0;
        win[1] = n1;
        win[2] = 0;
        win[3] = n2;
        rw=0.;
        rw=rw*rw*0.25;
    
	                    
	/* temporarily save gradient for wavenumber filtering */
        sprintf(jac,"%s_wavenumber.old.%i.%i",JACOBIAN,POS[1],POS[2]);
        FP3=fopen(jac,"wb");
                        
        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                fwrite(&grad[j][i],sizeof(float),1,FP3);
            }                       
        }                           
                                    
        fclose(FP3);                
                                    
        MPI_Barrier(MPI_COMM_WORLD);
                                  
        /* merge gradient file */ 
        sprintf(jac,"%s_wavenumber.old",JACOBIAN);
        if (MYID==0) mergemod(jac,3);  


if(MYID==0){    /* read the global model on node 0 and apply wavenumber damping */
        
        /* define temporary gradient matrix */
        v = matrix(0,n2,0,n1);
        
	printf("\n Smooth2 is applied to gradient (written by PE %d)\n",MYID); 
	
	fp_grad=fopen(jac,"rb");
	
	if (fp_grad==NULL) err(" Could not open gradient file ! ");
	
	/* load merged gradient */
	for (i=1;i<=n2;i++){
	   for (j=1;j<=n1;j++){
	        
	            fread(&gradtmp, sizeof(float), 1, fp_grad);
	            v[i][j] = gradtmp;		

            }
	}
	
	fclose(fp_grad);
	
	                                                                                                                                                                                             
	/* define the window function */
        for(ix=0; ix<n2; ++ix)
                for(iz=0; iz<n1; ++iz)
                        w[ix][iz] = 0;  
        for(ix=win[2]; ix<win[3]; ++ix)
                for(iz=win[0]; iz<win[1]; ++iz)
                        w[ix][iz] = 1;  

        if(win[0]>0 || win[1]<n1 || win[2]>0 || win[3]<n2){
        /* smooth the window function */
                for(iz=0; iz<n1; ++iz){
                        for(ix=0; ix<n2; ++ix){
                                d[ix] = 1.0+2.0*rw;
                                e[ix] = -rw;
                                f[ix] = w[ix][iz];
                        }
                        d[0] -= rw;
                        d[n2-1] -= rw;
                        tripd(d,e,f,n2);
                        for(ix=0; ix<n2; ++ix)
                                w[ix][iz] = f[ix];
                }
                for(ix=0; ix<n2; ++ix){
                        for(iz=0; iz<n1; ++iz){
                                d[iz] = 1.0+2.0*rw;
                                e[iz] = -rw;
                                f[iz] = w[ix][iz];
                }
                        d[0] -= rw;
                        d[n1-1] -= rw;
                        tripd(d,e,f,n1);
                        for(iz=0; iz<n1; ++iz)
                                w[ix][iz] = f[iz];
                }
        }

        /* solving for the smoothing velocity */
        for(iz=0; iz<n1; ++iz){
                for(ix=0; ix<n2-1; ++ix){
                        d[ix] = 1.0+r2*(w[ix][iz]+w[ix+1][iz]);
                        e[ix] = -r2*w[ix+1][iz];
                        f[ix] = v[ix][iz];
                }
                d[0] -= r2*w[0][iz];
                d[n2-1] = 1.0+r2*w[n2-1][iz];
                f[n2-1] = v[n2-1][iz];
                tripd(d,e,f,n2);
                for(ix=0; ix<n2; ++ix)
                        v[ix][iz] = f[ix];
        }
         for(ix=0; ix<n2; ++ix){
                for(iz=0; iz<n1-2; ++iz){
                        d[iz] = 1.0+r1*(w[ix][iz+1]+w[ix][iz+2]);
                        e[iz] = -r1*w[ix][iz+2];
                        f[iz] = v[ix][iz+1];
                }
                f[0] += r1*w[ix][1]*v[ix][0];
                d[n1-2] = 1.0+r1*w[ix][n1-1];
                f[n1-2] = v[ix][n1-1];
                tripd(d,e,f,n1-1);
                for(iz=0; iz<n1-1; ++iz)
                        v[ix][iz+1] = f[iz];
        }

      
	/* write damped gradient to temporary file */
	sprintf(jac,"%s_smooth2.old",JACOBIAN);
	FP3=fopen(jac,"wb");

	for (i=1;i<=n2;i++){
		for (j=1;j<=n1;j++){
		
                    gradtmp = v[i][j];
	            fwrite(&gradtmp,sizeof(float),1,FP3);	 
			 
		}
	}
	fclose(FP3);
	        
	/* free memory */
        free_matrix(v,0,n2,0,n1);
        
} /* end of if MYID==0*/

	 MPI_Barrier(MPI_COMM_WORLD);

         sprintf(jac,"%s_smooth2.old",JACOBIAN);
	 FP3=fopen(jac,"rb");

	 /* distribute spatial filtered gradient on computational nodes */
	 for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){
			
			fread(&gradtmp, sizeof(float),1,FP3);

			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				grad[jj][ii]=gradtmp;

			}
			
		}
	}

        fclose(FP3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(jac,"%s_wavenumber.old.%i.%i",JACOBIAN,POS[1],POS[2]);
        remove(jac);
                                

}
