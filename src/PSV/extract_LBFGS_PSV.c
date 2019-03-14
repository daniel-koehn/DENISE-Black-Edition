/*-------------------------------------------------------------------------------------
 * Extract gradients from Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS)
 * optimization (Nocedal and Wright, 2006) and store old gradients/models for the
 * PSV problem
 * 
 * Daniel Koehn
 * Kiel, 25.07.2016
 *
 * ------------------------------------------------------------------------------------ */

#include "fd.h"

void extract_LBFGS_PSV( int iter, float ** waveconv, float ** gradp, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float **ppi, float ** pu, float ** prho, float * r_LBFGS){

	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID, MYID_SHOT;
        extern float C_vp, C_vs, C_rho;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, h;
        int ki, itershift;
	FILE *FP3;
	
        itershift = 1;


/* extract updated Hessian-gradient product from r_LBFGS vector */
/* ------------------------------------------------------------ */

if(iter>1){

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
	sprintf(jac,"%s_p_vp.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
	sprintf(jac,"%s_p.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
        sprintf(jac,"%s_c.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
	sprintf(jac,"%s_p_vs.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
	sprintf(jac,"%s_p_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
        sprintf(jac,"%s_c_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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

	sprintf(jac,"%s_p_mrho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
	sprintf(jac,"%s_p_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
        sprintf(jac,"%s_c_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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
