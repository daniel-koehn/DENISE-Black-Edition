/*-------------------------------------------------------------------------------------
 * Extract gradients from Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS)
 * optimization (Nocedal and Wright, 2006) and store old gradients/models for the
 * SH problem
 * 
 * Daniel Koehn
 * Kiel, 14.10.2018
 *
 * ------------------------------------------------------------------------------------ */

#include "fd.h"

void extract_LBFGS_SH( int iter, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float ** waveconv_ts, float ** gradp_ts, float ** pu, float ** prho,  float ** ptaus, float * r_LBFGS){

	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID;
        extern float C_vs, C_rho, C_taus, C_vs_min, C_rho_min, C_taus_min;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, h;
        int ki, itershift;
	float tmp;
	FILE *FP3;
	
        itershift = 1;


/* extract updated Hessian-gradient product from r_LBFGS vector */
/* ------------------------------------------------------------ */

if(iter>1000){

     /* update gradients */
     h=1;
     
     /* Vs */
     for (i=1;i<=NX;i=i+IDX){   
        for (j=1;j<=NY;j=j+IDY){
                                 
	    waveconv_u[j][i] = r_LBFGS[h];
	    h++;  
                                                                  
	}
     }
                                                                               
     /* Density */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            waveconv_rho[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }

     /* Taus */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            waveconv_ts[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }

}
	
/* save old models Vs */
/* ------------------ */

    /* save old model */
	sprintf(jac,"%s_p_vs.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	       tmp = (pu[j][i] - C_vs_min) / C_vs;
               fwrite(&tmp,sizeof(float),1,FP3);
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
/* ------------------- */

	sprintf(jac,"%s_p_mrho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	       tmp = (prho[j][i] - C_rho_min) / C_rho;
               fwrite(&tmp,sizeof(float),1,FP3);
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
	
/* save old models Taus */
/* -------------------- */

    /* save old model */
	sprintf(jac,"%s_p_mts.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	       tmp = (ptaus[j][i] - C_taus_min) / C_taus;
               fwrite(&tmp,sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_mts.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_ts.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_ts[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_ts.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c_ts.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_ts[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_ts.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);	
	
}
