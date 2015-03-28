/* Output of the Hessian */
/*                                                 */
/* Daniel Koehn */
/* Kiel, the 10th of November 2013 */

#include "fd.h"

void hessian_out(float ** hessian_lam, float ** hessian_mu, float ** hessian_rho, float ** ppi, float ** pu, float ** prho){
  
extern int NX, NY, IDX, IDY, DTINV, INVMAT1, MYID, POS[4], FDORDER, SPATFILTER;
extern char JACOBIAN[STRING_SIZE];

/* local variables */
int i, j, k, l, ns_hess, ishot, irec, nd, NSRC_HESSIAN, NREC_HESSIAN, RECINC;
char jac[STRING_SIZE];

float lamss, muss;
float ** hessian_vp, ** hessian_vs, ** hessian_rhos;

FILE *FP4;

nd = FDORDER/2 + 1;

/* Diagonal elements of the Hessian*/
hessian_vp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_vs = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_rhos = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
           hessian_vp[j][i]=0.0;
           hessian_vs[j][i]=0.0;
           hessian_rhos[j][i]=0.0;
    }
}

/* construct Hessian for different material parameters */
for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
                                
                    
        if(INVMAT1==1){
          muss = prho[j][i] * pu[j][i] * pu[j][i];
          lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;}
                    
        if(INVMAT1==3){
          muss = pu[j][i];
          lamss = ppi[j][i];}
                      
        /* new Pseudo-Hessian with improved correction of the geometrical spreading */ 
        hessian_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * hessian_lam[j][i];           
        hessian_vp[j][i] = 2.0 * ppi[j][i] * prho[j][i] * hessian_lam[j][i];
        hessian_vs[j][i] = (- 4.0 * prho[j][i] * pu[j][i] * hessian_lam[j][i]) + 2.0 * prho[j][i] * pu[j][i] * hessian_mu[j][i];
        hessian_rhos[j][i]= ((((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * hessian_lam[j][i]) + (pu[j][i] * pu[j][i] *hessian_mu[j][i]) + hessian_rho[j][i]);      
                  
    }
}

/* apply wavenumber damping for Vp-, Vs- and density Hessian */
/*if(SPATFILTER==1){
  wavenumber(hessian);
  wavenumber(hessian_u);
  wavenumber(hessian_rho);
}*/

/* save Hessian for Vp */
/* ----------------------- */
sprintf(jac,"%s_hessian.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
        fwrite(&hessian_vp[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */   
sprintf(jac,"%s_hessian",JACOBIAN);
if (MYID==0) mergemod(jac,3);

/* save HESSIAN for mu */
/* ----------------------- */
sprintf(jac,"%s_hessian_u.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
        fwrite(&hessian_vs[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */   
sprintf(jac,"%s_hessian_u",JACOBIAN);
if (MYID==0) mergemod(jac,3);

/* save HESSIAN for rho */   
/* ----------------------- */
sprintf(jac,"%s_hessian_rho.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
       fwrite(&hessian_rhos[j][i],sizeof(float),1,FP4);
   }
}
    
fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */
sprintf(jac,"%s_hessian_rho",JACOBIAN);
if (MYID==0) mergemod(jac,3);

free_matrix(hessian_vp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_vs,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_rhos,-nd+1,NY+nd,-nd+1,NX+nd);

}
