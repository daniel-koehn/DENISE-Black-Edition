/*
 *   Model with slow formation and fluid-filled borehole
 *   last update 11.12.07, O. Hellwig
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DH;
	extern int   NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float rhov, muv, piv, vp, vs, y0, t, de, zplat1, zplat2, rplat;
	float *pts, ts, tp, sumu, sumpi, ws, *ri, *d, *ti, *dl, *vsl, z, r;
	float **checkp, **checks, **checkrho; 
	int   i, j, l, ii, jj, nk, k, nl;
	char filename_mu[STRING_SIZE];
	char filename_rho[STRING_SIZE]; 
	char filename_pi[STRING_SIZE];
				
	sprintf(filename_mu,"%s.mu",MFILE);
	sprintf(filename_rho,"%s.rho",MFILE);
	sprintf(filename_pi,"%s.pi",MFILE);
	
	/* parameters for borehole */  
        const float vp1=6100.0, vs1=3500.0, rho1=7850.0; /* steel */
        	
	/* parameters for layer 1 */
	const float vp2=1375.0, vs2=250.0, rho2=1800.0;
	
	/* parameters for air */
        const float vpa=300.0, vsa=1e-6, rhoa=1.25;
	
	/* parameters for water */    
        const float vpw=1500.0, vsw=1e-6, rhow=1000.0;         
	
	/* borehole radius */
	float R, a , b, dw; 
	const float rshift=10.0;
	
	/* parameters for the checkquerboard CTS test */
	int nyc, nxc;
	int idcx, idcy;
	float cxsign, cysign;
	
	nyc = 50;
	nxc = 50;
	
	cxsign = 1.0;
	cysign = 1.0;
	
	nk = 24;  /* number of knots for the Fino Monopile */
	nl = 7;  /* number of elastic layers (interface depth) */
	y0 = 22.0; /*  thickness of the air layer */
	dw = 21.8 + y0; /* depth of the sea floor in model coordinates */ 
	de = 0.01; /* small increment to add for definition of first point in the new layer */
	
	/* dimensions of the FINO platform */
	zplat1 = -(21.50-y0);
	zplat2 = -(20.5-y0);
	rplat  = 6.5;
	
	ri=vector(1,nk);
	d=vector(1,nk);
	dl=vector(1,nl);
	ti=vector(1,nk-1);
	vsl=vector(1,nl);
	
	checkp = matrix(1,nyc,1,nxc);
	checks = matrix(1,nyc,1,nxc);
	checkrho = matrix(1,nyc,1,nxc);
	
	dl[1] = dw;
	dl[2] = 6.0+dl[1];
	dl[3] = de+dl[2];
	dl[4] = 15.0+dl[2];
	dl[5] = de+dl[4];
	dl[6] = 30.0+dl[4];
        dl[7] = de+dl[6];
        
        vsl[1] = 150.0;
        vsl[2] = 220.0;
        vsl[3] = 250.0;
        vsl[4] = 310.0;
        vsl[5] = 350.0;
        vsl[6] = 380.0;
        vsl[7] = 430.0;

        /* Monopile model */
	
	ri[1] = 2.650;
	ri[2] = 2.650;
	ri[3] = 2.850;
	ri[4] = 3.050;
	ri[5] = 3.074;
	ri[6] = 3.146;
	ri[7] = 3.25;
	ri[8] = 3.0;
	ri[9] = 3.0;
	ri[10] = 3.0;
	ri[11] = 3.0;
	ri[12] = 3.0;
	ri[13] = 3.0;
	ri[14] = 3.078;
	ri[15] = 3.5;
	ri[16] = 3.747;
	ri[17] = 4.219;
	ri[18] = 4.455;
	ri[19] = 4.75;
	ri[20] = 4.75;
	ri[21] = 4.75;
	ri[22] = 4.75;
	ri[23] = 4.75;
	ri[24] = 4.75;
	
	d[1] = 21.5;
	d[2] = 20.0;
	d[3] = 15.4;
	d[4] = 11.956;
	d[5] = 11.0;
	d[6] = 8.2;
	d[7] = 4.2;	
	d[8] = 2.7;
	d[9] = 2.2;
	d[10] = 0.2;
	d[11] = -0.3;
	d[12] = -2.3;
	d[13] = -3.25;
	d[14] = -4.5;
	d[15] = -11.2;
	d[16] = -13.3;
	d[17] = -17.3;
	d[18] = -19.3;
	d[19] = -21.8;
	d[20] = -25.7;
	d[21] = -39.2;
	d[22] = -41.0;
	d[23] = -42.8;
	d[24] = -51.8;
	
	ti[1] = 50;
	ti[2] = 50;
	ti[3] = 50;
	ti[4] = 50;
	ti[5] = 50;
	ti[6] = 50;
	ti[7] = 50;
	
	ti[8] = 50;
	ti[9] = 50;
	ti[10] = 50;
	ti[11] = 50;
	ti[12] = 50;
	ti[13] = 50;
	ti[14] = 50;
	ti[15] = 50;
	ti[16] = 50;
	ti[17] = 50;
	ti[18] = 50;
	ti[19] = 50;
	ti[20] = 50;
	ti[21] = 50;
	ti[22] = 50;
	ti[23] = 50;
	
	/* shift model on the grid */
	for (k=1;k<=nk;k++){
	   d[k] = -(d[k]-y0);
	   ri[k] = ri[k]/2.0;
	   ti[k] = ti[k]/1000.0;
	   }
	
	/* define parameters for CTS checkerboard test */
        for (i=1;i<=nxc;i++){
	    for (j=1;j<=nyc;j++){
	       checkp[j][i]=100.0 ;
	       checks[j][i]=100.0 ;
	       checkrho[j][i]=0.0;
	    }
	}	
	
	/*-----------------------------------------------------------------------*/

	/* loop over global grid */
	idcx=1;
	for (i=1;i<=NXG;i++){
           idcy=1;
                cysign=-cysign;
		for (j=1;j<=NYG;j++){
	                
			  vs = 3821.0;
			  vp = 2206.0;
			  rhov = 2481.0;
	                
	                  /*if((i>nxc)&&(i<NXG-nxc)&&(j<NYG-nyc)){
	                  vs = 1333.0 + cxsign * cysign * checks[idcy][idcx];
	                  vp = 2310.0 + cxsign * cysign * checkp[idcy][idcx];
			  rhov = 1900.0 + cxsign * cysign * checkrho[idcy][idcx];
			  }
			  idcy++;*/
	                
	                /* calculate shape of the FINO Monopile */
	                /*for (k=1;k<=nk-1;k++){
	                
	                    if((z>=d[k])&&(z<=d[k+1])){
	                    
	                       if(ri[k]!=ri[k+1]){
	                       
	                       a = (d[k]-d[k+1])/(ri[k]-ri[k+1]);
	                       b = d[k] - a * ri[k];
	                       R = (z-b)/a;
	                       
	                       }
	                       
	                        if(ri[k]==ri[k+1]){
	                          R = ri[k];
	                        }
	                        
	                        t = ti[k];
	                    }
	                    
	                }*/ 

			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii = i-POS[1]*NX;
				jj = j-POS[2]*NY;
				
				u[jj][ii]    = vs;
				rho[jj][ii]  = rhov;
				pi[jj][ii]   = vp;
			}
			
		        if(idcy>nyc){
		           idcy = 1;
		           cysign = - cysign;
		        }	
		}
	   idcx++;	
	   if(idcx>nxc){
	      idcx = 1;
	      cxsign = - cxsign;
	   }	
	}	

	/* each PE writes his model to disk */
	writemod(filename_rho,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_rho,3);
	
	/* each PE writes his model to disk */
	writemod(filename_pi,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_pi,3);
	                        
	/* each PE writes his model to disk */
	writemod(filename_mu,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_mu,3);
	                        

	/*free_vector(pts,1,L);
	free_vector(ri,1,nk);
	free_vector(d,1,nk);
	free_vector(ti,1,nk-1);*/
}


