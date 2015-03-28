/*-------------------------------------------------------------
 *  Check FD-Grid for stability and grid dispersion.
 *  If the stability criterion is not fullfilled the program will
 *  terminate.                   
 *  last update  26/03/2012 
 *
 *  ----------------------------------------------------------*/


#include "fd.h"

void checkfd_ssg_visc(FILE *fp, float ** prho, float ** ppi, float ** pu, float ** ptaus, float ** ptaup, float *peta, float *hc){

	extern float DH, DT, TS;
        extern int NX, NY, MYID, INVMAT1, FW, L;

	/* local variables */

	float  c, cmax_p=0.0, cmin_p=1e9, cmax_s=0.0, cmin_s=1e9, fmax, gamma;
	float  cmax=0.0, cmin=1e9, dtstab, dhstab, cmax_r, cmin_r;
	float sumu, sumpi, ws, ts, ppi_ref, pu_ref;
	int nfw=iround(FW/DH);
	int i, j, l, ny1=1, nx, ny, nx_min, ny_min;
	
	
	ws=2.0*PI/TS; /*center frequency of source*/
	
	nx=NX; ny=NY; 

	/* low Q frame not yet applied as a absorbing boundary */
	/* if (!FREE_SURF) ny1=1+nfw;*/
	nfw=0;
	
		

	/* find maximum model phase velocity of shear waves at infinite
	      frequency within the whole model */
		for (i=1+nfw;i<=(nx-nfw);i++){
			for (j=ny1;j<=(ny-nfw);j++){
			
				sumu=0.0;
				for (l=1;l<=L;l++){
					ts=DT/peta[l];
					sumu=sumu+((ws*ws*ts*ts*ptaus[j][i])/(1.0+ws*ws*ts*ts));
				}
				
				if(INVMAT1==1){
					pu_ref=prho[j][i]*pu[j][i]*pu[j][i];}
				if(INVMAT1==3){
					pu_ref=pu[j][i];}
					
				
				/* minimum phase velocity of shear waves */
				c=sqrt(pu_ref/(prho[j][i]*(1.0+sumu)));
								
				if (cmin_s>c) cmin_s=c;
				
				
				/* maximum phase velocity of shear waves */
				c=sqrt(pu_ref*(1.0+L*ptaus[j][i])/(prho[j][i]*(1.0+sumu)));
				
				if (cmax_s<c) cmax_s=c;
			}
		}
	
	
	
	/* find maximum model phase velocity of P-waves at infinite
		 frequency within the whole model */
		for (i=1+nfw;i<=(nx-nfw);i++){
			for (j=ny1;j<=(ny-nfw);j++){
			        
				
				sumpi=0.0;
				for (l=1;l<=L;l++){
					ts=DT/peta[l];
					sumpi=sumpi+((ws*ws*ts*ts*ptaup[j][i])/(1.0+ws*ws*ts*ts));
				}
				
				if(INVMAT1==1){
					ppi_ref=prho[j][i]*ppi[j][i]*ppi[j][i];}
				if(INVMAT1==3){
					ppi_ref=ppi[j][i]+2*pu[j][i];}
					
				/* minimum phase velocity of P waves */
				c=sqrt(ppi_ref/(prho[j][i]*(1.0+sumpi)));
				
				if (cmin_p>c) cmin_p=c;
				
				
				/* maximum phase velocity of shear waves */
				c=sqrt(ppi_ref*(1.0+L*ptaup[j][i])/(prho[j][i]*(1.0+sumpi)));
				
				if (cmax_p<c) cmax_p=c;
			}
		}


	if (MYID==0){
		fprintf(fp,"\n\n\n **Message from checkfd (printed by PE %d):\n",MYID);
		fprintf(fp," Minimum and maximum P-wave and S-wave velocities within subvolumes: \n ");
		fprintf(fp," MYID\t Vp_min(f=0) \t Vp_max(f=inf) \t Vs_min(f=0) \t Vsmax(f=inf) \n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(fp," %d \t %e \t %e \t %e \t %e \n", MYID, cmin_p, cmax_p, cmin_s, cmax_s);

	if (cmax_s>cmax_p) cmax=cmax_s; 
	else cmax=cmax_p;
	if (cmin_s<cmin_p) cmin=cmin_s; 
	else cmin=cmin_p;

	/* find global maximum for Vp and global minimum for Vs*/
	MPI_Allreduce(&cmax,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&cmin,&cmin_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	cmax=cmax_r;
	cmin=cmin_r;	

	fmax=2.0/TS;
	dhstab = (cmin/(hc[0]*fmax));
	gamma = fabs(hc[1]) + fabs(hc[2]) + fabs(hc[3]) + fabs(hc[4]) + fabs(hc[5]) + fabs(hc[6]);
	dtstab = DH/(sqrt(2)*gamma*cmax);
	/*dtstab=DH/(sqrt(2.0)*cmax);*/

	/* find global minimum for NX and NY */
	MPI_Allreduce(&NX,&nx_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	MPI_Allreduce(&NY,&ny_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	

	if (MYID == 0) {

	fprintf(fp," Global values for entire model: \n");
	fprintf(fp," Vp_max= %e m/s \t Vs_min=%e m/s \n\n", cmax,cmin);
	fprintf(fp,"\n\n ------------------ CHECK FOR GRID DISPERSION --------------------\n");
	fprintf(fp," To satisfactorily limit grid dispersion the number of gridpoints \n");
	fprintf(fp," per minimum wavelength (of S-waves) should be 6 (better more).\n");
	fprintf(fp," Here the minimum wavelength is assumed to be minimum model phase velocity \n");
	fprintf(fp," (of S-waves) at maximum frequency of the source\n");
	fprintf(fp," devided by maximum frequency of the source.\n");
	fprintf(fp," Maximum frequency of the source is approximately %8.2f Hz\n",2.0/TS);
	fprintf(fp," The minimum wavelength (of S-waves) in the following simulation will\n");
	fprintf(fp," be %e meter.\n", cmin/fmax);
	fprintf(fp," Thus, the recommended value for DH is %e meter.\n", dhstab);
	fprintf(fp," You have specified DH= %e meter.\n\n", DH);
	if (DH>dhstab)
		warning(" Grid dispersion will influence wave propagation, choose smaller grid spacing (DH).");
	
	fprintf(fp," \n\n ----------------------- CHECK FOR STABILITY ---------------------\n");
	fprintf(fp," The following simulation is stable provided that\n\n");
	fprintf(fp," \t p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n\n");
	fprintf(fp," where cmax is the maximum phase velocity at infinite frequency\n");
	fprintf(fp," and gamma = sum(|FD coeff.|)\n");

	fprintf(fp," In the current simulation cmax is %8.2f m/s .\n\n",cmax);

	fprintf(fp," DT is the timestep and DH is the grid size.\n\n");
	fprintf(fp," In this simulation the stability limit for timestep DT is %e seconds .\n",dtstab);
	fprintf(fp," You have specified DT= %e s.\n", DT);
	if (DT>dtstab)
		err(" The simulation will get unstable, choose smaller DT. ");
	else fprintf(fp," The simulation will be stable.\n");

	fprintf(fp,"\n\n ----------------------- ABSORBING BOUNDARY ------------------------\n");
        if((FW>nx_min)||(FW>ny_min)){
	  err(" The width of the absorbing boundary is larger than one computational domain. Choose smaller FW or use less CPUs.");
	}

	fprintf(fp," Width (FW) of absorbing frame should be at least 10 gridpoints.\n");
	fprintf(fp," You have specified a width of %d gridpoints.\n",FW);
	if (FW<10) 
		warning(" Be aware of artificial reflections from grid boundaries ! \n");

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}

