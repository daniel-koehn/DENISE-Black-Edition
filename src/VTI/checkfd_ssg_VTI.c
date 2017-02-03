/*-----------------------------------------------------------------------
 *  Check FD-Grid for stability and grid dispersion for the VTI problem.
 *  If the stability criterion is not fullfilled the program will
 *  terminate.
 *
 *  Daniel Koehn           
 *        
 *  last update  02/02/2017
 *  --------------------------------------------------------------------- */


#include "fd.h"

void checkfd_ssg_VTI(FILE *fp, float ** prho, float ** c11, float ** c13, float ** c33, float ** c44, float *hc){

	extern float DH, DT, TS;
        extern int NX, NY, MYID, INVMAT1, FW;

	/* local variables */

	float  c, fmax, gamma;
	float  cmax=0.0, cmin=1e9, dtstab, dhstab, cmax_r, cmin_r;
	int nfw=iround(FW/DH);
	int i, j, ny1=1, nx, ny, nx_min, ny_min;


	nx=NX; ny=NY; 

	/* low Q frame not yet applied as a absorbing boundary */
	/* if (!FREE_SURF) ny1=1+nfw;*/
	nfw=0;
	

	/* find maximum model phase velocity of shear waves at infinite
	      frequency within the whole model */
		for (i=1+nfw;i<=(nx-nfw);i++){
			for (j=ny1;j<=(ny-nfw);j++){
			        
				/* calculate seismic velocities in different spatial directions */
				c = sqrt(c11[j][i]/prho[j][i]);				

				if (cmax<c) cmax=c;
				if (cmin>c) cmin=c;

				c = sqrt(c13[j][i]/prho[j][i]);				

				if (cmax<c) cmax=c;
				if (cmin>c) cmin=c;

				c = sqrt(c33[j][i]/prho[j][i]);				

				if (cmax<c) cmax=c;
				if (cmin>c) cmin=c;

				c = sqrt(c44[j][i]/prho[j][i]);				

				if (cmax<c) cmax=c;
				if (cmin>c) cmin=c;

			}
		}


	if (MYID==0){
		fprintf(fp,"\n\n\n **Message from checkfd (printed by PE %d):\n",MYID);
		fprintf(fp," Minimum and maximum P-wave and S-wave velocities within subvolumes: \n ");		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/* find global maximum for Vp and global minimum for Vs*/
	MPI_Allreduce(&cmax,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&cmin,&cmin_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	cmax=cmax_r;
	cmin=cmin_r;	

        if(MYID==0){printf("cmin = %e m/s \t cmax = %e m/s \n", cmin, cmax);}

	fmax = 1.0/TS;
	dhstab = (cmin/(hc[0]*fmax));
	gamma = fabs(hc[1]) + fabs(hc[2]) + fabs(hc[3]) + fabs(hc[4]) + fabs(hc[5]) + fabs(hc[6]);
	dtstab = DH/(sqrt(2)*gamma*cmax);
	/*dtstab=DH/(sqrt(2.0)*cmax);*/

	/* find global minimum for NX and NY */
	MPI_Allreduce(&NX,&nx_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	MPI_Allreduce(&NY,&ny_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	

	if (MYID == 0) {

	fprintf(fp," Global values for entire model: \n");
	fprintf(fp," V_max= %e m/s \t V_min=%e m/s \n\n", cmax,cmin);
	fprintf(fp,"\n\n ------------------ CHECK FOR GRID DISPERSION --------------------\n");
	fprintf(fp," To satisfactorily limit grid dispersion the number of gridpoints \n");
	fprintf(fp," per minimum wavelength (of S-waves) should be 6 (better more).\n");
	fprintf(fp," Here the minimum wavelength is assumed to be minimum model phase velocity \n");
	fprintf(fp," (of S-waves) at maximum frequency of the source\n");
	fprintf(fp," devided by maximum frequency of the source.\n");
	fprintf(fp," Maximum frequency of the source is approximately %8.2f Hz\n",1.0/TS);
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

