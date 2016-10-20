/*------------------------------------------------------------------------
 *  Write DENISE parameters                           
 *
 *  D. Koehn
 *  Kiel, 02.08.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){

	/* declaration of extern variables */
	extern int   NX, NY, NT, QUELLART, FDORDER, MAXRELERROR;
	extern int  SNAP, SNAP_FORMAT, SNAP_SHOT, L, SRCREC, TAPER;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, FPML, npower, k_max_PML;
	extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, FW;
	extern int  READMOD, READREC, BOUNDARY;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
	extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE];
	extern char SIGNAL_FILE[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern char  MFILE[STRING_SIZE], JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE];
	extern int NP, NPROCX, NPROCY, MYID;
	
	extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT1, MODE, PHYSICS, QUELLTYPB;
	extern int GRAD_METHOD, ORDER_SPIKE, RTM_SHOT;
	extern float FC_SPIKE_1, FC_SPIKE_2;
	extern int FILT_SIZE, MODEL_FILTER;
	extern int FILT_SIZE_GRAD, GRAD_FILTER;
	
	extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
	extern int SWS_TAPER_FILE;
	extern float SRTRADIUS, EXP_TAPER_GRAD_HOR;
	extern int MIN_ITER;
	extern char INV_MODELFILE[STRING_SIZE];
	extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
	
	extern int DTINV, RTMOD;
	extern int STEPMAX, TRKILL;
	extern float EPS_SCALE, SCALEFAC;
	extern char  TRKILL_FILE[STRING_SIZE];

	extern int NORMALIZE, NLBFGS, N_STREAMER;
        extern float REC_INCR_X, REC_INCR_Y;
	
	extern char MISFIT_LOG_FILE[STRING_SIZE];

	extern int GRAVITY, NGRAVB, NZGRAV, GRAV_TYPE;
	extern char DFILE[STRING_SIZE];
	extern int BACK_DENSITY;
	
	
	/* definition of local variables */
	int l;
	

	fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
	fprintf(fp,"\n");
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",NPROCX);
	fprintf(fp," Number of PEs in vertical direction (NPROCY): %d\n",NPROCY);
	fprintf(fp," Total number of PEs in use: %d\n",NP);
	fprintf(fp,"\n");
	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
	fprintf(fp," Grid-spacing (DH): %e meter\n", DH);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
	fprintf(fp," Timestep (DT): %e seconds\n", DT);
	fprintf(fp," Number of timesteps: %i \n",NT);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- FD ORDER -----------------------------\n");
	fprintf(fp," FDORDER = %d\n",FDORDER);
	fprintf(fp," MAXRELERROR = %d\n",MAXRELERROR);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");

	if (SRCREC){
		fprintf(fp," reading source positions, time delay, centre frequency \n");
		fprintf(fp," and initial amplitude from ASCII-file \n");
		fprintf(fp,"\t%s\n\n",SOURCE_FILE);}

	fprintf(fp," wavelet of source:");

	switch (QUELLART){
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," reading from \n\t %s\n",SIGNAL_FILE);
		break;
	case 4 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		break;
	case 5 :
		fprintf(fp," 1st derivative of Gaussian \n");
		break;
	case 6 :
                if(FC_SPIKE_1 > 0.0){
		  fprintf(fp," Bandpass filtered spike \n");
                  fprintf(fp," FC_SPIKE_1 = %f, FC_SPIKE_2 = %f, ORDER_SPIKE = %d \n",FC_SPIKE_1,FC_SPIKE_2,ORDER_SPIKE);
                }

                if(FC_SPIKE_1 < 0.0){
		  fprintf(fp," Lowpass filtered spike \n");
                  fprintf(fp," FC_SPIKE_2 = %f, ORDER_SPIKE = %d \n",FC_SPIKE_2,ORDER_SPIKE);
                }
		break;
	case 7 :
	        fprintf(fp," Klauder\n");
	        fprintf(fp," FC_SPIKE_1 = %f, FC_SPIKE_2 = %f, TS = %f \n",FC_SPIKE_1,FC_SPIKE_2,TS);
	        break;			                                	
	default :
		err(" Sorry, incorrect specification of source wavelet ! ");
	}

        fprintf(fp,"\n\n");

	/*fprintf(fp," Type of source:");
	switch (QUELLTYP){
	case 1 :
		fprintf(fp," explosive source \n");
		break;
	case 2 :
		fprintf(fp," point source with directive force in x-direction\n");
		break;
	case 3 :
		fprintf(fp," point source with directive force in (vertical) y-direction\n");
		break;
	case 4 :
		fprintf(fp," rotated point source with directive force in x- and y-direction\n");
		break;
        case 5 :
                fprintf(fp," source defined by moment tensor\n");
                break;	                                 	
	default :
		err(" Sorry, wrong source type specification ! ");
	}
	
	fprintf(fp,"\n");*/

	if (SEISMO){
		fprintf(fp," ------------------------- RECEIVER  --------------------------\n");
		if (READREC){
                        if(READREC==1){
			    fprintf(fp," reading receiver positions from single file \n");
			    fprintf(fp,"\t%s.dat\n\n",REC_FILE);
			}
                        if(READREC==2){
			    fprintf(fp," reading receiver positions from multiple files \n");
			    fprintf(fp,"\t%s_shot_X.dat\n\n",REC_FILE);
			}
			fprintf(fp," reference_point_for_receiver_coordinate_system:\n");
			fprintf(fp," x=%f \ty=%f\t z=%f\n",REFREC[1], REFREC[2], REFREC[3]);
		}

                if (N_STREAMER){
                        fprintf(fp," ------------------------- Towed streamer  --------------------------\n");
                        fprintf(fp," Assuming the first N_STREAMER = %d receivers belong to a streamer \n",N_STREAMER);
                        fprintf(fp," Shifting the streamer by ... \n");
                        fprintf(fp," REC_INCR_X=%f in x-direction \n",REC_INCR_X);
                        fprintf(fp," REC_INCR_Y=%f in y-direction \n",REC_INCR_Y);
                }
	}

	fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
	if (FREE_SURF) fprintf(fp," free surface at the top of the model ! \n");
	else fprintf(fp," no free surface at the top of the model ! \n");
	fprintf(fp,"\n");

	fprintf(fp," ------------------------- CPML ---------------------\n");
	if (FW>0.0){
		fprintf(fp," width of absorbing frame is %i gridpoints.\n",FW);
		fprintf(fp," CPML damping applied. \n");
		fprintf(fp," Damping velocity in the PML frame in m/s: %f .\n",DAMPING);
		fprintf(fp," Frequency within the PML frame in Hz: %f \n",FPML);
		fprintf(fp," npower: %f \n",npower);
		fprintf(fp," k_max: %f \n",k_max_PML); 
	}
	else fprintf(fp," absorbing frame not installed ! \n");


	switch (BOUNDARY){
	case 0 :
		fprintf(fp," No periodic boundary condition.\n");
		break;
	case 1 :
		fprintf(fp," Periodic boundary condition at left and right edges.\n");
		break;
	default :
		warning(" Wrong integer value for BOUNDARY specified ! ");
		warning(" No periodic boundary condition will be applied ");
		BOUNDARY=0;
		break;
	}

	if (READMOD){
		fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
		fprintf(fp," names of model-files: \n");
		fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
		fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);
		fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
		fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
		fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);
		for (l=1;l<=L;l++) fprintf(fp,"\t %1i. relaxation frequencies: %s.f%1i\n",l,MFILE,l);
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
	fprintf(fp," Number of relaxation mechanisms (L): %i\n",L);
	fprintf(fp," The L relaxation frequencies are at:  \n");
	for (l=1;l<=L;l++) fprintf(fp,"\t%f",FL[l]);
	fprintf(fp," Hz\n");
	fprintf(fp," Value for tau is : %f\n",TAU);


	if (SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(SNAP){
		case 1:
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			break;
		case 2:
			fprintf(fp," pressure field.\n");
			break;
		case 3:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			break;
		case 4:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			fprintf(fp," x- and y-component of particle velocity.\n");
			break;
		default:
			err(" sorry, incorrect value for SNAP ! \n");
		}

		fprintf(fp," \t write snapshots for shot SNAP_SHOT= %d \n", SNAP_SHOT);
		fprintf(fp," \t first (TSNAP1)= %8.5f s\n", TSNAP1);
		fprintf(fp," \t last (TSNAP2)=%8.5f s\n",TSNAP2);
		fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
		fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,NX);
		fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,NY);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
		switch (SNAP_FORMAT){
		case 1 :
			err(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			err(" Don't know the format for the Snapshot-data ! \n");
		}
	
		fprintf(fp,"\n\n");
	}
	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		if ((SEISMO==1) || (SEISMO==4)){
			fprintf(fp," seismograms of ");
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE_VX,SEIS_FILE_VY);
		}
		if ((SEISMO==2) || (SEISMO==4)){
			fprintf(fp," seismograms of pressure field (hydrophones).\n");
			fprintf(fp," output-file: \n ");
			fprintf(fp,"\t%s\n",SEIS_FILE_P);
		}
		if ((SEISMO==3) || (SEISMO==4)){
			fprintf(fp," seismograms of curl and div.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE_CURL,SEIS_FILE_DIV);
			
		}		
	
		switch (SEIS_FORMAT){
		case 1 :
			fprintf(fp," The data is written in IEEE SU-format . \n");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary IEEE (4 byte per float)");
			break;
		default:
			err(" Sorry. I don't know the format for the seismic data ! \n");
		}
		fprintf(fp," samplingrate of seismic data: %f s\n",NDT*DT);
		fprintf(fp," Number of samples per trace: %i \n", iround(NT/NDT));
		fprintf(fp," ----------------------------------------------------------\n");
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
        fprintf(fp," -----------------------  DENISE elastic specific parameters  ----------------------\n");
	fprintf(fp,"\n Maximum number of iterations: %d\n",ITERMAX);
	fprintf(fp," location of the measured seismograms : \n ");
	fprintf(fp,"\t%s\n\n",DATA_DIR);
	if (INVMAT1==1){
		fprintf(fp," INVMAT1=%d: Inversion parameters are vp, vs and rho.\n",INVMAT1);}
	if (INVMAT1==2){
		fprintf(fp," INVMAT1=%d: Inversion parameters are Zp, Zs and rho.\n",INVMAT1);}
	if (INVMAT1==3){
		fprintf(fp," INVMAT1=%d: Inversion parameters are lambda, mu and rho.\n",INVMAT1);}
	if (QUELLTYPB==1){
		fprintf(fp," QUELLTYPB=%d: Inversion of x and y component.\n\n",QUELLTYPB);}
	if (QUELLTYPB==2){
		fprintf(fp," QUELLTYPB=%d: Inversion of y component.\n\n",QUELLTYPB);}
	if (QUELLTYPB==3){
		fprintf(fp," QUELLTYPB=%d: Inversion of x component.\n\n",QUELLTYPB);}
        if (QUELLTYPB==4){
                fprintf(fp," QUELLTYPB=%d: Inversion of p component.\n\n",QUELLTYPB);}
        if (QUELLTYPB==5){
                fprintf(fp," QUELLTYPB=%d: Inversion of x and p component.\n\n",QUELLTYPB);}
        if (QUELLTYPB==6){
                fprintf(fp," QUELLTYPB=%d: Inversion of y and p component.\n\n",QUELLTYPB);}
        if (QUELLTYPB==7){
                fprintf(fp," QUELLTYPB=%d: Inversion of x, y and p component.\n\n",QUELLTYPB);}
	
	fprintf(fp," Shots used for step length estimation:\n");
	fprintf(fp,"\t TESTSHOT_START = %d \n",TESTSHOT_START);
	fprintf(fp,"\t TESTSHOT_END = %d \n",TESTSHOT_END);
	fprintf(fp,"\t TESTSHOT_INCR = %d \n\n",TESTSHOT_INCR);
	
	fprintf(fp," Cosine Taper used : \n ");
	fprintf(fp,"\t%d\n",TAPER); 
	
	fprintf(fp," Log file for misfit in each iteration step: \n ");
	fprintf(fp,"\t%s \n\n",MISFIT_LOG_FILE); 
	
	fprintf(fp," Output of inverted models to: \n ");
	fprintf(fp,"\t%s \n\n",INV_MODELFILE);
		
	fprintf(fp," --------------- Gradient tapering -------------------\n");
	if (SWS_TAPER_GRAD_VERT==1){
		fprintf(fp," SWS_TAPER_GRAD_VERT=%d: Vertical taper applied.\n",SWS_TAPER_GRAD_VERT);
		fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d)\n\n",GRADT1,GRADT2,GRADT3,GRADT4);} 
	else	fprintf(fp," SWS_TAPER_GRAD_VERT=%d: No vertical taper applied.\n\n",SWS_TAPER_GRAD_VERT);
	
	if (SWS_TAPER_GRAD_HOR==1){
		fprintf(fp," SWS_TAPER_GRAD_HOR=%d: Horizontal taper applied.\n",SWS_TAPER_GRAD_HOR);
		fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d, EXP_TAPER_GRAD_HOR=%f)\n\n",GRADT1,GRADT2,GRADT3,GRADT4,EXP_TAPER_GRAD_HOR);} 
	else	fprintf(fp," SWS_TAPER_GRAD_HOR=%d: No horizontal taper applied.\n\n",SWS_TAPER_GRAD_HOR);
	
	if (SWS_TAPER_GRAD_SOURCES==1){
		fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: Taper around the sources.\n",SWS_TAPER_GRAD_SOURCES);
		fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);} 
	else	fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: No taper around the sources applied.\n\n",SWS_TAPER_GRAD_SOURCES);
	
	if (SWS_TAPER_CIRCULAR_PER_SHOT==1){
		fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: Taper around the source for each shot.\n",SWS_TAPER_CIRCULAR_PER_SHOT);
		fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);} 
	else	fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: No taper around the sources applied.\n\n",SWS_TAPER_CIRCULAR_PER_SHOT);
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Gradient smoothing with 2D-Gaussian filter -------------------\n");
	if(MODEL_FILTER==1){
		fprintf(fp," GRAD_FILTER=%d: Jacobians are filtered.(FILT_SIZE_GRAD=%d)\n",MODEL_FILTER,FILT_SIZE_GRAD);}
	else 	fprintf(fp," MODEL_FILTER=%d: Jacobians are not filtered.\n",MODEL_FILTER);
		
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Limits of model parameters -------------------\n");
	fprintf(fp," VPLOWERLIM = %f \t\t VPUPPERLIM = %f \n",VPLOWERLIM,VPUPPERLIM);
	fprintf(fp," VSLOWERLIM = %f \t\t VSUPPERLIM = %f \n",VSLOWERLIM,VSUPPERLIM);
	fprintf(fp," RHOLOWERLIM = %f \t RHOUPPERLIM = %f \n",RHOLOWERLIM,RHOUPPERLIM);
	
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Optimization method -------------------\n");
	switch(GRAD_METHOD){
		case 1:
			fprintf(fp," GRAD_METHOD=%d: PCG\n",GRAD_METHOD);
			break;
		case 2:
			fprintf(fp," GRAD_METHOD=%d: LBFGS\n",GRAD_METHOD);
                        fprintf(fp," NLBFGS=%d \n",NLBFGS);
			break;
		case 3:
			fprintf(fp," GRAD_METHOD=%d: \n",GRAD_METHOD);
			break;
		default:
			err(" Sorry, incorrect value for GRAD_METHOD ! \n");
	}
	
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Model smoothing -------------------\n");
	if(MODEL_FILTER==1){
		fprintf(fp," MODEL_FILTER=%d: vp and vs models are filtered after each iteration step.\n",MODEL_FILTER);
		fprintf(fp," (FILT_SIZE=%d)\n",FILT_SIZE);}
	else 	fprintf(fp," MODEL_FILTER=%d: vp and vs models are not filtered after each iteration step.\n",MODEL_FILTER);
	

	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Trace kill -------------------\n");
	if (TRKILL){
		fprintf(fp," TRKILL=%d: Trace kill is applied \n",TRKILL);
		fprintf(fp," Reading trace kill matrix from file: %s \n\n",TRKILL_FILE);}
	else fprintf(fp," TRKILL=%d: No trace kill is applied \n",TRKILL);
	
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Trace normalization -------------------\n");
	if (NORMALIZE){
		fprintf(fp," NORMALIZE=%d: The measured and synthetic seismograms will be normalized.\n",NORMALIZE);
		fprintf(fp," before calculating the residuals. \n\n");}
	else fprintf(fp," NORMALIZE=%d: No normalization of measured and synthetic seismograms.\n",NORMALIZE);
	
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Reduce size of inversion grid -------------------\n");
	fprintf(fp," Every %d time sample is used for the calculation of the gradients.\n\n",DTINV);
	
	
	fprintf(fp,"\n\n");
	fprintf(fp," --------------- Step length estimation -------------------\n");
	fprintf(fp," EPS_SCALE = %f\n",EPS_SCALE);
	fprintf(fp," STEPMAX = %d\n",STEPMAX);
	fprintf(fp," SCALEFAC = %f\n",SCALEFAC);

        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Reverse Time Modelling -------------------\n");
        if (RTMOD==0){
           fprintf(fp," RTMOD=%d: No Reverse Time Modelling applied.\n",RTMOD);
        }
        if (RTMOD==1){
           fprintf(fp," RTMOD=%d: Reverse Time Modelling will be applied.\n",RTMOD);
        }
		
	fprintf(fp,"\n\n");  
	fprintf(fp," --------------- Gravity Modelling and Inversion -------------------\n");
	if (GRAVITY==0){
	   fprintf(fp," No Gravity Modelling and Inversion applied. \t GRAVITY=%d \n",GRAVITY);
	} 
	if (GRAVITY==1){
	   fprintf(fp," Gravity Modelling will be applied. \t\t\t GRAVITY=%d \n",GRAVITY);
	}
	if (GRAVITY==2){
	   fprintf(fp," Gravity Modelling and Inversion will be applied.\t GRAVITY=%d \n",GRAVITY);
	}
	
	fprintf(fp," Boundary in x-direction [gridpoints] \t\t\t NGRAVB = %d \n",NGRAVB);
	fprintf(fp," Boundary in z-direction [m] \t\t\t\t NZGRAV = %d \n",NZGRAV);
	
	if (GRAV_TYPE==1){
	   fprintf(fp," Modelling and Inversion of Gravity Data. \t\t GRAV_TYPE=%d \n",GRAV_TYPE);
	}
	if (GRAV_TYPE==2){
	   fprintf(fp," Modelling and Inversion of Gravity Gradient Data. \t\t GRAV_TYPE=%d \n",GRAV_TYPE);
	}
	
	if (BACK_DENSITY==1){
	   fprintf(fp," Initial Model is used as Background Density. \t\t BACK_DENSITY=%d \n",BACK_DENSITY);
	}
	if (BACK_DENSITY==2){
	   fprintf(fp," Self-defined Model is used as Background Density. \t BACK_DENSITY=%d \n",BACK_DENSITY);
	   fprintf(fp,"\t background density file:\n\t %s \n",DFILE);
	}

	                                                        
	fprintf(fp,"\n");
	fprintf(fp," **************************************************************\n");
	fprintf(fp,"\n");
}
