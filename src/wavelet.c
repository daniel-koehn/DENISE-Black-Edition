/*------------------------------------------------------------------------
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*   last update 19/01/02, T. Bohlen
*  ----------------------------------------------------------------------*/

#include "fd.h"
#include <complex.h>

float ** wavelet(float ** srcpos_loc, int nsrc, int ishot){


	/* extern variables */
	extern int QUELLART, NT, MYID;
	extern float  DT, TS, FC_SPIKE_1, FC_SPIKE_2;
	extern char SIGNAL_FILE[STRING_SIZE];
	extern FILE *FP;

	/*local variables */
	int nts, nt, k;
	float *psource=NULL, tshift, amp=0.0, a, fc, tau, t, ts, ag;
	float ** signals, k1, f0, t1;
	double complex amp1;
        char signal_wave[STRING_SIZE];

        sprintf(signal_wave,"%s_shot_%i.dat",SIGNAL_FILE,ishot);
	if (QUELLART==3) psource=rd_sour(&nts,fopen(signal_wave,"r"));
	
	signals=fmatrix(1,nsrc,1,NT);
	
	if(QUELLART==7){
	  k1 = (FC_SPIKE_2-FC_SPIKE_1)/TS;    /* rate of change of frequency with time */
	  f0 = (FC_SPIKE_2+FC_SPIKE_1)/2.0;   /* midfrequency of bandwidth */
	}
	
	for (nt=1;nt<=NT;nt++){
			t=(float)nt*DT;
			
			for (k=1;k<=nsrc;k++) {
				tshift=srcpos_loc[4][k];
				fc=srcpos_loc[5][k];
				a=srcpos_loc[6][k];
				ts=1.0/fc;
                                ag = PI*PI*fc*fc;
                                
				switch (QUELLART){
					case 1 : 
						tau=PI*(t-1.5*ts-tshift)/(1.5*ts); /* Ricker */
						amp=(((1.0-4.0*tau*tau)*exp(-2.0*tau*tau)));
					break;
					case 2 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=((sin(2.0*PI*(t-tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t-tshift)*fc)));

/*						amp=((sin(2.0*PI*(t+tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t+tshift)*fc)));
*/
					break;
					case 3 : 
						amp=psource[nt];
					break;  /* source wavelet from file SOURCE_FILE */
					case 4 : 
						/*tau=PI*(t-ts-tshift)/(1.5*ts);*/ /* Ricker */
						/*amp=((t-ts-tshift)*exp(-2.0*tau*tau));*/
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=pow(sin(PI*(t+tshift)/ts),3.0);
						break; /* sinus raised to the power of three */
					break; 
					case 5 : 
				                /* first derivative of a Gaussian */
				                ts=1.2/fc;
					        ag  = PI*PI*fc*fc;
			                        amp = - 2.0 * ag * (t-ts) * exp(-ag*(t-ts)*(t-ts));
					break;					
					case 6 : 
					        /* Bandlimited Spike */
						amp=0.0;
						if(nt==1+iround(tshift/DT)){
						amp = 1.0;}	
					break;
					case 7 : 
						/* Klauder wavelet */
						ts = 1.5/FC_SPIKE_1;
						amp=0.0;
						t1 = t-ts-tshift;
						amp=creal(cexp(0.0+2.0*PI*f0*t1*I)*sin(PI*k1*t1*(TS-t1))/(PI*k1*t1));
					break;                                                                                                                                              	
					default : 
						err("Which source-wavelet ? ");
					}
					
					
					signals[k][nt]=amp*a;
		}
	}
	
	fprintf(FP," Message from function wavelet written by PE %d \n",MYID);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
	fprintf(FP," have been assigned with a source signal. \n");
			
		
	if (QUELLART==3) free_vector(psource,1,NT);
	
	return signals;	

}
