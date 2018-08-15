/* Prepare_update_s_visc_SH
*
*  Daniel Koehn
*  Kiel, 03.12.2017 
*/

#include "fd.h"

void prepare_update_s_visc_SH(float *etajm, float *etaip, float *peta, float **fipjp, float **pujp, 
		float **puip, float **prho, float **ptaus, float **ptausipjp, float **f, float **g, 
		float *bip, float *bjm, float *cip, float *cjm, float ***dip, float ***d, float ***e) {

	extern int NX, NY, L, INVMAT1;
	extern float DT, *FL, ALPHA_VISC;
	int i, j, l;
	extern char  MFILE[STRING_SIZE];
	extern float TS;
		
	float sumu, ws, *pts;
	float mujp, muip;
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	}
	
	
	ws=2.0*PI*FL[1];
	
	sumu=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
	}
	
	ALPHA_VISC = sumu;

	for (l=1;l<=L;l++){
		etajm[l] = peta[l];
		etaip[l] = peta[l];
	}
	
	if (INVMAT1==1){
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){

				mujp=pujp[j][i]/(1.0+sumu*ptaus[j][i]);
				muip=puip[j][i]/(1.0+sumu*ptausipjp[j][i]);				
				
				fipjp[j][i] = muip*DT*(1.0+L*ptausipjp[j][i]);
				    f[j][i] = mujp*DT*(1.0+L*ptaus[j][i]);
			
				for (l=1;l<=L;l++){
					bip[l] = 1.0/(1.0+(etaip[l]*0.5));
					bjm[l] = 1.0/(1.0+(etajm[l]*0.5));
					cip[l] = 1.0-(etaip[l]*0.5);
					cjm[l] = 1.0-(etajm[l]*0.5);
					dip[j][i][l] = muip*etaip[l]*ptausipjp[j][i];
					d[j][i][l] = mujp*etajm[l]*ptaus[j][i];					
				}
			}
		}
	}
	
	if (INVMAT1==3){
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){

				mujp = pujp[j][i]/(1.0+sumu*ptaus[j][i]);
				muip = puip[j][i]/(1.0+sumu*ptausipjp[j][i]);
				
				fipjp[j][i] = muip*DT*(1.0+L*ptausipjp[j][i]);
				f[j][i] = mujp*DT*(1.0+L*ptaus[j][i]);

				for (l=1;l<=L;l++){
					bip[l] = 1.0/(1.0+(etaip[l]*0.5));
					bjm[l] = 1.0/(1.0+(etajm[l]*0.5));
					cip[l] = 1.0-(etaip[l]*0.5);
					cjm[l] = 1.0-(etajm[l]*0.5);
					dip[j][i][l] = muip*etaip[l]*ptausipjp[j][i];
					d[j][i][l] = mujp*etajm[l]*ptaus[j][i];

				}
			}
		}
	}	
	
	
	free_vector(pts,1,L);
}
