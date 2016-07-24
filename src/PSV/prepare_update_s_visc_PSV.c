/* Prepare_update_s_visc_PSV 
*
*  Daniel Koehn
*  Kiel, 24.07.2016 */

#include "fd.h"

void prepare_update_s_visc_PSV(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
		float **puipjp, float **ppi, float **prho, float **ptaus, float **ptaup,
		float **ptausipjp, float **f, float **g, float *bip, float *bjm,
		float *cip, float *cjm, float ***dip, float ***d, float ***e) {

	extern int NX, NY, L, INVMAT1;
	extern float DT, *FL;
	int i, j, l;
	extern char  MFILE[STRING_SIZE];
	extern float TS;
		
	float sumu, sumpi, ws, *pts;
	float mu, pi, muipjp;
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	}
	
	
	ws=2.0*PI*FL[1];
	/*ws=2.0*PI*(1.0/TS);*/
	
	sumu=0.0;
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
	}
	

	for (l=1;l<=L;l++){
		etajm[l] = peta[l];
		etaip[l] = peta[l];
	}
	
	if (INVMAT1==1){
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
				mu=(pu[j][i]*pu[j][i]*prho[j][i])/(1.0+sumu*ptaus[j][i]);
				pi=(ppi[j][i]*ppi[j][i]*prho[j][i])/(1.0+sumpi*ptaup[j][i]);
				muipjp=puipjp[j][i]/(1.0+sumu*ptausipjp[j][i]);
				
				fipjp[j][i] = muipjp*DT*(1.0+L*ptausipjp[j][i]);
				f[j][i] = mu*DT*(1.0+L*ptaus[j][i]);
				g[j][i] = pi*DT*(1.0+L*ptaup[j][i]);
				for (l=1;l<=L;l++){
					bip[l] = 1.0/(1.0+(etaip[l]*0.5));
					bjm[l] = 1.0/(1.0+(etajm[l]*0.5));
					cip[l] = 1.0-(etaip[l]*0.5);
					cjm[l] = 1.0-(etajm[l]*0.5);
					dip[j][i][l] = muipjp*etaip[l]*ptausipjp[j][i];
					d[j][i][l] = mu*etajm[l]*ptaus[j][i];
					e[j][i][l] = pi*etajm[l]*ptaup[j][i];
				}
			}
		}
	}
	
	if (INVMAT1==3){
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
				mu=pu[j][i]/(1.0+sumu*ptaus[j][i]);
				pi=(ppi[j][i]+2*pu[j][i])/(1.0+sumpi*ptaup[j][i]);
				muipjp=puipjp[j][i]/(1.0+sumu*ptausipjp[j][i]);
				
				fipjp[j][i] = muipjp*DT*(1.0+L*ptausipjp[j][i]);
				f[j][i] = mu*DT*(1.0+L*ptaus[j][i]);
				g[j][i] = pi*DT*(1.0+L*ptaup[j][i]);
				for (l=1;l<=L;l++){
					bip[l] = 1.0/(1.0+(etaip[l]*0.5));
					bjm[l] = 1.0/(1.0+(etajm[l]*0.5));
					cip[l] = 1.0-(etaip[l]*0.5);
					cjm[l] = 1.0-(etajm[l]*0.5);
					dip[j][i][l] = muipjp*etaip[l]*ptausipjp[j][i];
					d[j][i][l] = mu*etajm[l]*ptaus[j][i];
					e[j][i][l] = pi*etajm[l]*ptaup[j][i];
				}
			}
		}
	}	
	
	
	free_vector(pts,1,L);
}
