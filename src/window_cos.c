/************************************************/
/* Cosine window Andre Kurzmann */
/* Update M. Schaefer 08.09.2011 */
/************************************************/

#include "fd.h"

void window_cos(float **win, int npad, int nsrc, float t1, float t2, float t3, float t4) {

	int	n, k, it1, it2, it3, it4, flank1, flank2;
	extern float DT;
	extern int NT;
	
	float ** dummy;
	
	dummy=matrix(1,nsrc,1,NT);
	
	for (k=1; k<=nsrc; k++) {
		for (n=1; n<=npad; n++) {
		
					dummy[k][n]=win[k][n];
					}
				}
	
	it1=1;
	it2=iround(t2/DT);
	it3=iround(t3/DT);
	it4=NT;
		
	flank1 = it2-it1;
	flank2 = it4-it3;
	
	for (k=1; k<=nsrc; k++) {
		for (n=1; n<=npad; n++) {
			if (n<=it1)			win[k][n] = dummy[k][n]*0.0;
			else if ((n>it1) && (n<it2))	win[k][n] = dummy[k][n]*((1-cos((float)(n-it1)/flank1*PI))/2.0);
			else if ((n>=it2) && (n<it3))	win[k][n] = dummy[k][n]*1.0;
			else if ((n>=it3) && (n<it4))	win[k][n] = dummy[k][n]*(1-(1-cos((float)(n-it3)/flank2*PI))/2.0);
			else if (n>=it4)		win[k][n] = dummy[k][n]*0.0;
					}
			}
	free_matrix(dummy,1,nsrc,1,NT);

}
