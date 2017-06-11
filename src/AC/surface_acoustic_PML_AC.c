/*------------------------------------------------------------------------
 *   stress free surface condition
 *   
 *  D. Koehn,
 *  Kiel, 10.06.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_acoustic_PML_AC(int ndepth, float ** p){


	int i,j,m;
	int fdoh;
	extern int NX, FDORDER;
	
	fdoh = FDORDER/2;

	j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
        for (i=1;i<=NX;i++){
        p[j][i] = 0;
		for (m=1; m<=fdoh; m++) {
			p[j-m][i] = -p[j+m][i];
		}
	}

}
