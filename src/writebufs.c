/*------------------------------------------------------------------------
 *   write values of stress tensor components at the edges of the
 *   local grid into buffer arrays (which will be exchanged between
 *   processes)
 *   last update 12/02/02, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void writebufs(float ** sxx, float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top){


	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY;
	int i, j;



	if ((BOUNDARY) || (POS[1]!=0))	/* no boundary exchange at left edge of global grid */
	for (j=1;j<=NY;j++){
			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][1] =  sxx[j][1];
			bufferlef_to_rig[j][2] =  sxy[j][1];
			bufferlef_to_rig[j][3] =  sxx[j][2];
	}


	if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=1;j<=NY;j++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][1] =  sxx[j][NX];
			bufferrig_to_lef[j][2] =  sxy[j][NX];
			bufferrig_to_lef[j][3] =  sxy[j][NX-1];

	}

	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=NX;i++){
			/* storage of top of local volume into buffer */
			buffertop_to_bot[i][1]  = syy[1][i];
			buffertop_to_bot[i][2]  = sxy[1][i];
			buffertop_to_bot[i][3]  = syy[2][i];
	}


	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=NX;i++){			
			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][1]  = syy[NY][i];
			bufferbot_to_top[i][2]  = sxy[NY][i];
			bufferbot_to_top[i][3]  = sxy[NY-1][i];
			
	}


}
