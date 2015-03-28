/*------------------------------------------------------------------------
 *   read values of particle velocities (vx, vy) at the edges of the
 *   local grid from buffer arrays (which have been exchanged between
 *   processes)
 *   last update 12/02/02, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readbufv(float ** vx, float ** vy,  
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top) {

	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY;
	int i, j;



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=1;j<=NY;j++){

			vx[j][NX+1] = bufferlef_to_rig[j][1];
			vy[j][NX+1] = bufferlef_to_rig[j][2];
			vy[j][NX+2] = bufferlef_to_rig[j][3];
		
	}

						/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[1]!=0))	/* no boundary exchange at left edge of global grid */
	for (j=1;j<=NY;j++){
			vx[j][0] = bufferrig_to_lef[j][1];
			vy[j][0] = bufferrig_to_lef[j][2];
			vx[j][-1] = bufferrig_to_lef[j][3];
		
	}

	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=NX;i++){

			vx[NY+1][i] = buffertop_to_bot[i][1];
			vy[NY+1][i] = buffertop_to_bot[i][2];
			vx[NY+2][i] = buffertop_to_bot[i][3];
			
	}

	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=NX;i++){

			vx[0][i] = bufferbot_to_top[i][1];
			vy[0][i] = bufferbot_to_top[i][2];
			vy[-1][i] = bufferbot_to_top[i][3];
			
	}
}
