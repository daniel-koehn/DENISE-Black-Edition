/*------------------------------------------------------------------------
 *  Check if parameters MODE and PHYSICS are reasonable
 *
 *  D. Koehn
 *  Kiel, 02.02.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void check_mode_phys(){

	/* declaration of extern variables */
        extern int MODE, PHYSICS, MYID;
	
        if (MYID==0){

		printf("\n **Message from check_model_phys (printed by PE %d):\n",MYID);
		printf("\n");
		printf(" -----------------------  DENISE operation mode  ----------------------\n");		
		switch (MODE){
		case 0 :
			printf(" MODE=%d: Only forward modeling is applied.\n",MODE);
			break;
		case 1 :
			printf(" MODE=%d: Time-domain FWI is applied.\n",MODE);
			break;
		case 2 :
			printf(" MODE=%d: Elastic Reverse Time Migration RTM is applied.\n",MODE);
			break;
		case 3 :
			printf(" MODE=%d: Calculate FD-based gradient.\n",MODE);
			break;

		default:
			err(" Sorry, DENISE operation MODE unkown ! ");
		}
		
		printf("\n\n");

		printf(" -----------------------  DENISE Physics  ----------------------\n");		

		switch (PHYSICS){
		case 1 :
			printf(" PHYSICS=%d: Solve 2D isotropic elastic PSV problem.\n",PHYSICS);
			break;
		case 2 :
			printf(" PHYSICS=%d: Solve 2D acoustic problem.\n",PHYSICS);
			break;
		case 3 :
			printf(" PHYSICS=%d: Solve 2D elastic PSV VTI problem.\n",PHYSICS);
			break;
		case 4 :
			printf(" PHYSICS=%d: Solve 2D elastic PSV TTI problem.\n",PHYSICS);
			break;
		case 5 :
			printf(" PHYSICS=%d: Solve 2D isotropic elastic SH problem.\n",PHYSICS);
			break;

		default:
			err(" Ups, you are obviously in a parallel universe: PHYSICS unkown - call 07700 900461 ! ");
		}

		printf("\n\n");

        }
}
