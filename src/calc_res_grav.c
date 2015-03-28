/* Calculate gravity data residuals */
/*------------------------------------------------------------------------
 *
 * Daniel Koehn
 * Kiel, the 9th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double calc_res_grav(int ngrav, float *gz_mod, float *gz_res){

	int i, j, it;
        float * gz_data, l2, tmp;
        extern char GRAV_DATA_IN[STRING_SIZE];
        extern char GRAV_DATA_OUT[STRING_SIZE];
        extern int MYID;
	FILE *FP;

        gz_data = vector(1,ngrav);

        /* load gravity field data */
        FP=fopen(GRAV_DATA_IN,"r");         

        for(i=1;i<=ngrav;i++){
           fscanf(FP,"%d \t %e \n",&it, &tmp);
	   gz_data[i]=tmp;
        }

        fclose(FP);
        
        /* calculate gravity data residuals */
        l2=0.0;
        for(i=1;i<=ngrav;i++){
           
           gz_res[i] = gz_mod[i] - gz_data[i];

           if(MYID==0){
             printf("%d \t %e \n",i,gz_res[i]);
           }

           l2 += gz_res[i] * gz_res[i];

        }

       /* output of the gravity data residuals */
       if(MYID==0){

         FP=fopen(GRAV_DATA_OUT,"w");         

         for(i=1;i<=ngrav;i++){
         	fprintf(FP,"%d \t %e \n",i,gz_res[i]);
         }

         fclose(FP);

       }

       /* clean memory */
       free_vector(gz_data,1,ngrav);

       return l2;
}
