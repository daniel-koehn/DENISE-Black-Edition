/*------------------------------------------------------------------------
 *   Apply time damping (after Brossier (2009))                                 
 *   last update 31/08/11, D.Koehn
 *   modified    02/02/12, S.Heider
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void time_window(float **sectiondata, float *picked_times, int iter, int ntr_glob, int **recpos_loc, int ntr, int ns, int ishot){

/* declaration of variables */
extern float DT;
extern int REC1, REC2, MYID, TIMEWIN;
extern int POS[3];
extern float GAMMA, TWLENGTH_PLUS, TWLENGTH_MINUS, DT;
extern char PICKS_FILE[STRING_SIZE];
int READ_PICKED_TIMES;
char pickfile_char[STRING_SIZE];
float time, dump, dump1, taper, taper1;
float *pick_tmp, *pick_tmp1, *picked_times1;
int i, j, h;

FILE *fptime;

picked_times1 = vector(1,ntr);

/* read picked first arrival times from external pick file */
/* ------------------------------------------------------- */
if(TIMEWIN==1){

pick_tmp = vector(1,ntr_glob);
pick_tmp1 = vector(1,ntr_glob);

sprintf(pickfile_char,"%s%i.dat",PICKS_FILE,ishot);

fptime=fopen(pickfile_char,"r");
if (fptime == NULL) {
err(" picks_?.dat could not be opened !");
}

  for(i=1;i<=ntr_glob;i++){
    fscanf(fptime,"%f%f",&dump, &dump1);
    pick_tmp[i] = dump + TWLENGTH_PLUS;
    pick_tmp1[i] = dump1;
  }

fclose(fptime);

/* distribute picks on CPUs */
h=1;
  for(i=1;i<=ntr;i++){

    picked_times[h] = pick_tmp[recpos_loc[3][i]];
    picked_times1[h] = pick_tmp1[recpos_loc[3][i]];
    
    h++;}

free_vector(pick_tmp,1,ntr_glob);
free_vector(pick_tmp1,1,ntr_glob);
    
} /* end of if(TIMEWIN==1) */

/* read picked first arrival times from STA/LTA-picker file */
/* ------------------------------------------------------- */
if(TIMEWIN==4){

sprintf(pickfile_char,"%sshot%d.%d%d",PICKS_FILE,ishot,POS[1],POS[2]);

fptime=fopen(pickfile_char,"r");
if (fptime == NULL) {
err(" picks_?.dat could not be opened !");
}

  for(i=1;i<=ntr;i++){
    fscanf(fptime,"%e%e",&dump, &dump1);
    picked_times[i] = dump + TWLENGTH_PLUS;
    picked_times1[i] = dump1;
  }

fclose(fptime);
    
} /* end of if(TIMEWIN==4) */


/* Define constant time windows */
/* ---------------------------- */
if(TIMEWIN==2){
  h=1;
    for(i=1;i<=ntr;i++){
    
        picked_times[h] = -DT;
        picked_times1[h] = TWLENGTH_PLUS;
                
        h++;
    }                    
}

/* calculate RMS */
for(i=1;i<=ntr;i++){
      for(j=1;j<=ns;j++){
      
         if((TIMEWIN==1)||(TIMEWIN==2)){time = (float)(j * DT);}
         if((TIMEWIN==3)||(TIMEWIN==4)){time = (float)((ns-j+1) * DT);}

         dump = (time-picked_times[i]);
         taper = exp(-GAMMA*dump*dump);
         
         dump1 = (time-picked_times1[i]); 
         taper1 = exp(-TWLENGTH_MINUS*dump*dump);
	 
	 if(TIMEWIN==1){

           if(time>=picked_times[i]){
             sectiondata[i][j] = sectiondata[i][j] * taper;
           }          

           /*if(time<picked_times[i]){
             sectiondata[i][j] = sectiondata[i][j] * taper1;
           }*/    
                          
	 }
	 
	 if(TIMEWIN==2){
           if(time>=picked_times1[i]){
             sectiondata[i][j] = sectiondata[i][j] * taper;
           }  

           /*if(time<=picked_times[i]){
             sectiondata[i][j] = sectiondata[i][j] * taper1;
           }*/
         }
         
         if((TIMEWIN==3)||(TIMEWIN==4)){
          
            if(time>=picked_times[i]){
              sectiondata[i][j] = sectiondata[i][j] * taper;
            }  
         
         }  
     	   
      }     
}

free_vector(picked_times1,1,ntr);

} /* end of function time_window.c */
