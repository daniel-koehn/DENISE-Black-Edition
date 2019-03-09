/*------------------------------------------------------------------------
 *   write seismograms to files 
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_glob(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID_SHOT, RUN_MULTIPLE_SHOTS, MODE;	
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
	extern char  SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];

        char vxf[STRING_SIZE], vyf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE];
        int nsrc=1;		
	
	 /*if (RUN_MULTIPLE_SHOTS){*/
         /*       sprintf(vxf,"%s.shot%d.%d",SEIS_FILE_VX,ishot,MYID_SHOT);
                sprintf(vyf,"%s.shot%d.%d",SEIS_FILE_VY,ishot,MYID_SHOT);
                sprintf(curlf,"%s.shot%d.%d",SEIS_FILE_CURL,ishot,MYID_SHOT);
                sprintf(divf,"%s.shot%d.%d",SEIS_FILE_DIV,ishot,MYID_SHOT);
                sprintf(pf,"%s.shot%d.%d",SEIS_FILE_P,ishot,MYID_SHOT);*/
        /*}
        else{
                sprintf(vxf,"%s.%d",SEIS_FILE_VX,MYID_SHOT);
                sprintf(vyf,"%s.%d",SEIS_FILE_VY,MYID_SHOT);
                sprintf(curlf,"%s.%d",SEIS_FILE_CURL,MYID_SHOT);
                sprintf(divf,"%s.%d",SEIS_FILE_DIV,MYID_SHOT);
                sprintf(pf,"%s.%d",SEIS_FILE_P,MYID_SHOT);
                 
        }*/

	
	if(MODE>0){
	  sprintf(vxf,"%s.shot%d.it%d",SEIS_FILE_VX,ishot,iter);
	  sprintf(vyf,"%s.shot%d.it%d",SEIS_FILE_VY,ishot,iter);
	  sprintf(pf,"%s.shot%d.it%d",SEIS_FILE_P,ishot,iter);
          sprintf(curlf,"%s.shot%d.it%d",SEIS_FILE_CURL,ishot,iter);
          sprintf(divf,"%s.shot%d.it%d",SEIS_FILE_DIV,ishot,iter);
	}
		
	if(MODE==0){
	  sprintf(vxf,"%s.shot%d",SEIS_FILE_VX,ishot);
	  sprintf(vyf,"%s.shot%d",SEIS_FILE_VY,ishot);
	  sprintf(pf,"%s.shot%d",SEIS_FILE_P,ishot);
          sprintf(curlf,"%s.shot%d",SEIS_FILE_CURL,ishot);
          sprintf(divf,"%s.shot%d",SEIS_FILE_DIV,ishot);
	}
	
	
	switch (SEISMO){
	case 1 : /* particle velocities only */
		fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",MYID_SHOT,ntr,vxf);
		outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID_SHOT,ntr,vyf);
		outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;
	case 2 : /* pressure only */
		fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",MYID_SHOT,ntr,pf);
		outseis_glob(fp,fopen(pf,"w"),1, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;
	case 3 : /* curl and div only */
		fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",MYID_SHOT,ntr,divf);
		outseis_glob(fp,fopen(divf,"w"), 0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",MYID_SHOT,ntr,curlf);
		outseis_glob(fp,fopen(curlf,"w"), 0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
	
		break;	
	case 4 : /* everything */
		fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",MYID_SHOT,ntr,vxf);
		outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID_SHOT,ntr,vyf);
		outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		
		fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",MYID_SHOT,ntr,pf);
		outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);

		/* Divergence and Curl output is not working */
		/*fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",MYID_SHOT,ntr,divf);
		outseis_glob(fp,fopen(divf,"w"),0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",MYID_SHOT,ntr,curlf);
		outseis(fp,fopen(curlf,"w"),0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);*/
		break;
		
      }     
}
