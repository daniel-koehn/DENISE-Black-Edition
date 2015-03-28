/*------------------------------------------------------------------------
 *   write seismograms to files 
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID, RUN_MULTIPLE_SHOTS;	
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
	extern char  SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];

        char vxf[STRING_SIZE], vyf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE];
        int nsrc=1;		
	
	 /*if (RUN_MULTIPLE_SHOTS){*/
                sprintf(vxf,"%s.shot%d.%d",SEIS_FILE_VX,ishot,MYID);
                sprintf(vyf,"%s.shot%d.%d",SEIS_FILE_VY,ishot,MYID);
                sprintf(curlf,"%s.shot%d.%d",SEIS_FILE_CURL,ishot,MYID);
                sprintf(divf,"%s.shot%d.%d",SEIS_FILE_DIV,ishot,MYID);
                sprintf(pf,"%s.shot%d.%d",SEIS_FILE_P,ishot,MYID);
        /*}
        else{
                sprintf(vxf,"%s.%d",SEIS_FILE_VX,MYID);
                sprintf(vyf,"%s.%d",SEIS_FILE_VY,MYID);
                sprintf(curlf,"%s.%d",SEIS_FILE_CURL,MYID);
                sprintf(divf,"%s.%d",SEIS_FILE_DIV,MYID);
                sprintf(pf,"%s.%d",SEIS_FILE_P,MYID);
                 
        }*/

	
	if(iter>0){
	sprintf(vxf,"%s.shot%d_it%d.%d",SEIS_FILE_VX,ishot,iter,MYID);
	sprintf(vyf,"%s.shot%d_it%d.%d",SEIS_FILE_VY,ishot,iter,MYID);
	sprintf(pf,"%s.shot%d_it%d.%d",SEIS_FILE_P,ishot,iter,MYID);
	}
	
	if(iter==-1){
	sprintf(vxf,"%s.shot%d_it%d.%d",SEIS_FILE_VX,ishot,iter,MYID);
	sprintf(vyf,"%s.shot%d_it%d.%d",SEIS_FILE_VY,ishot,iter,MYID);
	sprintf(pf,"%s.shot%d_it%d.%d",SEIS_FILE_P,ishot,iter,MYID);
	}
	
	
	switch (SEISMO){
	case 1 : /* particle velocities only */
		fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",MYID,ntr,vxf);
		outseis(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;
	case 2 : /* pressure only */
		fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",MYID,ntr,SEIS_FILE_P);
		outseis(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;
	case 3 : /* curl and div only */
		fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",MYID,ntr,SEIS_FILE_DIV);
		outseis(fp,fopen(divf,"w"), 0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",MYID,ntr,SEIS_FILE_CURL);
		outseis(fp,fopen(curlf,"w"), 0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
	
		break;	
	case 4 : /* everything */
		fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",MYID,ntr,SEIS_FILE_VX);
		outseis(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID,ntr,SEIS_FILE_VY);
		outseis(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		
		fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",MYID,ntr,SEIS_FILE_P);
		outseis(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);

		fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",MYID,ntr,SEIS_FILE_DIV);
		outseis(fp,fopen(divf,"w"),0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",MYID,ntr,SEIS_FILE_CURL);
		outseis(fp,fopen(curlf,"w"),0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);	
		break;
		
      }     
}
