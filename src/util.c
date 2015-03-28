#define NR_END 1
#define FREE_ARG char *

#include "fd.h"


void err(char err_text[]){
        extern int MYID;

        fprintf(stdout,"Message from PE %d\n",MYID);
        fprintf(stdout,"R U N - T I M E  E R R O R: \n");
        fprintf(stdout,"%s\n",err_text);
        fprintf(stdout,"...now exiting to system.\n");

        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
}

void warning(char warn_text[]){
	/* standard warnings handler */
	fprintf(stdout,"W A R N I N G   M E S S A G E: \n");
	fprintf(stdout,"%s\n",warn_text);
}


double maximum(float **a, int nx, int ny){
	/* find absolute maximum of array a[1...nx][1...ny] */
	double maxi=0.0;
	int i, j;


	for (j=1;j<=ny;j++)
		for (i=1;i<=nx;i++)
			if (fabs((double) a[i][j])>maxi) maxi=fabs((double)a[i][j]);
	return maxi;
}

float minimum_m(float **mat, int nx, int ny)
{
        float minm;
        int i,j;

        minm = mat[1][1];
        for (i=1;i<=nx;i++)
                for (j=1;j<=ny;j++)
                {
                        if (i*j==1) continue;
                        if (mat[j][i] < minm)
                        {
                                minm = mat[j][i];
                        }
                }
        return minm;
}

float maximum_m(float **mat, int nx, int ny)
{
        float maxm;
        int i,j;

        maxm = mat[1][1];
        for (i=1;i<=nx;i++)
                for (j=1;j<=ny;j++)
                {
                        if (i*j==1) continue;
                        if (mat[j][i] > maxm)
                        {
                                maxm = mat[j][i];
                        }
                }
        return maxm;
}




float *vector(int nl, int nh){
	/* allocate a float vector with subscript range v[nl..nh] and initializing
		   this vector, eg. vector[nl..nh]=0.0 */
	float *v;
	int i;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) err("allocation failure in function vector()");
	for (i=0;i<(nh-nl+1+NR_END);i++) v[i]=0.0;
	return v-nl+NR_END;
}


int *ivector(int nl, int nh){
	/* allocate an int vector with subscript range v[nl..nh] and initializing
		   this vector, eg. ivector[nl..nh]=0 */
	int *v;
	int i;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) err("allocation failure in function ivector()");
	for (i=0;i<(nh-nl+1+NR_END);i++) v[i]=0;
	return v-nl+NR_END;
}

unsigned short int *usvector(int nl, int nh){
	/* allocate an short int vector with subscript range v[nl..nh] and initializing
		   this vector, eg. ivector[nl..nh]=0 */
	unsigned short int *v;
	int i;

	v=(unsigned short int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned short int)));
	if (!v) err("allocation failure in function usvector()");
	for (i=0;i<(nh-nl+1+NR_END);i++) v[i]=0;
	return v-nl+NR_END;
}


unsigned char *cvector(int nl, int nh){
	/* allocate an unsigned char vector with subscript range v[nl..nh] */
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) err("allocation failure in function cvector()");
	return v-nl+NR_END;
}


unsigned long *lvector(int nl, int nh){
	/* allocate an unsigned long vector with subscript range v[nl..nh] and
		  initializing this vector, eg. vector[nl..nh]=0.0 */
	unsigned long *v;
	int i;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned long)));
	if (!v) err("allocation failure in function lvector()");
	for (i=0;i<(nh-nl+1+NR_END);i++) v[i]=0;
	return v-nl+NR_END;
}

double *dvector(int nl, int nh){
	/* allocate a double vector with subscript range v[nl..nh] and initializing
		   this vector, eg. vector[nl..nh]=0.0 */
	double *v;
	int i;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) err("allocation failure in function dvector()");
	for (i=0;i<(nh-nl+1+NR_END);i++) v[i]=0.0;
	return v-nl+NR_END;
}

float **fmatrix(int nrl, int nrh, int ncl, int nch){
	/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) err("allocation failure 1 in function fmatrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) err("allocation failure 2 in function fmatrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* initializing matrix */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++) m[i][j]=0.0;

	/* return pointer to array of pointer to rows */
	return m;
}

float **matrix(int nrl, int nrh, int ncl, int nch){
	/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) err("allocation failure 1 in function matrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) err("allocation failure 2 in function matrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* initializing matrix */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++) m[i][j]=0.0;

	/* return pointer to array of pointer to rows */
	return m;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch){
	/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) err("allocation failure 1 in function matrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) err("allocation failure 2 in function dmatrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* initializing matrix */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++) m[i][j]=0.0;

	/* return pointer to array of pointer to rows */
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch){
	/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t) ((nrow+NR_END)*sizeof(int*)));
	if (!m) err("allocation failure 1 in function imatrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) err("allocation failure 2 in function imatrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* initializing matrix */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++) m[i][j]=0;

	/* return pointer to array of pointer to rows */
	return m;
}

unsigned short int **usmatrix(int nrl, int nrh, int ncl, int nch){
	/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	unsigned short int **m;

	/* allocate pointers to rows */
	m=(unsigned short int **) malloc((size_t) ((nrow+NR_END)*sizeof(unsigned short int*)));
	if (!m) err("allocation failure 1 in function usmatrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl]=(unsigned short int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned short int)));
	if (!m[nrl]) err("allocation failure 2 in function usmatrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* initializing matrix */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++) m[i][j]=0;

	/* return pointer to array of pointer to rows */
	return m;
}

float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh){
	/* allocate a float 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh]=0.0 */
	int i,j,d, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t) ((nrow+NR_END)*sizeof(float**)));
	if (!t) err("allocation failure 1 in function f3tensor() ");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) err("allocation failure 2 in function f3tensor() ");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) err("allocation failure 3 in function f3tensor() ");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for (i=nrl+1;i<=nrh;i++){
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* initializing 3tensor */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++)
			for (d=ndl;d<=ndh;d++) t[i][j][d]=0.0;

	/* return pointer to array of pointer to rows */
	return t;
}


int ***i3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh){
	/* allocate a integer 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh]=0.0 */
	int i,j,d, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;

	/* allocate pointers to pointers to rows */
	t=(int ***) malloc((size_t) ((nrow+NR_END)*sizeof(int**)));
	if (!t) err("allocation failure 1 in function i3tensor() ");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
	if (!t[nrl]) err("allocation failure 2 in function i3tensor() ");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
	if (!t[nrl][ncl]) err("allocation failure 3 in function i3tensor() ");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for (i=nrl+1;i<=nrh;i++){
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* initializing 3tensor */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++)
			for (d=ndl;d<=ndh;d++) t[i][j][d]=0;

	/* return pointer to array of pointer to rows */
	return t;
}


float ****f4tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh, int nvl, int nvh){
	/* allocate a float 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]=0.0 */
	int i,j,d,v, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1, nval=nvh-nvl+1;
	float ****t;

	/* allocate pointers to pointers to rows */
	t=(float ****) malloc((size_t) ((nrow+NR_END)*sizeof(float**)));
	if (!t) err("allocation failure 1 in function f4tensor() ");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) err("allocation failure 2 in function f4tensor() ");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float **) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) err("allocation failure 3 in function f4tensor() ");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	/* allocate values and set pointers to them */
	t[nrl][ncl][ndl]=(float *) malloc((size_t)((nrow*ncol*ndep*nval+NR_END)*sizeof(float)));
	if (!t[nrl][ncl][ndl]) err("allocation failure 4 in function f4tensor() ");
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= nvl;

	/*for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;*/
	for (d=ndl+1;d<=ndh;d++) t[nrl][ncl][d]=t[nrl][ncl][d-1]+nval;

	for (i=nrl+1;i<=nrh;i++){
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		t[i][ncl][ndl]=t[i-1][ncl][ndl]+ncol*ndep*nval;
		for (j=ncl+1;j<=nch;j++){
			t[i][j]=t[i][j-1]+ndep;
			t[i][j][ndl]=t[i][j-1][ndl]+ndep*nval;
			for (d=ndl+1;d<=ndh;d++){
				t[i][j][d]=t[i][j][d-1]+nval;
			}
		}
	}

	/* initializing 4tensor */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++)
			for (d=ndl;d<=ndh;d++)
				for (v=nvl;d<=nvh;v++) t[i][j][d][v]=0.0;

	/* return pointer to array of pointer to rows */
	return t;
}





void free_vector(float *v, int nl, int nh){
	/* free a float vector allocated with vector() */
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, int nl, int nh){
	/* free a int vector allocated with vector() */
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(char *v, int nl, int nh){
	/* free a char vector allocated with vector() */
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, int nl, int nh){
        /* free a char vector allocated with vector() */
	free((FREE_ARG) (v+nl-NR_END));
        }
                
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch){
	/* free a float matrix allocated by matrix() */
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch){
	/* free a integer matrix allocated by imatrix() */
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_usmatrix(unsigned short int **m, int nrl, int nrh, int ncl, int nch){
	/* free a integer matrix allocated by imatrix() */
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh){
	/* free a float matrix allocated by f3tensor() */
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_i3tensor(int ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh){
	/* free a float matrix allocated by i3tensor() */
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void zero(float *A, int u_max){
	int u=0;
	while(++u<=u_max) *(A++)=0.0;
}

void normalize_data(float **data, int ntr, int ns){

	float *max_min_trace=NULL;
	float max, min, tmp;
	int i,j;


	max_min_trace = vector(1,ntr); 
	for(i=1;i<=ntr;i++){

      	    max=0.0;
      	    min=0.0;

     	    for(j=2;j<=ns;j++){

      	       /* Estimate max and min */
               if(data[i][j]>max){max = data[i][j];}
               if(data[i][j]<min){min = data[i][j];}

	       if (max>(fabs(min))){tmp=max;}
	       else{tmp=fabs(min);}

            }
 
	    max_min_trace[i]=tmp;	

	    /* set maximum_values=0.0 to 1.0 */
	    if(max_min_trace[i]==0.0){max_min_trace[i]=1.0;}
   

	    for(j=2;j<=ns;j++){
		data[i][j] = data[i][j]/max_min_trace[i];
	    }
	
	}

        free_vector(max_min_trace,1,ntr);
}

void normalize_data_rel(float **data, float **data_obs, int ntr, int ns){

	float *max_min_trace=NULL, *max_min_obs_trace=NULL;
	float max, min, max_obs, min_obs, tmp;
	int i,j;

	max_min_trace = vector(1,ntr); 
	max_min_obs_trace = vector(1,ntr); 

	for(i=1;i<=ntr;i++){

      	    max=0.0;
      	    min=0.0;

      	    max_obs=0.0;
      	    min_obs=0.0;
        
     	    for(j=2;j<=ns;j++){

      	       /* Estimate max and min of model and field data */
               if(data[i][j]>max){max = data[i][j];}
               if(data[i][j]<min){min = data[i][j];}

               if(data_obs[i][j]>max_obs){max_obs = data_obs[i][j];}
               if(data_obs[i][j]<min_obs){min_obs = data_obs[i][j];}

            }

	    if (max>(fabs(min))){tmp=max;}
	    else{tmp=fabs(min);}

            max_min_trace[i]=tmp;

	    if (max_obs>(fabs(min_obs))){tmp=max_obs;}
	    else{tmp=fabs(min_obs);}

	    max_min_obs_trace[i]=tmp;	

	    /* set maximum_values=0.0 to 1.0 */
	    if(max_min_trace[i]==0.0){max_min_trace[i]=1.0;}
            if(max_min_obs_trace[i]==0.0){max_min_obs_trace[i]=1.0;}

	    for(j=2;j<=ns;j++){
		data_obs[i][j] = data_obs[i][j]*(max_min_trace[i]/max_min_obs_trace[i]);
	    }
	
	}

        free_vector(max_min_trace,1,ntr);
        free_vector(max_min_obs_trace,1,ntr);

}

void normalize_data_rms(float **data, float **data_obs, int ntr, int ns){

	float rms, rms_obs;
	int i,j;
 
        rms=0.0;
      	rms_obs=0.0;

	for(i=1;i<=ntr;i++){        
     	    for(j=2;j<=ns;j++){

      	       /* Estimate rms-values of model and field data */
               rms += fabs(data[i][j]);
               rms_obs += fabs(data_obs[i][j]);
               
            }

	}

        rms=rms/(ntr*(ns-1));
        rms_obs=rms_obs/(ntr*(ns-1));

        /* set rms-values=0.0 to 1.0 */
	if(rms==0.0){rms=1.0;}
        if(rms_obs==0.0){rms_obs=1.0;}

	for(i=1;i<=ntr;i++){        

	    for(j=2;j<=ns;j++){
		data_obs[i][j] = data_obs[i][j]*(rms/rms_obs);
	    }

	}

}

void normalize_STF(float *data, int ns){

	float max;
	int j;

        max=fabs(data[1]);

     	for(j=2;j<=ns;j++){

      	   /* Estimate maximum */
           if(fabs(data[j])>max){max = fabs(data[j]);}
           

	}
   
	for(j=1;j<=ns;j++){
	   data[j] = data[j]/max;
	}
	
}

void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Input:
d[]     the diagonal of A
e[]     the superdiagonal of A
b[]     the rhs of Ax=b

Output:
b[]     b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
*****************************************************************************/
{
        int k;
        float temp;
       
        /* decomposition */
        for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
        }

        /* substitution */
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];
       
        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k];
       
}

