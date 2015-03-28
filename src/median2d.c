/*****************************************************************************/
/*         calculate the median of a 2D matrix                               */
/*****************************************************************************/


#include "fd.h"

int minpos(float *a, int n);
float vecmax(float *a, int n);


float median2d(float **mat, int ny, int nx){

	int i, j, k;
	float *t, *s, maxglob, med;
	

	t = vector(1,nx*ny);
	s = vector(1,nx*ny);
	
	k = 0;
	for (j=1;j<=ny;j++)
		for (i=1;i<=nx;i++)
		{
			k++;
			t[k] = mat[j][i];
		}
		
	maxglob = vecmax(t,k);
	
	for (j=1;j<=k;j++)
	{
		for (i=1;i<=k;i++)
		{
			s[j] = t[minpos(t,k)];
			t[minpos(t,k)] = maxglob;
		}
	}
	
	med = s[(k+1)/2];
	
	free_vector(t,1,k);
	free_vector(s,1,k);
	
	return med;
}


int minpos(float *a, int n){
	float mini;
	int j;
	int minp = 1;

	mini = a[1];
	for (j=2;j<=n;j++)
	{	
		if (a[j]<mini)
		{
			mini = a[j];
			minp = j;
		}
	}
		
	return minp;
}


float vecmax(float *a, int n){
	float maxi;
	int j;

	maxi = a[1];
	for (j=2;j<=n;j++)
	{	
		if (a[j]>maxi)
		{
			maxi = a[j];
		}
	}
		
	return maxi;
}
