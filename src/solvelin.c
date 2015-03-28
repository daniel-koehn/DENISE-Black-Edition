/*
 *------------------------------------------------------------------------
 *
 *   Solve linear equation systems
 *
 *------------------------------------------------------------------------
*/


#include "fd.h"

void solvelin(float  **AA, float *bb, float *x, int e, int method)
{

	/* local variables */
	int k, m, n, rows, columns;
	float a, c, **A, *b;

	
	rows = e;
	columns = e;

	A = matrix(1,rows,1,columns);
	b = vector(1,rows);

	/* temporary variables */
	for (k=1;k<=rows;k++) {
		for (n=1;n<=columns;n++)  A[k][n] = AA[k][n];
		b[k] = bb[k];
	}


	switch (method)
	{
	case 1:	/* Gauﬂ algorithm */
	{
		for (k=1;k<=rows-1;k++)
			for (n=k;n<=rows-1;n++)
			{
				a = A[n+1][k]/A[k][k];
				for (m=1;m<=columns;m++) A[n+1][m] = A[n+1][m] - a*A[k][m];
				b[n+1] = b[n+1] - a*b[k];
			}
		
		for (k=rows;k>=1;k--)
		{
			c = b[k];
			for (m=columns;m>=k+1;m--) c = c - A[k][m]*x[m];
			x[k] = c/A[k][k];
		}
		break;
	} /* END of case Gauﬂ */
		
	} /* END of switch (method) */
	
	
	free_matrix(A,1,rows,1,columns);
	free_vector(b,1,rows);

	return;
}




