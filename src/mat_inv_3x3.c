/*-------------------------------------------------------------
 *
 * Invert 3x3 square matrix 
 *  
 * Daniel Koehn
 *
 * Kiel, 31.10.2018
 *-------------------------------------------------------------
*/

#include "fd.h"

void mat_inv_3x3(float **A, float **Ainv) {
	
	/* local variables */
	float det = 0;
	float a, b, c, d, e, f, g, h, i;
		
	a = A[1][1]; b = A[1][2]; c = A[1][3];
	d = A[2][1]; e = A[2][2]; f = A[2][3]; 
	g = A[3][1]; h = A[3][2]; i = A[3][3];
	
	/* calculate determinant */
	det = a * e * i + b * f * g + c * d * h - c * e * g - a * f * h - b * d * i;
	
	/* calculate inverse matrix */
	Ainv[1][1] = (e * i - f * h) / det; Ainv[1][2] = (c * h - b * i) / det; Ainv[1][3] = (b * f - c * e) / det;
	Ainv[2][1] = (f * g - d * i) / det; Ainv[2][2] = (a * i - c * g) / det; Ainv[2][3] = (c * d - a * f) / det;
	Ainv[3][1] = (d * h - e * g) / det; Ainv[3][2] = (b * g - a * h) / det; Ainv[3][3] = (a * e - b * d) / det;  
		
}
