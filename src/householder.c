//---------------------------------------------------------------------------
//   {Householdertransformation fÜr lineare Systeme}
//   {Lösung überbestimmter linearer Gleichungssysteme mit
//    Householdertransformation}

 //====================================================================//
 //  Die Funktion house dient zur Loesung eines linearen Gleichungs-   //
 //  systems:  mat // x = b.                                           //
 //                                                                    //
 //  Dabei sind: m die Zeilenzahl der Matrix,                          //
 //              n die Spaltenzahl der Matrix,                         //
 //              mat die m x n Koeffizientenmatrix des Systems,        //
 //              m >= n, mit rang (mat) = n,                           //
 //              b die rechte Seite des Systems (m-Vektor),            //
 //              x der Loesungsvektor des Gleichungssystems (n-Vektor).//
 //                                                                    //
 //  house benutzt die Householder Transformation zur Loesung des      //
 //  linearen, ueberbestimmten Gleichungssystems.                      //
 //  x ist die Loesung des Minimierungsproblems mat // x - b im Sinne  //
 //  der euklidschen Norm. Diese Loesung des Minimierungproblems muss  //
 //  nicht notwendigerweise auch eine Loesung des Gleichungssystems    //
 //  sein (Pruefung !).                                                //
 //====================================================================//
 //                                                                    //
 //   Eingabeparameter:                                                //
 //   ================                                                 //
 //      m        int m;  ( m > 0 )                                    //
 //               Anzahl der Zeilen von mat                            //
 //      n        int n;  ( n > 0 )                                    //
 //               Anzahl der Spalten von mat                           //
 //      mat      double                                               //
 //               Matrix des Gleichungssystems:                        //
 //                             mat[i][j], i = 0(1)m-1, j = 0(1)n-1.   //
 //      b        double   b[m];                                       //
 //               Rechte Seite des Gleichungssystems.                  //
 //                                                                    //
 //    mat und b werden in der Rechnung ueberschrieben !               //
 //                                                                    //
 //   Ausgabeparameter:                                                //
 //   ================                                                 //
 //      b        double   b[n];                                       //
 //               Loesungsvektor des Systems.                          //
 //                                                                    //
 //   Rueckgabewert:                                                   //
 //   =============                                                    //
 //      = 0      alles ok                                             //
 //      = 1      m oder n < 1 oder m < n                              //
 //      = 2      mat hat nicht Rang n.                                //
 //      = 3      Matrix hat Rang < n (numerisch).                     //
 //      = 4      zu wenig Speicher.                                   //
 //                                                                    //
 //====================================================================//

#include "fd.h"

int householder(int m, int n, float **mat, float *b){

  int i, j, k;
  float f, r, alpha, ak, eps, tmp, sum, norm, maxnorm;
  float * d;

  d = vector(1,n);
  
  if ((m < 1) || (n < 1) || (m < n)) return (1);

  //eps = (double)(2.0 * MACH_EPS);
  eps = (double)(2.0 * 1);

              
  for (i = 0; i < n; i++)            /*  Householder Transformation   */
  {
    r = 0.0;
    for (k = i; k < m; k++)
      r += mat[k][i] * mat[k][i];

    if (r == 0.0)                   /*  Matrix hat nicht Hoechstrang */
    {
     // vmfree(vmblock);
      return (2);
    }

    if (mat[i][i] >= 0.0)
      alpha = sqrt(r);
    else
      alpha = - sqrt(r);
                               
    ak = 1.0 / (r + alpha * mat[i][i]);
    mat[i][i] += alpha;

    d[i] = - alpha;

    maxnorm = 0.0;
    for (k = i + 1; k < n; k++)
    {
      norm = f = 0.0;
      for (j = i; j < m; j++)
      {
        tmp = mat[j][k];
        f += tmp * mat[j][i];
        norm += tmp * tmp;
      }

      if (norm > maxnorm)
        maxnorm = norm;

      f *= ak;
      for (j = i; j < m; j++)
        mat[j][k] -= f * mat[j][i];
    }

 /*   if (abs(alpha) < eps * sqrt(maxnorm))      // Loesbar ?
    {
      //vmfree(vmblock);
      return (3);
    }
*/

    for(f = 0, j = i; j < m; j++) // Rechte Seite transformieren
      f += b[j] * mat[j][i];

    f *= ak;
    for (j = i; j < m; j++)
      b[j] -= f * mat[j][i];

  } // ende for i

  for (i = n - 1; i >= 0; i--)        // Rueckwaertselimination       
  {
    sum = 0;
    for (k = i + 1; k < n; k++)
      sum += mat[i][k] * b[k];


    b[i] = (b[i] - sum) / d[i];
  }
                            
  
  free_vector(d,1,n);

  return (0);
       	
}



