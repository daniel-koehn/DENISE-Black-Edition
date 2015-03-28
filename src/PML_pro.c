/*
* Define damping profiles for CPML boundary condition
* This C-PML implementation is adapted from the 2nd order isotropic CPML code by Dimitri Komatitsch and based in part on formulas given in Roden and Gedney (2000). 
* Additionally the code is based on the following references:
* 
* @ARTICLE{KoMa07,
* author = {Dimitri Komatitsch and Roland Martin},
* title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
*          at grazing incidence for the seismic wave equation},
* journal = {Geophysics},
* year = {2007},
* volume = {72},
* number = {5},
* pages = {SM155-SM167},
* doi = {10.1190/1.2757586}}
*
* @ARTICLE{MaKoEz08,
* author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
* title = {An unsplit convolutional perfectly matched layer improved at grazing
* incidence for seismic wave equation in poroelastic media},
* journal = {Geophysics},
* year = {2008},
* volume = {73},
* pages = {T51-T61},
* number = {4},
* doi = {10.1190/1.2939484}}
*
* @ARTICLE{MaKo09,
* author = {Roland Martin and Dimitri Komatitsch},
* title = {An unsplit convolutional perfectly matched layer technique improved
* at grazing incidence for the viscoelastic wave equation},
* journal = {Geophysical Journal International},
* year = {2009},
* volume = {179},
* pages = {333-344},
* number = {1},
* doi = {10.1111/j.1365-246X.2009.04278.x}}
*
* @ARTICLE{MaKoGe08,
* author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
* title = {A variational formulation of a stabilized unsplit convolutional perfectly
* matched layer for the isotropic or anisotropic seismic wave equation},
* journal = {Computer Modeling in Engineering and Sciences},
* year = {2008},
* volume = {37},
* pages = {274-304},
* number = {3}}
*
* The original CPML technique for Maxwell's equations is described in:
*
* @ARTICLE{RoGe00,
* author = {J. A. Roden and S. D. Gedney},
* title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
*          of the {CFS}-{PML} for Arbitrary Media},
* journal = {Microwave and Optical Technology Letters},
* year = {2000},
* volume = {27},
* number = {5},
* pages = {334-339},
* doi = {10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A}}
*/

#include "fd.h"

void PML_pro(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
            float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
            float * d_y, float * K_y, float * alpha_prime_y, float * a_y, float * b_y, 
            float * d_y_half, float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half)
{

	/* extern variables */

	extern float DH, DAMPING, DT, FPML;
	extern int FREE_SURF, NX, NY, BOUNDARY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3], FW;
	extern FILE *FP;
	
	/* local variables */
	int i, j, ifw, ii, jj, xb, yb, xe, ye, h;
	int USE_PML_XMIN, USE_PML_XMAX, USE_PML_YMIN, USE_PML_YMAX;
	float amp, a, *coeff, fc, cp;
      
      extern float npower, k_max_PML;
      /*const float npower = 2.0;*/ /*  power to compute d0 profile */
      /*const float k_max_PML = 1.0;*/ /* from Gedney page 8.11 */
      const float alpha_max_PML = 2.0 * PI * (FPML/2.0); /* from festa and Vilotte */

      float thickness_PML_x, thickness_PML_y, xoriginleft, xoriginright, yoriginbottom, yorigintop;
      float Rcoef , d0_x, d0_y, xval, yval, abscissa_in_PML, abscissa_normalized;

      char profile[STRING_SIZE];

      /* define profile of absorption in PML region */

      /* thickness of the PML layer in meters */
      thickness_PML_x = FW*DH;
      thickness_PML_y = FW*DH;

      /* reflection coefficient (INRIA report section 6.1) */
      Rcoef = 0.001;

      /* compute d0 from INRIA report section 6.1 */
      d0_x = - (npower + 1) * DAMPING * log(Rcoef) / (2.0 * thickness_PML_x);
      d0_y = - (npower + 1) * DAMPING * log(Rcoef) / (2.0 * thickness_PML_y);
	
      /* damping in the X direction */
      /* -------------------------- */

      /* origin of the PML layer (position of right edge minus thickness, in meters) */
      xoriginleft = thickness_PML_x;
      xoriginright = (NXG-1) * DH - thickness_PML_x;

      /* left boundary */

      for (i=1;i<=FW;i++){
         
          K_x[i] = 1.0;
          K_x_half[i] = 1.0;
          xval = DH * (i-1);

            /* define damping profile at the grid points */
            abscissa_in_PML = xoriginleft - xval;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x[i] = d0_x * pow(abscissa_normalized,npower);

            /* this taken from Gedney page 8.2 */
               K_x[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_x[i] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = xoriginleft - (xval + DH/2.0);

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x_half[i] = d0_x * pow(abscissa_normalized,npower);

               /* this taken from Gedney page 8.2 */
               K_x_half[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_x_half[i] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

       /* just in case, for -5 at the end */
       if(alpha_prime_x[i] < 0.0){ alpha_prime_x[i] = 0.0;}
       if(alpha_prime_x_half[i] < 0.0) {alpha_prime_x_half[i] = 0.0;}

       b_x[i] = exp(- (d_x[i] / K_x[i] + alpha_prime_x[i]) * DT);
       b_x_half[i] = exp(- (d_x_half[i] / K_x_half[i] + alpha_prime_x_half[i]) * DT);

       /* avoid division by zero outside the PML */
       if(abs(d_x[i]) > 1.0e-6){ a_x[i] = d_x[i] * (b_x[i] - 1.0) / (K_x[i] * (d_x[i] + K_x[i] * alpha_prime_x[i]));}
       if(abs(d_x_half[i]) > 1.0e-6){ a_x_half[i] = d_x_half[i] * (b_x_half[i] - 1.0) / (K_x_half[i] * (d_x_half[i] + K_x_half[i] * alpha_prime_x_half[i]));}

      } /* end of left boundary */ 

      /* right boundary */
           
      for (i=NXG-FW+1;i<=NXG;i++){
             
             h=i-NXG+2*FW;
             
             K_x[h] = 1.0;
	     K_x_half[h] = 1.0;
	     xval = DH * (i-1);                                 
              
            /* define damping profile at the grid points */
            abscissa_in_PML = xval - xoriginright;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x[h] = d0_x * pow(abscissa_normalized,npower);

            /* this taken from Gedney page 8.2 */
               K_x[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_x[h] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = xval + DH/2.0 - xoriginright;

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x_half[h] = d0_x * pow(abscissa_normalized,npower);

               /* this taken from Gedney page 8.2 */
               K_x_half[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_x_half[h] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

            /* just in case, for -5 at the end */
               if(alpha_prime_x[h] < 0.0){ alpha_prime_x[h] = 0.0;}
               if(alpha_prime_x_half[h] < 0.0) {alpha_prime_x_half[h] = 0.0;}

               b_x[h] = exp(- (d_x[h] / K_x[h] + alpha_prime_x[h]) * DT);
               b_x_half[h] = exp(- (d_x_half[h] / K_x_half[h] + alpha_prime_x_half[h]) * DT);

            /* avoid division by zero outside the PML */
               if(abs(d_x[h]) > 1.0e-6){ a_x[h] = d_x[h] * (b_x[h] - 1.0) / (K_x[h] * (d_x[h] + K_x[h] * alpha_prime_x[h]));}
               if(abs(d_x_half[h]) > 1.0e-6){ a_x_half[h] = d_x_half[h] * (b_x_half[h] - 1.0) / (K_x_half[h] * (d_x_half[h] + K_x_half[h] * alpha_prime_x_half[h]));}

       } /* end of right boundary */   



	
      /* damping in the Y direction */
      /* -------------------------- */

      /* origin of the PML layer (position of right edge minus thickness, in meters) */
      yoriginbottom = thickness_PML_y;
      yorigintop = (NYG-1) * DH - thickness_PML_y;

      for (i=1;i<=FW;i++){
          
          K_y[i] = 1.0;
          K_y_half[i] = 1.0;
          yval = DH * (i-1);

          /* left boundary */

            /* define damping profile at the grid points */
            abscissa_in_PML = yoriginbottom - yval;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_y;
               d_y[i] = d0_y * pow(abscissa_normalized,npower);

            /* this taken from Gedney page 8.2 */
               K_y[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_y[i] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = yoriginbottom - (yval + DH/2.0);

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_y;
               d_y_half[i] = d0_y * pow(abscissa_normalized,npower);

               /* this taken from Gedney page 8.2 */
               K_y_half[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_y_half[i] = alpha_max_PML * (1.0 - abscissa_normalized);
            }
          b_y[i] = exp(- (d_y[i] / K_y[i] + alpha_prime_y[i]) * DT);
          b_y_half[i] = exp(- (d_y_half[i] / K_y_half[i] + alpha_prime_y_half[i]) * DT);

          /* avoid division by zero outside the PML */
          if(abs(d_y[i]) > 1.0e-6){ a_y[i] = d_y[i] * (b_y[i] - 1.0) / (K_y[i] * (d_y[i] + K_y[i] * alpha_prime_y[i]));}
          if(abs(d_y_half[i]) > 1.0e-6){ a_y_half[i] = d_y_half[i] * (b_y_half[i] - 1.0) / (K_y_half[i] * (d_y_half[i] + K_y_half[i] * alpha_prime_y_half[i]));}

      } /* end of left boundary */ 

      /* top boundary */
      for (i=NYG-FW+1;i<=NYG;i++){  
           
            h=i-NYG+2*FW;
            
            K_y[h] = 1.0;
            K_y_half[h] = 1.0;
            yval = DH * (i-1);

            /* define damping profile at the grid points */
            abscissa_in_PML = yval - yorigintop;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_y;
               d_y[h] = d0_y * pow(abscissa_normalized,npower);

            /* this taken from Gedney page 8.2 */
               K_y[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_y[h] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = yval + DH/2.0 - yorigintop;

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_y;
               d_y_half[h] = d0_x * pow(abscissa_normalized,npower);

               /* this taken from Gedney page 8.2 */
               K_y_half[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_normalized,npower);
               alpha_prime_y_half[h] = alpha_max_PML * (1.0 - abscissa_normalized);
            }

          b_y[h] = exp(- (d_y[h] / K_y[h] + alpha_prime_y[h]) * DT);
          b_y_half[h] = exp(- (d_y_half[h] / K_y_half[h] + alpha_prime_y_half[h]) * DT);

          /* avoid division by zero outside the PML */
          if(abs(d_y[h]) > 1.0e-6){ a_y[h] = d_y[h] * (b_y[h] - 1.0) / (K_y[h] * (d_y[h] + K_y[h] * alpha_prime_y[h]));}
          if(abs(d_y_half[h]) > 1.0e-6){ a_y_half[h] = d_y_half[h] * (b_y_half[h] - 1.0) / (K_y_half[h] * (d_y_half[h] + K_y_half[h] * alpha_prime_y_half[h]));}
      
       } /* end of top boundary */   



      
      
      /*if(MYID==0){
      sprintf(profile,"Kx.dat");
      FP=fopen(profile,"w");
      
    
      for (i=1;i<=2*FW+1;i=i++){
       fprintf(FP,"%e \n",K_y[i]);
      }
      fclose(FP);
      
      sprintf(profile,"Kx_half.dat");                                  
      FP=fopen(profile,"w");  
                    

      for (i=1;i<=2*FW+1;i++){
	fprintf(FP,"%e \n",K_y_half[i]);                          
      }
      fclose(FP); 

      sprintf(profile,"ax.dat");                                  
      FP=fopen(profile,"w");  
                    
     
      for (i=1;i<=2*FW+1;i++){
	fprintf(FP,"%e \n",a_y[i]);                          
      }
      fclose(FP); 
      
      sprintf(profile,"bx.dat");                                  
      FP=fopen(profile,"w");  
                    
   
      for (i=1;i<=2*FW+1;i++){
	fprintf(FP,"%e \n",b_y[i]);                          
      }
      fclose(FP);
              
      }*/
      

	/* compute coefficients for left and right grid boundaries (x-direction) */
	/*if ((!BOUNDARY) && (POS[1]==0))
	{
		yb=1; ye=NY; 
		for (i=1;i<=ifw;i++)
		{
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][i]=coeff[i];
		}
	}
			
	if ((!BOUNDARY) && (POS[1]==NPROCX-1))
	{
		yb=1; ye=NY;
		for (i=1;i<=ifw;i++){
			ii=NX-i+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][ii]=coeff[i];
		}
	}*/
	

	/* compute coefficients for top and bottom grid boundaries (y-direction) */

	/*if ((POS[2]==0) && (!(FREE_SURF)))
	{
		xb=1; xe=NX;
		for (j=1;j<=ifw;j++)
		{
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[j][i]=coeff[j];
		}
	}

	if (POS[2]==NPROCY-1)
	{
		xb=1; xe=NX;
		for (j=1;j<=ifw;j++)
		{
			jj=NY-j+1;
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[jj][i]=coeff[j];
		}
	}*/


/*	sprintf(modfile,"absorb_coeff.bin");

	writemod(modfile,absorb_coeff,3); 

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); 
*/

}



