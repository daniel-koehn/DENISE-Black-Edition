/*------------------------------------------------------------------------
 *  fd.h - header file for DENISE
 *
 *  Daniel Koehn
 *  Kiel, 24.07.2016
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 100
#define STRING_SIZE 74
#define STRING_SIZE2 256
#define REQUEST_COUNT 4

/* ---------------------------------- */
/* declaration of PSV data-structures */
/* ---------------------------------- */

/* PSV (visco)-elastic wavefield variables */
struct wavePSV{
   float  ** pvx, ** pvy, **  pvxp1, **  pvyp1, **  pvxm1, **  pvym1;
   float  ** psxx, **  psxy, **  psyy, ** ux, ** uy, ** uxy, ** uyx, ** uttx, ** utty;
   float *** pr, ***pp, ***pq;
} wavePSV; 

/* PSV PML variables*/
struct wavePSV_PML{
   float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half, * alpha_prime_x_half; 
   float * a_x_half, * b_x_half, * d_y, * K_y, * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half; 
   float * alpha_prime_y_half, * a_y_half, * b_y_half, ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y; 
   float ** psi_sxy_x, ** psi_vxx, ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;
   float  **  absorb_coeff;
} wavePSV_PML;

/* PSV material parameters */
struct matPSV{
   float  **prho, **prip, **prjp, **ppi, **pu, **puipjp;
   float **ptaus, **ptaup, *etaip, *etajm, *peta, **ptausipjp, **fipjp, ***dip, *bip, *bjm;
   float *cip, *cjm, ***d, ***e, **f, **g;
} matPSV;

/* PSV FWI variables */
struct fwiPSV{
   float  **  prho_old, **pu_old, **ppi_old;
   float  ** Vp0, ** Vs0, ** Rho0;
   float  **waveconv, **waveconv_lam, **waveconv_mu, **waveconv_rho, **waveconv_rho_s, **waveconv_u;
   float **waveconv_shot, **waveconv_u_shot, **waveconv_rho_shot;
   float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u;
   float  *forward_prop_x, *forward_prop_y, *forward_prop_rho_x, *forward_prop_u, *forward_prop_rho_y;
} fwiPSV;

/* PSV seismogram variables */
struct seisPSV{
   float ** sectionvx, ** sectionvy, ** sectionp, ** sectioncurl, ** sectiondiv;
   float ** fulldata, ** fulldata_vx, ** fulldata_vy;
   float ** fulldata_p, ** fulldata_curl,  ** fulldata_div;
} seisPSV;

/* PSV seismogram variables for FWI */
struct seisPSVfwi{
   float ** sectionvxdata, ** sectionvxdiff, ** sectionvxdiffold, ** sectionvydiffold;
   float ** sectionvydiff, ** sectionvydata, ** sectionread;
   float ** sectionpdata, ** sectionpdiff, ** sectionpdiffold;
   float energy;
   float L2;
} seisPSVfwi;

/* Acquisition geometry */
struct acq{
   int * recswitch, ** recpos, ** recpos_loc;
   float ** srcpos, **srcpos_loc, ** srcpos1;
   float ** srcpos_loc_back, ** signals;
} acq;

/* PSV MPI variables */
struct mpiPSV{
   float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top;
} mpiPSV;

/* ---------------------------------- */
/* declaration of VTI data-structures */
/* ---------------------------------- */

/* VTI material parameters */
struct matVTI{
   float  **prho, **prip, **prjp, **c11, **c13, **c33, **c44, **c44h;
} matVTI;

/* ------------- */
/* PSV functions */
/* ------------- */

void alloc_fwiPSV(struct fwiPSV *fwiPSV);

void alloc_matPSV(struct matPSV *matPSV);

void alloc_mpiPSV(struct mpiPSV *mpiPSV);

void alloc_seisPSV(int ntr, int ns, struct seisPSV *seisPSV);

void alloc_seisPSVfull(struct seisPSV *seisPSV, int ntr_glob);

void alloc_seisPSVfwi(int ntr, int ntr_glob, int ns, struct seisPSVfwi *seisPSVfwi);

void alloc_PSV(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML);

void ass_gradPSV(struct fwiPSV *fwiPSV, struct matPSV *matPSV, int iter);

float calc_mat_change_test_PSV(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float  **  rhonp1, float **  pi, float **  pinp1, float **  u, float **  unp1, 
int iter, int epstest, float eps_scale, int itest);

void calc_res_PSV(struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob,  int ntr, int nsrc_glob, float ** srcpos, int ishot, int ns, int iter,
                  int swstestshot);

void dealloc_PSV(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML);

void exchange_s_PSV(float ** sxx, float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top,
MPI_Request * req_send, MPI_Request * req_rec);

void exchange_v_PSV(float ** vx, float ** vy,  
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top,
MPI_Request * req_send, MPI_Request * req_rec);

void extract_LBFGS_PSV( int iter, float ** waveconv, float ** gradp, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float **ppi, float ** pu, float ** prho, float * r_LBFGS);

void extract_PCG_PSV(float * PCG_old, float ** waveconv, float ** waveconv_u, float ** waveconv_rho);

void FD_PSV();

void FWI_PSV();

float grad_obj_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int iter, int nsrc, int ns, int ntr, int ntr_glob, int nsrc_glob, 
         int nsrc_loc, int ntr_loc, int nstage, float **We, float **Ws, float **Wr, float ** taper_coeff, int hin, int *DTINV_help, 
         MPI_Request * req_send, MPI_Request * req_rec);

void matcopy_PSV(float ** prho, float ** ppi, float ** pu, float ** ptaup,
float ** ptaus);

void matcopy_elastic_PSV(float ** prho, float ** ppi, float ** pu);

void mem_fwiPSV(int nseismograms,int ntr, int ns, int fdo3, int nd, float buffsize, int ntr_glob);

void mem_PSV(int nseismograms,int ntr, int ns, int fdo3, int nd, float buffsize);

void model_freq_out_PSV(float ** ppi, float  **  rho, float **  pu, int iter, float freq);

float obj_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int nsrc, int nsrc_loc, int nsrc_glob, int ntr, int ntr_glob, 
int ns, int itest, int iter, float **Ws, float **Wr, int hin, int *DTINV_help, float eps_scale, MPI_Request * req_send, MPI_Request * req_rec);

void outseis_PSVfor(struct seisPSV *seisPSV, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int iter, FILE *FP);

void outseis_PSVres(struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int nstage, FILE *FP);

void physics_PSV();

void precond_PSV(struct fwiPSV *fwiPSV, struct acq *acq, int nsrc, int ntr_glob, float ** taper_coeff, FILE *FP_GRAV);

void prepare_update_s_visc_PSV(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
float **puipjp, float **ppi, float **prho, float **ptaus, float **ptaup,
float **ptausipjp, float **f, float **g, float *bip, float *bjm,
float *cip, float *cjm, float ***dip, float ***d, float ***e);

void psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, int hin, int *DTINV_help, int mode, MPI_Request * req_send, MPI_Request * req_rec);

void readmod_visc_PSV(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta);

void readmod_elastic_PSV(float  **  rho, float **  pi, float **  u);

void RTM_PSV();

void RTM_PSV_out(struct fwiPSV *fwiPSV);

void RTM_PSV_out_shot(struct fwiPSV *fwiPSV, int ishot);

float step_length_est_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int iter, int nsrc, int ns, int ntr, int ntr_glob, float * epst1, 
         float * L2t, int nsrc_glob, int nsrc_loc, int *step1, int *step3, int nxgrav, int nygrav, int ngrav, float **gravpos, float *gz_mod, int NZGRAV, int ntr_loc, 
         float **Ws, float **Wr, int hin, int *DTINV_help, MPI_Request * req_send, MPI_Request * req_rec);

void stf_psv(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matPSV *matPSV, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, struct seisPSV *seisPSV, 
             struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, int nsrc, int ns, int ntr, int ntr_glob, int iter, float **Ws, 
             float **Wr, int hin, int *DTINV_help, MPI_Request * req_send, MPI_Request * req_rec);

void surface_elastic_PML_PSV(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float  **  pi, float  **  u, float ** rho, 
			     float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx);

void surface_visc_PML_PSV(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float ***p, float ***q, float  **  ppi, float  **  pu, 
			  float **prho, float **ptaup, float **ptaus, float *etajm, float *peta, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx);

void store_LBFGS_PSV(float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float ** waveconv, float ** gradp, float ** waveconv_u, 
float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float **ppi, float ** pu, float ** prho, int nxnyi, 
int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);

void store_PCG_PSV(float * PCG_old, float ** waveconv, float ** waveconv_u, float ** waveconv_rho);

void update_s_elastic_PML_PSV(int nx1, int nx2, int ny1, int ny2,
float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float **rho, float *hc, int infoout,
float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int sws);

void update_s_visc_PML_PSV(int nx1, int nx2, int ny1, int ny2,
float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
float **   sxy, float ** pi, float ** u, float ** uipjp, float **rho, float *hc, int infoout,
float ***r, float ***p, float ***q, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, 
float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int mode);

void update_v_PML_PSV(int nx1, int nx2, int ny1, int ny2, int nt,
float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float **  uttx, float **  utty,float ** sxx, float ** syy,
float ** sxy, float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_sxx_x, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_syx_x);

void zero_denise_elast_PSV(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float ** vxm1, float ** vym1, 
float ** vxym1, float ** vxp1, float ** vyp1, float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, 
float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs);

void zero_denise_visc_PSV(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float ** vxm1, float ** vym1, 
float ** vxym1, float ** vxp1, float ** vyp1, float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, 
float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ***pr, float ***pp, float ***pq);

/* ------------- */
/* VTI functions */
/* ------------- */

void alloc_matVTI(struct matVTI *matVTI);

void av_c44(float ** c44, float ** c44h);

void checkfd_ssg_VTI(FILE *fp, float ** prho, float ** c11, float ** c13, float ** c33, float ** c44, float *hc);

void FD_VTI();

void matcopy_elastic_VTI(float ** rho, float ** pi, float ** u);

void model_elastic_VTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44);

void physics_VTI();

void readmod_elastic_VTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44);

void seismo_ssg_vti(int lsamp, int ntr, int **recpos, float **sectionvx, 
float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
float **vx, float **vy, float **sxx, float **syy, float *hc);

void snap_vti(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx, float **syy, float *hc);

void update_s_elastic_PML_VTI(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** c11,  float ** c13, float ** c33, float ** c44h, float ** absorb_coeff, float *hc, int infoout,
        float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int mode);

void vti(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matVTI *matVTI, struct fwiPSV *fwiPSV, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct seisPSVfwi *seisPSVfwi, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, int hin, int *DTINV_help, int mode, MPI_Request * req_send, MPI_Request * req_rec);

/* ----------------- */
/* General functions */
/* ----------------- */

void window_cos(float **win, int npad, int nsrc, float it1, float it2, float it3, float it4);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, MPI_Comm newcomm_nodentr);

int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch);

void absorb(float ** absorb_coeff);

void smooth_model(float ** pinp1, float ** unp1, float ** rho, int iter);

void taper_grad(float ** waveconv, float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int sws);

void taper_grad_shot(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void spat_filt(float ** waveconv, int iter, int sws);

void apply_tdfilt(float **section, int ntr, int ns, int order, float fc2, float fc1);

void av_mat(float **  pi, float **  u, 
float **  ppijm, float **  puip, float ** pujm);

void av_mue(float ** u, float ** uipjp, float ** rho);

void av_rho(float **rho, float **rip, float **rjp);

void av_tau(float **taus, float **tausipjp);

float median2d(float **mat, int ny, int nx);

void calc_envelope(float ** datatrace, float ** envelope, int ns, int ntr);

void calc_hilbert(float ** datatrace, float ** envelope, int ns, int ntr);

double calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos, int **recpos_loc, float **srcpos, int nsrc_glob, int ishot, int iter);

double calc_res_grav(int ngrav, float *gz_mod, float *gz_res);

double calc_misfit(float **sectiondiff, int ntr, int ns, int LNORM, float L2, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot);

float calc_opt_step(float *  L2t, float * epst, int sws);

double calc_energy(float **sectiondata, int ntr, int ns, float energy, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot);

void checkfd_ssg_elastic(FILE *fp, float ** prho, float ** ppi, float ** pu, float *hc);

void checkfd_ssg_visc(FILE *fp, float ** prho, float ** ppi, float ** pu, float ** ptaus, float ** ptaup, float * peta, float *hc);

void check_mode_phys();

void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, 
MPI_Request *req_send, MPI_Request *req_rec);

void conv_FD(float * temp_TS, float * temp_TS1, float * temp_conv, int ns);

void copy_mat(float ** A, float ** B);

void descent(float ** grad, float ** gradm);

float dotp(float * vec1, float *vec2, int n1, int n2, int sw);

float exchange_L2(float L2, int sw, int bcast_l2);

void eprecond(float ** W, float ** vx, float ** vy);

void eprecond1(float ** We, float ** Ws, float ** Wr);

void extend_mod(float  **rho_grav, float  **rho_grav_ext, int nxgrav, int nygrav);

void gauss_filt(float ** waveconv);

void grav_grad(int ngrav, float **gravpos, float **grad_grav, float *gz_res);

void grav_mod(float  **rho, int ngrav, float **gravpos, float *gz, int NXGRAV, int NYGRAV, int NZGRAV);

void read_back_density(float ** rho_back);

float *holbergcoeff(void);

int householder(int m, int n, float **mat, float *b);

void info(FILE *fp);

void init_grad(float ** A);

void initproc(void);

void interpol(int ni1, int ni2, float **  intvar, int cfgt_check);

void LBFGS(int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);
           
double LU_decomp(double  **A, double *x, double *b,int n);

float maximum_m(float **mat, int nx, int ny);

void median_src(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int iter, int sws);

float minimum_m(float **mat, int nx, int ny);

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta);

void model_elastic(float  **  rho, float **  pi, float **  u);
			  
void merge(int nsnap, int type);

void mergemod(char modfile[STRING_SIZE], int format);

void msource(int nt, float ** sxx, float ** syy, float ** sxy, float **  srcpos_loc, float ** signals, int nsrc, int sw);

void norm(float **waveconv);

void note(FILE *fp);

void  outseis(FILE *fp, FILE *fpdata, int comp, float **section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  outseis_glob(FILE *fp, FILE *fpdata, int comp, float **section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  outseis_vector(FILE *fp, FILE *fpdata, int comp, float *section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  inseis(int comp, float **section, int ntr, int ns, int sws, int iter);

void  taper(float **sectionpdiff, int ntr, int ns);

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void PCG(float * PCG_new, float * PCG_old, float * PCG_dir, int PCG_class);

void PML_pro(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
float * d_y, float * K_y, float * alpha_prime_y, float * a_y, float * b_y, 
float * d_y_half, float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half);

void psource(int nt, float ** sxx, float ** syy,
float **  srcpos_loc, float ** signals, int nsrc, int sw);

float *rd_sour(int *nts,FILE* fp_source);

void read_density_glob(float ** rho_grav, int sws);

float readdsk(FILE *fp_in, int format);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy);

float **read_grav_pos(int *ngrav);

void read_par(FILE *fp_in);

void read_par_inv(FILE *fp,int nstage,int stagemax);

int **receiver(FILE *fp, int *ntr, int ishot);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy);

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos_loc, int nsrc,int ns, int iter);

void saveseis_glob(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos_loc, int nsrc,int ns, int iter);

void snap(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx,
	float **syy, float **u, float **pi, float *hc);

void snapmerge(int nsnap);

float **sources(int *nsrc);

void solvelin(float  **AA, float *bb, float *x, int e, int method);

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx, 
float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu); 

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx, 
float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu,
float **prho, float *hc);

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc);

float **splitsrc_back(int **recpos,int *nsrc_loc, int nsrc);

void stalta(float **sectionp, int ntr, int nst, float *picked_times, int ishot);

void stf(float **sectionvy_obs, float **sectionvy, int ntr_glob, int ishot, int ns, int iter, int nshots, float **signals, int **recpos, float **srcpos);

void  timedomain_filt(float ** data, float fc, int order, int ntr, int ns, int method);
void  timedomain_filt_vector(float * data, float fc, int order, int ntr, int ns, int method);

void time_window(float **sectiondata, float * picked_times, int iter, int ntr_glob, int **recpos_loc, int ntr, int ns, int ishot);

void time_window_stf(float **sectiondata, int iter, int ntr_glob, int ns, int ishot);

void tripd(float *d, float *e, float *b, int n);

float ** wavelet(float ** srcpos_loc, int nsrc, int ishot);

float ** wavelet_stf(int nsrc, int ishot, float ** signals_stf);

void  wavenumber(float ** grad);

void write_par(FILE *fp);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float ** array, int format);

void zero_LBFGS(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS);

void zero_PCG(float * PCG_old, float * PCG_new, float * PCG_dir, int PCG_vec);
		 
void FLnode(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta);

void smooth_grad(float ** waveconv);

void  smooth2(float ** grad);

/* declaration of functions for parser*/

/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(FILE *fp, char input_file[STRING_SIZE],char ** varname_list,char ** value_list);

void print_objectlist_screen(FILE *fp, int number_readobject,char ** varname_list,char ** value_list);

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]);

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]);

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int * int_buffer,
		char ** varname_list,char ** value_list);

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float * double_buffer,
		char ** varname_list,char ** value_list);

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
		char ** varname_list,char ** value_list);

int is_string_blankspace(char string_in[STRING_SIZE]);

void remove_blankspaces_around_string(char string_in[STRING_SIZE] );

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int * number_read_object,
		char ** varname_list,char ** value_list );

/* utility functions */
void err(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh); 
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, 
int ndh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);


