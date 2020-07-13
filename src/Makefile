# Makefile for DENISE Black-Edition

#--------------------------------------------------------
# edit here:

# source code for model generation
MODEL = model_ainos_visc.c
MODEL_EL = model_ainos.c
MODEL_AC = model_acoustic.c
MODEL_VTI = zinc_vti.c
MODEL_TTI = zinc_tti.c
EXEC= ../bin

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4 

# On CRAY T3E
# CC=cc

# On Linux NEC-cluster with Intel-MPI
CC=mpiicc
LFLAGS=-lm -lcseife -lfftw3 -ftz -lstdc++
CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
SFLAGS=-L./../libcseife -L$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/lib
IFLAGS=-I./../libcseife -I$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/include -I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
# CC=mpicc -Minfo=all -ta=multicore #tesla:cc60,managed
# LFLAGS=-lm -lcseife -lfftw3 -lstdc++
# CFLAGS=-O3 -w -D_FORTIFY_SOURCE=0
# SFLAGS=-L./../libcseife -L/scratch/dragon/amd/kazeiv/gpuhackathon2018/fftw_sw/fftw-3.3.7/lib/
# IFLAGS=-I./../libcseife -I./../include -I/scratch/dragon/amd/kazeiv/gpuhackathon2018/fftw_sw/fftw-3.3.7/include/

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
#CC=mpicc
#LFLAGS=-lm -lcseife -lfftw3 -lstdc++
#CFLAGS=-O3 -w -fno-stack-protector -D_FORTIFY_SOURCE=0
#SFLAGS=-L./../libcseife
#IFLAGS=-I./../libcseife -I./../include

# On HLRN system
#CC=mpcc
#LFLAGS=-lm

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ipo
#LFLAGS=-lmpi -lm -i-static

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $< $(IFLAGS)

VPATH = snapmerge:PSV:VTI:TTI:AC:SH:models
        
SNAPMERGE_SCR = \
	merge.c \
	read_par.c \
	readdsk.c \
	snapmerge.c \
	writedsk.c \
	util.c

DENISE := $(wildcard *.c)

# add PSV source files
DENISE += alloc_fwiPSV.c \
	  alloc_matPSV.c \
	  alloc_mpiPSV.c \
	  alloc_PSV.c \
 	  alloc_seisPSV.c \
	  alloc_seisPSVfull.c \
	  alloc_seisPSVfwi.c \
	  ass_gradPSV.c \
	  calc_mat_change_test_PSV.c \
	  calc_res_PSV.c \
	  dealloc_PSV.c \
	  exchange_s_PSV.c \
	  exchange_v_PSV.c \
	  extract_LBFGS_PSV.c \
	  extract_PCG_PSV.c \
	  FD_PSV.c \
	  FWI_PSV.c \
	  grad_obj_psv.c \
	  matcopy_PSV.c \
	  matcopy_elastic_PSV.c \
	  mem_fwiPSV.c \
	  mem_PSV.c \
	  $(MODEL) \
	  $(MODEL_EL) \
	  model_freq_out_PSV.c \
	  model_it_out_PSV.c \
	  obj_psv.c \
	  outseis_PSVfor.c \
	  outseis_PSVres.c \
	  physics_PSV.c \
	  precond_PSV.c \
	  prepare_update_s_visc_PSV.c \
	  psv.c \
	  readmod_elastic_PSV.c \
	  readmod_visc_PSV.c \
	  RTM_PSV.c \
	  RTM_PSV_out.c \
	  RTM_PSV_out_shot.c \
	  step_length_est_psv.c \
	  stf_psv.c \
	  store_LBFGS_PSV.c \
	  store_PCG_PSV.c \
	  surface_visc_PML_PSV.c \
	  surface_elastic_PML_PSV.c \
	  update_s_elastic_PML_PSV.c \
	  update_s_visc_PML_PSV.c \
	  update_v_PML_PSV.c \
	  zero_denise_elast_PSV.c \
	  zero_denise_visc_PSV.c

# add VTI source files
DENISE += alloc_matVTI.c \
	  ass_gradVTI.c \
	  checkfd_ssg_VTI.c \
	  FD_VTI.c \
	  grad_obj_VTI.c \
	  matcopy_elastic_VTI.c \
	  $(MODEL_VTI) \
	  physics_VTI.c \
	  readmod_elastic_VTI.c \
	  RTM_VTI.c \
	  seismo_ssg_VTI.c \
	  snap_ssg_VTI.c \
	  update_s_elastic_PML_VTI.c \
	  VTI.c

# add TTI source files
DENISE += alloc_matTTI.c \
	  checkfd_ssg_TTI.c \
	  FD_TTI.c \
          grad_obj_TTI.c \
	  $(MODEL_TTI) \
	  physics_TTI.c \
	  readmod_elastic_TTI.c \
	  rot_el_tensor_TTI.c \
	  RTM_TTI.c \
	  update_s_elastic_PML_TTI.c \
	  TTI.c

# add AC source files
DENISE += ac.c \
	  alloc_AC.c \
	  alloc_fwiAC.c \
	  alloc_matAC.c \
	  ass_gradAC.c \
	  calc_mat_change_test_AC.c \
	  checkfd_acoustic.c \
	  dealloc_AC.c \
	  exchange_p_AC.c \
	  exchange_v_AC.c \
	  extract_LBFGS_AC.c \
	  extract_PCG_AC.c \
	  FD_AC.c \
	  FWI_AC.c \
	  grad_obj_ac.c \
	  matcopy_acoustic_AC.c \
	  mem_fwiAC.c \
	  $(MODEL_AC) \
	  model_freq_out_AC.c \
	  model_it_out_AC.c \
	  obj_ac.c \
	  physics_AC.c \
	  precond_AC.c \
	  psource_AC.c \
	  readmod_AC.c \
	  RTM_AC.c \
	  RTM_AC_out.c \
	  RTM_AC_out_shot.c \
	  seismo_AC.c \
	  snap_AC.c \
	  step_length_est_ac.c \
	  stf_ac.c \
	  store_LBFGS_AC.c \
	  store_PCG_AC.c \
	  surface_acoustic_PML_AC.c \
	  update_s_acoustic_PML_AC.c \
	  update_v_PML_AC.c \
	  zero_denise_acoustic_AC.c

# add SH source files
DENISE += alloc_fwiSH.c \
	  alloc_matSH.c \
	  alloc_seisSH.c \
	  alloc_seisSHfull.c \
	  alloc_seisSHfwi.c \
	  alloc_SH.c \
	  ass_gradSH.c \
	  apply_inv_hessSH.c \
	  av_mu_SH.c \
	  init_grad_coeff.c \
	  inv_rho_SH.c \
	  calc_mat_change_test_SH.c \
	  calc_res_SH.c \
	  checkfd_elast_SH.c \
	  checkfd_visc_SH.c \
	  dealloc_SH.c \
	  eprecond_SH.c \
	  exchange_s_SH.c \
	  exchange_v_SH.c \
	  extract_LBFGS_SH.c \
	  extract_PCG_SH.c \
	  FD_SH.c \
	  FD_grad_SH.c \
	  FWI_SH.c \
	  grad_obj_sh.c \
	  matcopy_elastic_SH.c \
	  matcopy_SH.c \
	  mem_SH.c \
	  model_freq_out_SH.c \
	  model_it_out_SH.c \
	  obj_sh.c \
	  outseis_SHfor.c \
	  outseis_SHres.c \
	  physics_SH.c \
	  precond_SH.c \
	  prepare_update_s_visc_SH.c \
	  readmod_elastic_SH.c \
	  readmod_visc_SH.c \
	  RTM_SH_out_shot.c \
	  saveseis_glob_SH.c \
	  sh.c \
          step_length_est_sh.c \
	  stf_sh.c \
	  store_LBFGS_SH.c \
	  store_PCG_SH.c \
          store_pseudo_hess_SH.c \
	  update_s_elastic_PML_SH.c \
	  update_s_visc_PML_SH.c \
	  update_v_PML_SH.c \
	  zero_denise_elast_SH.c \
	  zero_denise_visc_SH.c

SNAPMERGE_OBJ = $(SNAPMERGE_SCR:.c=.o)
DENISE_OBJ = $(DENISE:.c=.o)      

snapmerge:	$(SNAPMERGE_OBJ)
	$(CC) $(SFLAGS) $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LFLAGS)

denise:		$(DENISE_OBJ)
	$(CC) $(SFLAGS) $(DENISE_OBJ) -o $(EXEC)/denise $(LFLAGS)

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	rm $(EXEC)/denise 

