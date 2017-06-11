# Makefile for fdveps

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
#CC=mpiicc
# LFLAGS=-lm -lcseife -lfftw3 -lstdc++
# CFLAGS=-O3 
#LFLAGS=-lm -lcseife -lfftw3 -ftz -lstdc++
#CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
#SFLAGS=-L./../libcseife -L$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/lib
#IFLAGS=-I./../libcseife -I$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/include -I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
CC=mpicc
LFLAGS=-lm -lcseife -lfftw3 -lstdc++
CFLAGS=-O3 -w -fno-stack-protector -D_FORTIFY_SOURCE=0
SFLAGS=-L./../libcseife
IFLAGS=-I./../libcseife -I./../include

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

VPATH = snapmerge:PSV:VTI:TTI:AC:models
        
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
	  alloc_matAC.c \
	  checkfd_acoustic.c \
	  dealloc_AC.c \
	  exchange_p_AC.c \
	  exchange_v_AC.c \
	  FD_AC.c \
	  matcopy_acoustic_AC.c \
	  $(MODEL_AC) \
	  physics_AC.c \
	  psource_AC.c \
	  readmod_AC.c \
	  seismo_AC.c \
	  snap_AC.c \
	  surface_acoustic_PML_AC.c \
	  update_s_acoustic_PML_AC.c \
	  update_v_PML_AC.c \
	  zero_denise_acoustic_AC.c

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
