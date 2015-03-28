# Makefile for fdveps

#--------------------------------------------------------
# edit here:

# source code for model generation

#MODEL = hh.c
MODEL = model_ainos_visc.c
MODEL_EL = model_ainos.c
EXEC= ../bin

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4 

# On CRAY T3E
# CC=cc

# On Linux RZ-cluster with Intel-MPI
CC=mpiicc
LFLAGS=-lm -lcseife -lfftw3 -lstdc++
CFLAGS=-O3 
SFLAGS=-L./../libcseife -L$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/lib
IFLAGS=-I./../libcseife -I$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
#CC=mpicc
#LFLAGS=-lm -lcseife -lfftw3 -lstdc++
#CFLAGS=-O3 -fno-stack-protector
#SFLAGS=-L./../libcseife
#IFLAGS=-I./../libcseife

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

AUTO_PICK_SRC = \
        auto_pick.c \
        read_par.c \
        exchange_par.c \
        inseis.c \
        calc_res.c \
        util.c \
        receiver.c \
        splitrec.c \
        stalta.c \
        mer.c \
        timedomain_filt.c \
        time_window.c \
        sources.c \
        write_par.c \
        initproc.c
        
SNAPMERGE_SCR = \
	snapmerge.c \
	merge.c \
	read_par.c \
	readdsk.c \
	writedsk.c \
	util.c 


FDVEPS_UTIL = \
	absorb.c \
	exchange_par.c \
	info.c \
	initproc.c \
	$(MODEL) \
	merge2.c \
	mergemod.c \
	note.c \
	outseis.c \
	read_checkpoint.c\
	readdsk.c \
	read_par.c \
	readmod.c \
    	receiver.c \
	rd_sour.c \
	save_checkpoint.c\
	saveseis.c \
	sources.c \
	splitsrc.c \
	splitrec.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c 


FDVEPS_SRC_SSG= \
	fdveps_ssg.c \
	holbergcoeff.c\
	checkfd_hc.c \
	psource.c \
	av_mue.c\
	av_tau.c\
	comm_ini.c\
	exchange_v.c \
	exchange_s.c \
	snap_ssg.c \
	seismo_ssg.c \
	surface.c \
	update_v_ssg.c \
	update_s_visc_ssg.c \
	matcopy.c \
	$(FDVEPS_UTIL) 

FDVEPS_SRC_RSG= \
	fdveps_rsg.c \
	checkfd_rsg.c \
	psource_rsg.c \
	exchange_rsg.c \
	snap_rsg.c \
	seismo_rsg.c \
	update_v_rsg_2nd.c \
	update_s_rsg_2nd.c \
	matcopy.c \
	$(MODEL_EL) \
	$(FDVEPS_UTIL) 

FDELAST_SRC_RSG_2ND= \
	fdelast_rsg.c \
	checkfd_rsg_elastic_2nd.c \
	psource_rsg.c \
	exchange_rsg.c \
	snap_rsg.c \
	seismo_rsg.c \
	update_v_rsg_2nd.c \
	update_s_elastic_rsg_2nd.c \
	matcopy_elastic.c \
	$(MODEL_EL) \
	readmod_elastic.c \
	$(FDVEPS_UTIL) 

FDELAST_SRC_SSG_4TH= \
	fdelast_ssg_old.c \
	checkfd_ssg_elastic_4th.c \
	psource.c \
	exchange_v.c \
	exchange_s.c \
	snap_ssg_4th.c \
	seismo_ssg_4th.c \
	surface_elastic_4th.c \
	update_v_ssg_4th.c \
	update_s_elastic_ssg_4th.c \
	absorb.c \
	av_mue.c \
	exchange_par.c \
	info.c \
	initproc.c \
	$(MODEL_EL) \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	note.c \
	outseis.c \
	output_source_signal.c \
	readdsk.c \
	read_par.c \
	read_checkpoint.c\
	readmod_elastic.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	save_checkpoint.c\
	sources.c \
	splitsrc.c \
	splitrec.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c 

FDELAST_SRC_SSG_2ND= \
	fdelast_ssg.c \
	checkfd_ssg_elastic_2nd.c \
	psource.c \
	exchange_v.c \
	exchange_s.c \
	snap_ssg_2nd.c \
	seismo_ssg_2nd.c \
	surface_elastic_2nd.c \
	update_v_ssg_2nd.c \
	update_s_elastic_ssg_2nd.c \
	absorb.c \
	av_mue.c \
	exchange_par.c \
	info.c \
	initproc.c \
	$(MODEL_EL) \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	note.c \
	outseis.c \
	output_source_signal.c \
	readdsk.c \
	read_par.c \
	read_checkpoint.c\
	readmod_elastic.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	save_checkpoint.c\
	sources.c \
	splitsrc.c \
	splitrec.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c 


DENISE= \
	denise.c \
	dotp.c \
	stf.c \
	window_cos.c \
	calc_envelope.c \
	calc_hilbert.c \
	calc_mat_change.c \
        calc_mat_change_test.c \
        calc_res.c \
        calc_res_grav.c \
	calc_misfit.c \
        calc_opt_step.c \
	calc_opt_step_test.c \
	calc_energy.c \
	checkfd_ssg_elastic.c \
	checkfd_ssg_visc.c \
	conv_FD.c \
	psource.c \
	hessian_out.c \
	holbergcoeff.c \
        householder.c \
	comm_ini.c \
	exchange_v.c \
	exchange_s.c \
	exchange_L2.c \
	eprecond.c \
	eprecond1.c \
	forward_mod.c \
        grav_grad.c \
        grav_mod.c \
	snap_ssg.c \
	seismo_ssg.c \
	step_length_est.c \
	surface_elastic_2nd.c \
	surface_elastic_PML.c \
	surface_PML.c \
	update_v_ssg.c \
	update_v_PML.c \
	prepare_update_s.c \
	update_s_elastic_ssg.c \
	update_s_elastic_PML.c \
	update_s_visc_PML.c \
	absorb.c \
	av_mue.c \
	av_rho.c \
	av_tau.c \
	median.c \
	exchange_par.c \
        extend_mod.c \
	info.c \
	inseis.c \
	initproc.c \
	interpol.c \
	LBFGS1.c \
	LU_decomp.c \
	smooth2.c \
	smooth_model.c \
	smooth_grad.c \
	$(MODEL) \
	$(MODEL_EL) \
	matcopy.c \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	median_src.c \
	max_grad.c \
	model_freq_out.c \
	note.c \
	norm.c \
	outseis.c \
	outseis_vector.c \
	outseis_glob.c \
	catseis.c \
	output_source_signal.c \
	PCG.c \
	PML_pro.c \
        read_density_glob.c \
	readdsk.c \
        read_grav_pos.c \
	read_par.c \
	read_par_inv.c \
	read_checkpoint.c\
	readmod.c \
	readmod_elastic.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	saveseis_glob.c \
	save_checkpoint.c\
	sources.c \
	solvelin.c \
	spat_filt.c \
	splitsrc.c \
	splitsrc_back.c \
	splitrec.c \
	stalta.c \
	taper.c \
	taper_grad.c \
	taper_grad_shot.c \
	timedomain_filt.c \
	timedomain_filt_vector.c \
	time_window.c \
	time_window_stf.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	wavenumber.c \
	writemod.c \
	write_par.c \
	writedsk.c \
	zero_fdveps.c \
	zero_fdveps_visc.c \
	zero_LBFGS.c \
	zero_hessian.c

DENISE_ES_SRC_SSG= \
	denise_es.c \
	calc_mat_change.c \
        calc_mat_change_test.c \
        calc_res.c \
        calc_opt_step.c \
	calc_opt_step_test.c \
	checkfd_ssg_elastic.c \
	psource.c \
	holbergcoeff.c\
	comm_ini.c\
	exchange_v.c \
	exchange_s.c \
	exchange_L2.c \
	exchange_par.c \
	exchange_mod_es.c \
	fft.c \
	fft_filt.c \
	snap_ssg.c \
	seismo_ssg.c \
	surface_elastic_2nd.c \
	update_v_ssg.c \
	update_s_elastic_ssg.c \
	absorb.c \
	av_mue.c \
	av_rho.c \
	median.c \
	info.c \
	inseis.c \
	initproc.c \
	interpol.c \
	LU_decomp.c \
	$(MODEL_EL) \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	max_grad.c \
	note.c \
	outseis.c \
	output_source_signal.c \
	readdsk.c \
	read_par.c \
	read_checkpoint.c\
	readmod_elastic.c \
	readmod_elastic_es.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	save_checkpoint.c\
	sources.c \
	solvelin.c \
	spat_filt.c \
	splitsrc.c \
	splitsrc_back.c \
	splitrec.c \
	stalta.c \
	taper.c \
	taper_grad.c \
	time_window.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c \
	zero_fdveps.c
	
FDELAST_SRC_SSG_2ND_PML= \
	fdelast_ssg_PML.c \
	checkfd_ssg_elastic_2nd.c \
	psource_PML.c \
	exchange_v_PML.c \
	exchange_s_PML.c \
	snap_ssg_2nd_PML.c \
	seismo_ssg_2nd_PML.c \
	surface_elastic_2nd.c \
	update_v_ssg_2nd_PML.c \
	update_s_elastic_ssg_2nd_PML.c \
	absorb_PML.c \
	av_mue.c \
	exchange_par_PML.c \
	info.c \
	initproc.c \
	$(MODEL_EL) \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	note.c \
	outseis.c \
	output_source_signal.c \
	readdsk.c \
	read_par.c \
	read_checkpoint.c\
	readmod_elastic.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	save_checkpoint.c\
	sources.c \
	splitsrc.c \
	splitrec.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c 
		
AFDELAST_SRC_SSG_2ND= \
	fdelast_ssg_adapt.c \
	checkfd_ssg_elastic_2nd.c \
	psource.c \
	exchange_v.c \
	exchange_s.c \
	snap_ssg_2nd.c \
	seismo_ssg_2nd.c \
	surface_elastic_2nd.c \
	update_v_ssg_2nd_adapt.c \
	update_s_elastic_ssg_2nd_adapt.c \
	absorb.c \
	av_mue.c \
	exchange_par.c \
	info.c \
	initproc.c \
	interpol.c \
	$(MODEL_EL) \
	matcopy_elastic.c \
	merge2.c \
	mergemod.c \
	note.c \
	outseis.c \
	readdsk.c \
	read_par.c \
	read_checkpoint.c\
	readmod_elastic.c \
    	receiver.c \
	rd_sour.c \
	saveseis.c \
	save_checkpoint.c\
	sources.c \
	splitsrc.c \
	splitrec.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c 

AUTO_PICK_OBJ = $(AUTO_PICK_SRC:%.c=%.o)	
SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)

FDVEPS_OBJ_SSG = $(FDVEPS_SRC_SSG:%.c=%.o)
FDVEPS_OBJ_RSG = $(FDVEPS_SRC_RSG:%.c=%.o)

FDELAST_OBJ_SSG_4TH = $(FDELAST_SRC_SSG_4TH:%.c=%.o)
FDELAST_OBJ_SSG_2ND = $(FDELAST_SRC_SSG_2ND:%.c=%.o)
DENISE_OBJ = $(DENISE:%.c=%.o)
DENISE_ES_OBJ_SSG = $(DENISE_ES_SRC_SSG:%.c=%.o)
FDELAST_OBJ_RSG_2ND = $(FDELAST_SRC_RSG_2ND:%.c=%.o)
AFDELAST_OBJ_SSG_2ND = $(AFDELAST_SRC_SSG_2ND:%.c=%.o)
FDELAST_OBJ_SSG_2ND_PML = $(FDELAST_SRC_SSG_2ND_PML:%.c=%.o)

auto_pick:	$(AUTO_PICK_OBJ)
	$(CC) $(SFLAGS) $(AUTO_PICK_OBJ) -o ../bin/autopick $(LFLAGS)        

snapmerge:	$(SNAPMERGE_OBJ)
	$(CC) $(SFLAGS) $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LFLAGS)

fdveps_ssg:		$(FDVEPS_OBJ_SSG)
	$(CC) $(LFLAGS) $(FDVEPS_OBJ_SSG) -o $(EXEC)/fdveps_ssg

fdveps_rsg:		$(FDVEPS_OBJ_RSG)
	$(CC) $(LFLAGS) $(FDVEPS_OBJ_RSG) -o $(EXEC)/fdveps_rsg

fdelast_ssg_4th:		$(FDELAST_OBJ_SSG_4TH)
	$(CC) $(LFLAGS) $(FDELAST_OBJ_SSG_4TH) -o $(EXEC)/fdelast_ssg_4th

fdelast_ssg_2nd:		$(FDELAST_OBJ_SSG_2ND)
	$(CC) $(LFLAGS) $(FDELAST_OBJ_SSG_2ND) -o $(EXEC)/fdelast_ssg_2nd

denise:		$(DENISE_OBJ)
	$(CC) $(SFLAGS) $(DENISE_OBJ) -o $(EXEC)/denise $(LFLAGS)

fdelast_rsg_2nd:		$(FDELAST_OBJ_RSG_2ND)
	$(CC) $(LFLAGS) $(FDELAST_OBJ_RSG_2ND) -o $(EXEC)/fdelast_rsg_2nd

afdelast_ssg_2nd:		$(AFDELAST_OBJ_SSG_2ND)
	$(CC) $(LFLAGS) $(AFDELAST_OBJ_SSG_2ND) -o $(EXEC)/fdelast_ssg_2nd_adapt	

fdelast_ssg_2nd_pml:		$(FDELAST_OBJ_SSG_2ND_PML)
	$(CC) $(LFLAGS) $(FDELAST_OBJ_SSG_2ND_PML) -o $(EXEC)/fdelast_ssg_2nd_pml

denise_es:		$(DENISE_ES_OBJ_SSG)
	$(CC) $(LFLAGS) $(DENISE_ES_OBJ_SSG) -o $(EXEC)/denise_es

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
