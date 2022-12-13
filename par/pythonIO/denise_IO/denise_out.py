# DENISE_out functions
# 
# Daniel Koehn
# Kiel, 24.06.2019

# Import Libraries 
# ----------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LightSource, Normalize
from matplotlib.pyplot import gca
from pylab import rcParams
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter

# Define dictionary for DENISE FD modelling/FWI parameters
para = {
    "NX": 1,
    "NY": 1,
    "DH": 0.1
}

def calc_max_freq(vp,vs,para):
	
	# define gridpoints per minimum wavelength for Taylor and Holberg operators 
	gridpoints_per_wavelength = np.matrix('23. 8. 6. 5. 5. 4. ; 49.7 8.32 4.77 3.69 3.19 2.91 ; 22.2 5.65 3.74 3.11 2.80 2.62 ; 15.8 4.8 3.39 2.90 2.65 2.51 ; 9.16 3.47 2.91 2.61 2.45 2.36')

	# estimate minimum s-wave velocity != 0 in the model 
	minvs = np.min(vs[np.nonzero(vs)])
	print("Minimum Vs: ", minvs, " m/s")

	# estimate minimum p-wave velocity != 0 in the model 
	minvp = np.min(vp[np.nonzero(vp)])
	print("Minimum Vp: ", minvp, " m/s")

	# estimate minimum seismic velocity in model
	minvel = minvs
	if(minvp < minvs):
		minvel = minvp

	# calculate maximum frequency of the source wavelet based on grid dispersion criterion
	print("no of gridpoints per minimum wavelength = ", gridpoints_per_wavelength[para["max_relative_error"],(int)(para["FD_ORDER"]/2)-1])
	fmax = minvel / (gridpoints_per_wavelength[para["max_relative_error"],(int)(para["FD_ORDER"]/2)-1] * para["DH"])
	print("maximum source wavelet frequency = ", fmax, " Hz")

	return fmax
	
def check_domain_decomp(para):

	# check if domain decomposition in x-direction is correct
	decomp_x = para["NX"] % para["NPROCX"]
	
	if decomp_x == 0 :
		print("Domain decomposition in x-direction is correct NX % NPROCX = 0")
	if decomp_x != 0 :
		print("Domain decomposition in x-direction is not correct NX % NPROCX = ", decomp_x)

	# check if domain decomposition in y-direction is correct
	decomp_y = para["NY"] % para["NPROCY"]
	
	if decomp_y == 0 :
		print("Domain decomposition in y-direction is correct NY % NPROCY = 0")
	if decomp_y != 0 :
		print("Domain decomposition in y-direction is not correct NY % NPROCY = ", decomp_y)

	return

def check_src_rec_pml(x,y,para,sws):

	# check if sources & receivers are located in the computational domain and not the PMLs
	flag = 0    # flag source and receiver positions in PML boundary frame
	n = len(x)  # number of sources or receivers
	for i in range(n):
		
		# left boundary 
		if x[i] <= para["FW"] * para["DH"] and x[i] >= para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located in the left PML boundary frame")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located in the left PML boundary frame")
				flag+=1

		if x[i] < para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located outside FD grid (left boundary)")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located outside FD grid (left boundary)")
				flag+=1
				
		# right boundary 
		if x[i] >= (para["NX"] - para["FW"]) * para["DH"] and x[i] <= para["NX"] * para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located in the right PML boundary frame")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located in the right PML boundary frame")
				flag+=1

		if x[i] > para["NX"] * para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located outside FD grid (right boundary)")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located outside FD grid (right boundary)")
				flag+=1


		# top boundary 
		if para["FREE_SURF"] != 1 and y[i] <= para["FW"] * para["DH"] and y[i] >= para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located in the top PML boundary frame")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located in the top PML boundary frame")
				flag+=1

		if y[i] < para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located outside FD grid (top boundary)")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located outside FD grid (top boundary)")
				flag+=1


		# bottom boundary 
		if y[i] >= (para["NY"] - para["FW"]) * para["DH"] and y[i] <= para["NY"] * para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located in the bottom PML boundary frame")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located in the bottom PML boundary frame")
				flag+=1

		if y[i] > para["NY"] * para["DH"]:
			if(sws==1):
				print("receiver # " + str(i) + " is located outside FD grid (bottom boundary)")
				flag+=1				

			if(sws==2):
				print("source # " + str(i) + " is located outside FD grid (bottom boundary)")
				flag+=1


	if flag == 0 and sws == 1:
		print("No receivers in the PML boundary frame - test passed")

	if flag == 0 and sws == 2:
		print("No sources in the PML boundary frame - test passed")

	return

def check_stability(vp,vs,para):
	
	# define FD operator weights for Taylor and Holberg operators: 
	
	# Taylor coefficients
	if(para["max_relative_error"] == 0):
		hc = np.matrix('1.0 0.0 0.0 0.0 0.0 0.0 ; 9.0/8.0 -1.0/24.0 0.0 0.0 0.0 0.0 ; 75.0/64.0 -25.0/384.0 3.0/640.0 0.0 0.0 0.0 ; 1225.0/1024.0 -245.0/3072.0 49.0/5120.0 -5.0/7168.0 0.0 0.0 ; 19845.0/16384.0 -735.0/8192.0 567.0/40960.0 -405.0/229376.0 35.0/294912.0 0.0 ; 160083.0/131072.0 -12705.0/131072.0 22869.0/1310720.0 -5445.0/1835008.0 847.0/2359296.0 -63.0/2883584.0')	

	# Holberg coefficients (E = 0.1%)
	if(para["max_relative_error"] == 1):
		hc = np.matrix('1.0010 0.0 0.0 0.0 0.0 0.0 ; 1.1382 -0.046414 0.0 0.0 0.0 0.0 ; 1.1965 -0.078804 0.0081781 0.0 0.0 0.0 ; 1.2257 -0.099537 0.018063 -0.0026274 0.0 0.0 ; 1.2415 -0.11231 0.026191 -0.0064682 0.001191 0.0 ; 1.2508 -0.12034 0.032131 -0.010142 0.0029857 -0.00066667')		
	
	# Holberg coefficients (E = 0.5%)
	if(para["max_relative_error"] == 2):
		hc = np.matrix('1.0050 0.0 0.0 0.0 0.0 0.0 ; 1.1534 -0.052806 0.0 0.0 0.0 0.0 ; 1.2111 -0.088313 0.011768 0.0 0.0 0.0 ; 1.2367 -0.10815 0.023113 -0.0046905 0.0 0.0 ; 1.2496 -0.11921 0.031130 -0.0093272 0.0025161 0.0 ; 1.2568 -0.12573 0.036423 -0.013132 0.0047484 -0.0015979')
	
	# Holberg coefficients (E = 1.0%)
	if(para["max_relative_error"] == 3):
		hc = np.matrix('1.0100 0.0 0.0 0.0 0.0 0.0 ; 1.1640 -0.057991 0.0 0.0 0.0 0.0 ; 1.2192 -0.094070 0.014608 0.0 0.0 0.0 ; 1.2422 -0.11269 0.026140 -0.0064054 0.0 0.0 ; 1.2534 -0.12257 0.033755 -0.011081 0.0036784 0.0 ; 1.2596 -0.12825 0.038550 -0.014763 0.0058619 -0.0024538')
	
	# Holberg coefficients (E = 3.0%)
	if(para["max_relative_error"] == 4):
		hc = np.matrix('1.0300 0.0 0.0 0.0 0.0 0.0 ; 1.1876 -0.072518 0.0 0.0 0.0 0.0 ; 1.2341 -0.10569 0.022589 0.0 0.0 0.0 ; 1.2516 -0.12085 0.032236 -0.011459 0.0 0.0 ; 1.2596 -0.12829 0.038533 -0.014681 0.0072580 0.0 ; 1.2640 -0.13239 0.042217 -0.017803 0.0081959 -0.0051848')
	
	# estimate maximum s-wave velocity != 0 in the model 
	maxvs = np.max(vs[np.nonzero(vs)])
	print("Maximum Vs: ", maxvs, " m/s")

	# estimate maximum p-wave velocity != 0 in the model 
	maxvp = np.max(vp[np.nonzero(vp)])
	print("Maximum Vp: ", maxvp, " m/s")

	# estimate maximum seismic velocity in model
	maxvel = maxvp
	if(maxvp < maxvs):
		maxvel = maxvs
	
	# calculate dt according to CFL criterion	
	fdcoeff = (int)(para["FD_ORDER"] / 2) - 1 
	gamma = np.sum(np.abs(hc[fdcoeff,:]))
	dt = para["DH"] / (np.sqrt(2.) * gamma * maxvel)
	
	print("According to the Courant-Friedrichs-Lewy (CFL) criterion")
	print("the maximum time step is DT = " + '{:.2e}'.format(dt) + " s"  )
	#print("gamma = ", gamma)
	
	return dt

def check_steplength(nshot,para):
	
	# check steplength parameters
	if(para["TESTSHOT_END"] > nshot):
		print("TESTSHOT_END = " + str(para["TESTSHOT_END"]) + " > " + "nshot = " + str(nshot) + "\n")
		print("Choose TESTSHOT_END <= nshot")
	
	elif(para["TESTSHOT_START"] < 1):
		print("TESTSHOT_START = " + str(para["TESTSHOT_START"]) + " < " + " 1 \n")
		print("Choose TESTSHOT_START >= 1")

	else:
		print("Step length check - passed")

	return
	
def do_plot(n, model, cm, an, title, vpmin, vpmax, x, y, font):
    
    ax=plt.subplot(3, 1, n)
    #ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    #ax.set_yticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5])
    
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    #plt.rc('text', usetex=True)
    #rc('text', usetex=True)
    
    # plt.pcolor(x, y, vp, cmap=cm, vmin=vpmin)
    plt.imshow(model, cmap=cm, interpolation='none', extent=[np.min(x),np.max(x),np.min(y),np.max(y)], vmin=vpmin, vmax=vpmax)
    #a = gca()
    #a.set_xticklabels(a.get_xticks(), font)
    #a.set_yticklabels(a.get_yticks(), font)
    plt.axis('scaled')
    plt.ylabel('Depth [km]', fontdict=font)
    if n==3:
        plt.xlabel('Distance [km]', fontdict=font)
    plt.gca().invert_yaxis()
    cbar=plt.colorbar(aspect=8, pad=0.02)
    cbar.set_label(title, fontdict=font, labelpad=10)
    plt.text(0.1, 0.32,an,fontdict=font,color='white')
    plt.tight_layout()

def plot_model(vp,vs,rho,x,y,cmap,vpmin,vpmax,vsmin,vsmax,rhomin,rhomax):

	FSize = 20
	font = {'color':  'black',
			'weight': 'normal',
			'size': FSize}
	mpl.rc('xtick', labelsize=FSize) 
	mpl.rc('ytick', labelsize=FSize) 
	rcParams['figure.figsize'] = 12, 11

	plt.close('all')
	plt.figure()
	do_plot(1, vp, cmap, '(a)', r"$\rm{V_p [m/s]}$", vpmin, vpmax, x, y, font)
	do_plot(2, vs, cmap, '(b)', r"$\rm{V_s [m/s]}$", vsmin, vsmax, x, y, font)
	do_plot(3, rho, cmap, '(c)', r"$\rm{\rho [kg/m^3]}$", rhomin, rhomax, x, y, font)
	
	#plt.savefig('test.png', format='png', dpi=100)
	plt.savefig('model.pdf', bbox_inches='tight', format='pdf')
	plt.show()
	
	return
	
def plot_acq(vp,xrec,yrec,xsrc,ysrc,x,y,cmap,vpmin,vpmax):

	FSize = 20
	font = {'color':  'black',
			'weight': 'normal',
			'size': FSize}
	mpl.rc('xtick', labelsize=FSize) 
	mpl.rc('ytick', labelsize=FSize)

	plt.figure(figsize=(15, 4))

	plt.imshow(vp, cmap, interpolation='none', extent=[np.min(x),np.max(x),np.min(y),np.max(y)])
	plt.plot(xrec,yrec,'cv',markersize=5)
	plt.plot(xsrc,ysrc,'r*',markersize=10)

	#a = gca()
	#a.set_xticklabels(a.get_xticks(), font)
	#a.set_yticklabels(a.get_yticks(), font)

	plt.ylabel('Depth [m]', fontdict=font)
	plt.xlabel('Distance [m]', fontdict=font)
	plt.gca().invert_yaxis()
	plt.axis('tight')
	cbar=plt.colorbar(aspect=8, pad=0.01)
	cbar.set_label(r"$\rm{V_p [m/s]}$", fontdict=font, labelpad=10)
	plt.tight_layout()
	plt.savefig('Marmousi_model_acq.pdf', bbox_inches='tight', format='pdf')
	plt.show()
	
	return
	
def write_denise_para(para):

	# create and open DENISE parameter file 
	fp = open(para["filename"], mode='w')
	
	# write FD modelling/FWI parameters
	fp.write("#-----------------------------------------------------------------\n")
	fp.write("# PARAMETER FILE FOR DENISE BLACK-EDITION\n")
	fp.write("#-----------------------------------------------------------------\n")
	fp.write("# description:\n")
	fp.write("# description/name of the model: " + para["descr"] + "\n")
	fp.write("#\n")
	fp.write("# ------------------ DENISE Mode ---------------------------------\n")
	fp.write("# Operation mode:\n")
	fp.write("(forward_modelling_only=0;FWI=1;RTM=2)_(MODE) = " + str(para["MODE"]) + "\n")
	fp.write("#\n")
	fp.write("# ---------------- DENISE Physics -----------------------------\n")
	fp.write("(2D-PSV=1;2D-AC=2;2D-VTI=3;2D-TTI=4;2D-SH=5)_(PHYSICS) = " + str(para["PHYSICS"]) + "\n")
	fp.write("#\n")
	fp.write("#-------------- Domain Decomposition -----------------------------\n")
	fp.write("number_of_processors_in_x-direction_(NPROCX) = " + str(para["NPROCX"]) +  "\n")
	fp.write("number_of_processors_in_y-direction_(NPROCY) = " + str(para["NPROCY"]) +  "\n")
	fp.write("#\n")
	fp.write("#-------------------- FD order -----------------------------------\n")
	fp.write("# Order of ssg FD coefficients (values: 2, 4, ..., 12)\n")
	fp.write("FD_ORDER = " + str(para["FD_ORDER"]) +  "\n")
	fp.write("# Maximum relative group velocity error E\n")
	fp.write("# (minimum number of grid points per shortest wavelength is defined by FD_ORDER and E)\n")
	fp.write("# values: 0 = Taylor coefficients\n")
	fp.write("# 1 = Holberg coeff.: E = 0.1 %\n")
	fp.write("# 2 = E = 0.5 %\n")
	fp.write("# 3 = E = 1.0 %\n")
	fp.write("# 4 = E = 3.0 %\n")
	fp.write("max_relative_error = " + str(para["max_relative_error"]) + "\n")
	fp.write("#-------------------- 2-D Grid -----------------------------------\n")
	fp.write("number_of_gridpoints_in_x-direction_(NX) = " + str(para["NX"]) + "\n")
	fp.write("number_of_gridpoints_in_y-direction_(NY) = " + str(para["NY"]) + "\n")
	fp.write("distance_between_gridpoints(in_m)_(DH) = " + str(para["DH"]) + "\n")
	fp.write("#\n")
	fp.write("# Note that y denotes the depth ! Only at y = 0 a free surface boundary condition is applied \n")
	fp.write("#\n")
	fp.write("#-------------------Time Stepping -------------------------------\n")
	fp.write("time_of_wave_propagation_(in_sec)_(TIME) = " + str(para["TIME"]) + "\n")
	fp.write("timestep_(in_seconds)_(DT) = " + '{:.1e}'.format(para["DT"]) + "\n")
	fp.write("#\n")
	fp.write("#--------------------Source---------------------------------------\n")
	fp.write("# Shape_of_source-signal:\n")
	fp.write("(ricker=1;fumue=2;from_SOURCE_FILE=3;SIN**3=4;Gaussian_deriv=5;Spike=6;Klauder=7)_(QUELLART) = " + str(para["QUELLART"]) + "\n")
	fp.write("SIGNAL_FILE = " + para["SIGNAL_FILE"] + "\n")
	fp.write("duration_of_Klauder_wavelet_(in_seconds)_(TS) = " + str(para["TS"]) + "\n")
	fp.write("read_source_positions_from_SOURCE_FILE_(yes=1)_(SRCREC) = 1\n")
	fp.write("SOURCE_FILE = " + para["SOURCE_FILE"] + "\n")
	fp.write("run_multiple_shots_defined_in_SOURCE_FILE_(yes=1)_(RUN_MULTIPLE_SHOTS) = " + str(para["RUN_MULTIPLE_SHOTS"]) + "\n")
	fp.write("corner_frequency_of_highpass_filtered_spike_(FC_SPIKE_1) = " + str(para["FC_SPIKE_1"]) + "\n")
	fp.write("corner_frequency_of_lowpass_filtered_spike_(FC_SPIKE_2) = " + str(para["FC_SPIKE_2"]) + "\n")
	fp.write("order_of_Butterworth_filter_(ORDER_SPIKE) = " + str(para["ORDER_SPIKE"]) + "\n")
	fp.write("write_source_wavelet_(yes=1)_(WRITE_STF) = " + str(para["WRITE_STF"]) + "\n")
	fp.write("#\n")
	fp.write("#--------------------- Model -------------------------------------\n")
	fp.write("read_model_parameters_from_MFILE(yes=1)(READMOD) = 1 \n")
	fp.write("MFILE = " + para["MFILE"] + "\n")
	fp.write("write_model_files_(yes=1)_(WRITEMOD) = 1 \n")
	fp.write("#\n")
	fp.write("#---------------------Q-approximation-----------------------------\n")
	fp.write("Number_of_relaxation_mechanisms_(L) = " + str(para["L"]) + "\n")
	fp.write("L_Relaxation_frequencies_(FL) = " + str(para["FL"]) + "\n")
	fp.write("Tau_(TAU) = 1.0\n")
	fp.write("#\n")
	fp.write("#----------------------Free Surface-------------------------------\n")
	fp.write("free_surface_(yes=1)(FREE_SURF) = " + str(para["FREE_SURF"]) + "\n")
	fp.write("#\n")
	fp.write("#--------------------PML Boundary---------------------------\n")
	fp.write("width_of_absorbing_frame_(in_gridpoints)_(No<=0)_(FW) = " + str(para["FW"]) + "\n")
	fp.write("Damping_velocity_in_CPML_(in_m/s)_(DAMPING) = " + str(para["DAMPING"]) + "\n")
	fp.write("Frequency_within_the_PML_(Hz)_(FPML) = " + str(para["FPML"]) + "\n")
	fp.write("npower = " + str(para["npower"]) + "\n")
	fp.write("k_max_PML = " + str(para["k_max_PML"]) + "\n")
	fp.write("# apply_periodic_boundary_condition_at_edges_(BOUNDARY):\n")
	fp.write("(no=0)_(left_and_right=1) = 0\n")
	fp.write("#\n")
	fp.write("#----------------------Snapshots----------------------------------\n")
	fp.write("output_of_snapshots_(SNAP)(yes>0) = " + str(para["SNAP"]) + "\n")
	fp.write("# output of particle velocities: SNAP=1\n")
	fp.write("# output of pressure field: SNAP=2\n")
	fp.write("# output of curl and divergence energy: SNAP=3\n")
	fp.write("# output of both particle velocities and energy : SNAP=4\n")
	fp.write("write_snapshots_for_shot_no_(SNAP_SHOT) = " + str(para["SNAP_SHOT"]) + "\n")
	fp.write("first_snapshot_(in_sec)_(TSNAP1) = " + str(para["TSNAP1"]) + "\n")
	fp.write("last_snapshot_(in_sec)_(TSNAP2) = " + str(para["TSNAP2"]) + "\n")
	fp.write("increment_(in_sec)_(TSNAPINC) = " + str(para["TSNAPINC"]) + "\n")
	fp.write("increment_x-direction_(IDX) = " + str(para["IDX"]) + "\n")
	fp.write("increment_y-direction_(IDY) = " + str(para["IDY"]) + "\n")
	fp.write("data-format_(SNAP_FORMAT)(ASCII(2);BINARY(3)) = 3\n")
	fp.write("basic_filename_(SNAP_FILE) = " + para["SNAP_FILE"] + "\n")
	fp.write("#\n")
	fp.write("#----------------------Receiver-----------------------------------\n")
	fp.write("output_of_seismograms_(SEISMO) = " + str(para["SEISMO"]) + "\n")
	fp.write("# SEISMO=0: no seismograms\n")
	fp.write("# SEISMO=1: particle-velocities\n")
	fp.write("# SEISMO=2: pressure (hydrophones)\n")
	fp.write("# SEISMO=3: curl and div\n")
	fp.write("# SEISMO=4: everything\n")
	fp.write("read_receiver_positions_from_file_(single_file=1/multiple_files=2)_(READREC) = " + str(para["READREC"]) + "\n")
	fp.write("REC_FILE = " + para["REC_FILE"] + "\n")
	fp.write("reference_point_for_receiver_coordinate_system_(REFREC) = 0.0 , 0.0\n")
	fp.write("#\n")
	fp.write("#-------------------- Towed streamer -------------------------------\n")
	fp.write("# parameters for towed streamer acquisition\n")
	fp.write("The_first_(N_STREAMER)_receivers_in_REC_FILE_belong_to_streamer = 0\n")
	fp.write("Cable_increment_per_shot_(REC_INCR_X) = 80\n")
	fp.write("Cable_increment_per_shot_(REC_INCR_Y) = 0\n")
	fp.write("#\n")
	fp.write("#-------------------- Seismograms --------------------------------\n")
	fp.write("samplingrate_(in_timesteps!)_(NDT) = " + str(para["NDT"]) + "\n")
	fp.write("data-format_(SU(1);ASCII(2);BINARY(3)) = 1\n")
	fp.write("# output files for seismograms\n")
	fp.write("# particle velocities (if SEISMO=1 or SEISMO=4)\n")
	fp.write("filename_for_Vx_(SEIS_FILE_VX) = " + para["SEIS_FILE_VX"] + "\n")
	fp.write("filename_for_Vy_(SEIS_FILE_VY) = " + para["SEIS_FILE_VY"] + "\n")
	fp.write("# curl and div of wavefield (if SEISMO=3 or SEISMO=4)\n")
	fp.write("filename_for_curl_(SEIS_FILE_CURL) = " + para["SEIS_FILE_CURL"] + "\n")
	fp.write("filename_for_div_(SEIS_FILE_DIV) = " + para["SEIS_FILE_DIV"] + "\n")
	fp.write("# pressure wavefield (hydrophones) (if SEISMO=2 or SEISMO=4)\n")
	fp.write("filename_for_pressure_(SEIS_FILE_P) = " + para["SEIS_FILE_P"] + "\n")
	fp.write("#\n")
	fp.write("# ----------------------------------------------------------------\n")
	fp.write("# each PE is printing log-information to LOG_FILE.MYID\n")
	fp.write("log-file_for_information_about_progress_of_program_(LOG_FILE) = " + para["LOG_FILE"] + "\n")
	fp.write("info_of_processing_element_zero_to_stdout_(yes=1/no=0)_(LOG) = 1\n")
	fp.write("# ----------------------------------------------------------------\n")
	fp.write("# DENISE FWI specific parameters\n")
	fp.write("number_of_TDFWI_iterations_(ITERMAX) = " + str(para["ITERMAX"]) + "\n")
	fp.write("output_of_gradient_(JACOBIAN) = " + para["JACOBIAN"] + "\n")
	fp.write("seismograms_of_measured_data_(DATA_DIR) = " + para["DATA_DIR"] + "\n")
	fp.write("cosine_taper_(yes=1/no=0)_(TAPER) = 0\n")
	fp.write("taper_length_(in_rec_numbers)_(TAPERLENGTH) = " + str(para["TAPERLENGTH"]) + "\n")
	fp.write("gradient_taper_geometry_(GRADT1,GRADT2,GRADT3,GRADT4) = " + str(para["GRADT1"]) + ", " + str(para["GRADT2"]) + ", " + str(para["GRADT3"]) + ", " + str(para["GRADT4"]) + "\n")
	fp.write("type_of_material_parameters_to_invert_(Vp,Vs,rho=1/Zp,Zs,rho=2/lam,mu,rho=3)_(INVMAT1) = " + str(para["INVMAT1"]) + "\n")
	fp.write("gradient_formulation_(GRAD_FORM) = " + str(para["GRAD_FORM"]) + "\n")
	fp.write("adjoint_source_type_(x-y_components=1/y_comp=2/x_comp=3/p_comp=4/x-p_comp=5/y-p_comp=6/x-y-p_comp=7)_(QUELLTYPB) = " + str(para["QUELLTYPB"]) + "\n")
	fp.write("testshots_for_step_length_estimation_(TESTSHOT_START,TESTSHOT_END,TESTSHOT_INCR) = " + str(para["TESTSHOT_START"]) + ", " + str(para["TESTSHOT_END"]) + ", " + str(para["TESTSHOT_INCR"]) + "\n")
	fp.write("#\n")
	fp.write("# ----- Definition of gradient taper geometry ----- #\n")
	fp.write("# Vertical taper\n")
	fp.write("apply_vertical_taper_(yes=1)_(SWS_TAPER_GRAD_VERT) = " + str(para["SWS_TAPER_GRAD_VERT"]) + "\n")
	fp.write("# Horizontal taper\n")
	fp.write("apply_horizontal_taper_(yes=1)_(SWS_TAPER_GRAD_HOR) = " + str(para["SWS_TAPER_GRAD_HOR"]) + "\n")
	fp.write("exponent_of_depth_scaling_for_preconditioning_(EXP_TAPER_GRAD_HOR) = " + str(para["EXP_TAPER_GRAD_HOR"]) + "\n")
	fp.write("# Circular taper around all sources (not at receiver positions)\n")
	fp.write("apply_cylindrical_taper_(yes=1)_(SWS_TAPER_GRAD_SOURCES) = " + str(para["SWS_TAPER_GRAD_SOURCES"]) + "\n")
	fp.write("apply_cylindrical_taper_per_shot_(yes=1)_(SWS_TAPER_CIRCULAR_PER_SHOT) = " + str(para["SWS_TAPER_CIRCULAR_PER_SHOT"]) + "\n")
	fp.write("(1=error_function,2=log_function)_(SRTSHAPE) = " + str(para["SRTSHAPE"]) + "\n")
	fp.write("radius_in_m_(SRTRADIUS) = " + str(para["SRTRADIUS"]) + "\n")
	fp.write("# --> minimum for SRTRADIUS is 5x5 gridpoints\n")
	fp.write("filtsize_in_gridpoints_(FILTSIZE) = 1\n")
	fp.write("read_taper_from_file_(yes=1)_(SWS_TAPER_FILE) = " + str(para["SWS_TAPER_FILE"]) + "\n")
	fp.write("taper_file_basename_(TFILE) = " + para["TFILE"] + "\n")
	fp.write("#\n")
	fp.write("# ----- Output of inverted models ----- #\n")
	fp.write("write_inverted_model_after_each_iteration_(yes=1)_(INV_MOD_OUT) = " + str(para["INV_MOD_OUT"]) + "\n")
	fp.write("output_of_models_(INV_MODELFILE) = " + str(para["INV_MODELFILE"]) + "\n")
	fp.write("#\n")
	fp.write("# ----- Upper and lower limits for model parameters ----- #\n")
	fp.write("upper_limit_for_vp/lambda_(VPUPPERLIM) = " + str(para["VPUPPERLIM"]) + "\n")
	fp.write("lower_limit_for_vp/lambda_(VPLOWERLIM) = " + str(para["VPLOWERLIM"]) + "\n")
	fp.write("upper_limit_for_vs/mu_(VSUPPERLIM) = " + str(para["VSUPPERLIM"]) + "\n")
	fp.write("lower_limit_for_vs/mu_(VSLOWERLIM) = " + str(para["VSLOWERLIM"]) + "\n")
	fp.write("upper_limit_for_rho_(RHOUPPERLIM) = " + str(para["RHOUPPERLIM"]) + "\n")
	fp.write("lower_limit_for_rho_(RHOLOWERLIM) = " + str(para["RHOLOWERLIM"]) + "\n")
	fp.write("upper_limit_for_qs_(QSUPPERLIM) = " + str(para["QSUPPERLIM"]) + "\n")
	fp.write("lower_limit_for_qs_(QSLOWERLIM) = " + str(para["QSLOWERLIM"]) + "\n")
	fp.write("#\n")
	fp.write("# ----- Optimization-Method ------ #\n")
	fp.write("gradient_method_(PCG=1/LBFGS=2)_(GRAD_METHOD) = " + str(para["GRAD_METHOD"]) + "\n")
	fp.write("PCG_BETA_(Fletcher_Reeves=1/Polak_Ribiere=2/Hestenes_Stiefel=3/Dai_Yuan=4) = " + str(para["PCG_BETA"]) + "\n")
	fp.write("save_(NLBFGS)_updates_during_LBFGS_optimization = " + str(para["NLBFGS"]) + "\n")
	fp.write("#\n")
	fp.write("#----- Definition of smoothing the models vp and vs ----- #\n")
	fp.write("apply_spatial_filtering_(1=yes)_(MODEL_FILTER) = 0\n")
	fp.write("filter_length_in_gridpoints_(FILT_SIZE) = 5\n")
	fp.write("#\n")
	fp.write("#----- Reduce size of inversion grid ------#\n")
	fp.write("use_only_every_DTINV_time_sample_for_gradient_calculation_(DTINV) = " + str(para["DTINV"]) + "\n")
	fp.write("#\n")
	fp.write("#----- Step length estimation ------# \n")
	fp.write("maximum_model_change_of_maximum_model_value_(EPS_SCALE) = " + str(para["EPS_SCALE"]) + "\n")
	fp.write("maximum_number_of_attemps_to_find_a_step_length_(STEPMAX) = " + str(para["STEPMAX"]) + "\n")
	fp.write("SCALEFAC = " + str(para["SCALEFAC"]) + "\n")
	fp.write("#\n")
	fp.write("#----- Trace killing -----#\n")
	fp.write("apply_trace_killing_(yes=1)_(TRKILL) = " + str(para["TRKILL"]) + "\n")
	fp.write("TRKILL_FILE = " + para["TRKILL_FILE"] + "\n")
	fp.write("#\n")
	fp.write("#----- Time damping -----#\n")
	fp.write("files_with_picked_times_(PICKS_FILE) = " + para["PICKS_FILE"] + "\n")
	fp.write("#\n")
	fp.write("#----- FWI LOG FILE -----#\n")
	fp.write("log_file_for_misfit_evolution_(MISFIT_LOG_FILE) = " + para["MISFIT_LOG_FILE"] + "\n")
	fp.write("#\n")
	fp.write("# ----- Minimum number of iteration per frequency ----- #\n")
	fp.write("MIN_ITER = 0\n")
	fp.write("#\n")
	fp.write("# ----- Definition of smoothing the Jacobians with 2D-Gaussian ----- #\n")
	fp.write("apply_spatial_filtering_(yes=1)_(GRAD_FILTER) = 0\n")
	fp.write("filter_length_in_gridpoints_(FILT_SIZE_GRAD) = 10\n")
	fp.write("#\n")
	fp.write("# ----- FWT double-difference time-lapse mode ---------------------------- #\n")
	fp.write("activate_time_lapse_mode_(yes=1)_(TIMELAPSE) = 0\n")
	fp.write("# if TIMELAPSE == 1, DATA_DIR should be the directory containing the data differences\n")
	fp.write("# between time t0 and t1 \n")
	fp.write("seismograms_of_synthetic_data_at_t0_(DATA_DIR_T0) = su/CAES_spike_time_0/DENISE_CAES\n")
	fp.write("#\n")
	fp.write("# ----- Reverse Time Modelling ------------------------- #\n")
	fp.write("apply_reverse_time_modelling_(yes=1)_(RTMOD) = 0\n")
	fp.write("#\n")
	fp.write("# ----- Gravity Modelling/Inversion -----#\n")
	fp.write("# no gravity modelling and inversion: GRAVITY=0\n")
	fp.write("# activate only gravity modelling: GRAVITY=1\n")
	fp.write("# activate gravity modelling and inversion: GRAVITY=2\n")
	fp.write("active_gravity_modelling_(=1)_and_inversion_(=2)_(GRAVITY) = 0\n")
	fp.write("# boundaries in x-direction in gridpoints\n")
	fp.write("boundary_gridpoints_in_x_(NGRAVB) = 500\n")
	fp.write("# boundaries in z-direction in meter\n")
	fp.write("boundary_meter_in_z_(NZGRAV) = 200000\n")
	fp.write("# model and invert gravity data: GRAV_TYPE=1\n")
	fp.write("# model and invert gravity gradient data: GRAV_TYPE=2\n")
	fp.write("use_of_gravity_(=1)_or_gravity_gradient_(=2)_data_(GRAV_TYPE) = 1\n")
	fp.write("# use initial density model (=1) ,only reasonable for inversion, as background density or self-defined model (=2)\n")
	fp.write("chosen_background_density_model_(BACK_DENSITY) = 2\n")
	fp.write("# if BACK_DENSITY = 2, define your model \n")
	fp.write("filename_for_background_density_(DFILE) = gravity/background_density.rho\n")
	fp.write("#\n")
	fp.write("# ----- RTM parameters ---------------------------- #\n")
	fp.write("output_of_RTM_result_for_each_shot_(yes=1)_(RTM_SHOT) = 0\n")
	fp.write("#\n")
	
	# close DENISE parameter file
	fp.close()

	return

def write_denise_workflow_header(para):

	# create and open DENISE workflow file 
	fp = open(para["filename_workflow"], mode='w')

	# write header to DENISE workflow file
	fp.write("PRO \t TIME_FILT \t FC_low \t FC_high \t ORDER \t TIME_WIN \t GAMMA \t TWIN- \t TWIN+ \t INV_VP_ITER \t INV_VS_ITER \t INV_RHO_ITER \t INV_QS_ITER \t SPATFILTER \t WD_DAMP \t WD_DAMP1 \t EPRECOND \t LNORM \t ROWI \t STF_INV \t OFFSETC_STF \t EPS_STF \t NORMALIZE \t OFFSET_MUTE \t OFFSETC \t SCALERHO \t SCALEQS \t ENV \t GAMMA_GRAV \t N_ORDER \n")

	# close DENISE workflow file
	fp.close()
	
	return	


def write_denise_workflow(para):

	# open DENISE workflow file 
	fp = open(para["filename_workflow"], mode='a')

	# write header to DENISE workflow file
	fp.write(str(para["PRO"]) + "\t" + str(para["TIME_FILT"]) + "\t" + str(para["FC_LOW"]) + "\t" + str(para["FC_HIGH"]) + "\t" + str(para["ORDER"]) + "\t" + str(para["TIME_WIN"]) + "\t" + str(para["GAMMA"]) + "\t" + str(para["TWIN-"]) + "\t" + str(para["TWIN+"]) + "\t" + str(para["INV_VP_ITER"]) + "\t" + str(para["INV_VS_ITER"]) + "\t" + str(para["INV_RHO_ITER"]) + "\t" + str(para["INV_QS_ITER"]) + "\t" + str(para["SPATFILTER"]) + "\t" + str(para["WD_DAMP"]) + "\t" + str(para["WD_DAMP1"]) + "\t" + str(para["EPRECOND"]) + "\t" + str(para["LNORM"]) + "\t" + str(para["ROWI"]) + "\t" + str(para["STF"]) + "\t" + str(para["OFFSETC_STF"]) + "\t" + str(para["EPS_STF"]) + "\t" + "0" + "\t" + str(para["OFFSET_MUTE"]) + "\t" + str(para["OFFSETC"]) + "\t" + str(para["SCALERHO"]) + "\t" + str(para["SCALEQS"]) + "\t" + str(para["ENV"]) + "\t" + "0" + "\t" + str(para["N_ORDER"]) + "\n")

	# close DENISE workflow file
	fp.close()
	
	return	


