# Common modules
import numpy as np
import cPickle as pickle
import argparse
import os
import time
import matplotlib.pyplot as plt

# Own modules
import ReadHaloData as RHD
import OtherFunctions as OF

class AnalyseCritPts():
	def __init__(self, model, npart, sigma):
		self.npart = npart
		self.sigma = sigma
		self.model = model
		if model == 'lcdm':
			Dispersemodel = 'LCDM'
		elif model == 'symmA':
			model = 'symm_A'
			Dispersemodel = 'SymmA'
		elif model == 'symmB':
			model = 'symm_B'
			Dispersemodel = 'SymmB'
		elif model == 'symmC':
			model = 'symm_C'
			Dispersemodel = 'SymmC'
		elif model == 'symmD':
			model = 'symm_D'
			Dispersemodel = 'SymmD'
		else:
			Dispersemodel = model

		#self.read_filament_data(Dispersemodel, npart, sigma)
		self.HaloID, self.HaloPos, self.HaloMass, self.VirRadius = RHD.Read_halo_data(model)

	def read_filament_data(self, model, npart, sigma):
		""" Reads filament data from DisPerSE """
		print "Reading filament data"
		cachedir_foldername_extra = model + 'npart'+str(npart) + 'nsig'+str(sigma)
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/FilamentData/' + cachedir_foldername_extra + '/'
		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'Conditions_check.p'
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'
		# Read critical points from pickle and creates 3D coordinates of critical point positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(CP_coordinate_data)
		CP_3DPos = []
		for j in range(len(self.CritPointXpos)):
			CP_3DPos.append(np.column_stack((self.CritPointXpos[j], self.CritPointYpos[j], self.CritPointZpos[j])))
		self.CP_3DPos = np.array(CP_3DPos)
		Maxima_cp = self.CP_type == 3
		Saddle_cp = self.CP_type == 2
		self.CP_3DPos_max = self.CP_3DPos[Maxima_cp]
		self.CP_3DPos_saddle = self.CP_3DPos[Saddle_cp]

	def Unpack_filament_data(self, CP_coord):
		""" Unpacks critial point data from the read filament data module. """
		self.CritPointXpos = CP_coord[0]
		self.CritPointYpos = CP_coord[1]
		self.CritPointZpos = CP_coord[2]
		self.CP_type = np.array(CP_coord[3])
		self.CP_persistent_pair = CP_coord[4] 
		self.Critpts_filamentID = CP_coord[5] 
		self.CP_id_of_connecting_filament = CP_coord[6]
		self.Number_filaments_connecting_to_CP = np.array(CP_coord[7])

	def CP_in_halo(self, HaloCord, CPCord, radius):
		""" Checks whether critical points are within a halo with a given radius """
		norm_axis = 1
		if len(CPCord) == 1:
			norm_axis = 0
		Diff_len = np.linalg.norm(HaloCord - CPCord, axis=norm_axis)**2
		InHaloes = Diff_len <= radius**2
		return InHaloes

	def Check_cps(self):
		""" Checking how many critical points there are within each halo """

		cachedir_NumCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/NumCPInHalo/'
		cachedir_CPIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/CPIDsInHalo/'
		cachedir_MaximaCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/MaximaCPInHalo/'
		cachedir_SaddleCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/SaddleCPInHalo/'
		if not os.path.isdir(cachedir_NumCP):
			os.makedirs(cachedir_NumCP)
		if not os.path.isdir(cachedir_CPIDs):
			os.makedirs(cachedir_CPIDs)
		if not os.path.isdir(cachedir_MaximaCP):
			os.makedirs(cachedir_MaximaCP)
		if not os.path.isdir(cachedir_SaddleCP):
			os.makedirs(cachedir_SaddleCP)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_NumCP = cachedir_NumCP + Common_filename
		cachefile_CPIDs = cachedir_CPIDs + Common_filename
		cachefile_MaximaCP = cachedir_MaximaCP + Common_filename
		cachefile_SaddleCP = cachedir_SaddleCP + Common_filename
		if os.path.isfile(cachefile_NumCP):
			print "Reading CPs in halos of model " + self.model + "  ..."
			NumCP_in_halo = np.load(cachefile_NumCP)
			CP_ids_in_halo = np.load(cachefile_CPIDs)
		else:
			print "Calculating CPs in halos"
			timer = time.time()
			NumCP_in_halo = []
			CP_ids_in_halo = []
			for i in range(len(self.HaloPos)):
				InHalo = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos, self.VirRadius[i])
				OK_cps = np.where(InHalo)[0]
				NumCP_in_halo.append(len(OK_cps))
				CP_ids_in_halo.append(OK_cps.astype(np.int32))
			NumCP_in_halo = np.asarray(NumCP_in_halo)
			CP_ids_in_halo = np.asarray(CP_ids_in_halo)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_NumCP, NumCP_in_halo)
			np.save(cachefile_CPIDs, CP_ids_in_halo)

		if: os.path.isfile(cachefile_MaximaCP):
			print "Reading maxima and saddle cps in halos ..."
			MaximaCP_in_halo = np.load(cachefile_MaximaCP)
			SaddleCP_in_halo = np.load(cachefile_SaddleCP)
		else:
			print "Calculating CPs (maxima and saddle seperate) in halos"
			timer = time.time()
			MaximaCP_in_halo = []
			SaddleCP_in_halo = []
			for i in range(len(self.HaloPos)):
				InHalo_max = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos_max, self.VirRadius[i])
				InHalo_saddle = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos_saddle, self.VirRadius[i])
				OK_cps_m = np.where(InHalo_max)[0]
				OK_cps_s = np.where(InHalo_saddle)[0]
				MaximaCP_in_halo.append(len(OK_cps_m))
				SaddleCP_in_halo.append(len(OK_cps_s))
			MaximaCP_in_halo = np.asarray(MaximaCP_in_halo)
			SaddleCP_in_halo = np.asarray(SaddleCP_in_halo)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_MaximaCP, MaximaCP_in_halo)
			np.save(cachefile_SaddleCP, SaddleCP_in_halo)
		return NumCP_in_halo, CP_ids_in_halo, MaximaCP_in_halo, SaddleCPInHalo

class Plot_results():
	def __init__(self, models, Nsigma, foldername, filetype='png'):
		if filetype == 'png':
			foldername += 'PNG/'
		elif filetype == 'pdf':
			foldername += 'PDF/'
		self.filetype = '.' + filetype
		self.raw_filetype = filetype
		sigma_name_folder = 'Sigma'+str(Nsigma) + '/'
		foldername += sigma_name_folder
		self.Nsigma = Nsigma
		self.results_dir = os.path.join(savefile_directory, foldername)
		if not os.path.isdir(self.results_dir):
			os.makedirs(self.results_dir)

		self.All_legends = []
		self.Symm_legends = []
		self.fofr_legends = []
		self.Plot_colors_symm = ['k', 'orange', 'g', 'r', 'c']
		self.Plot_colors_fofr = ['k', 'purple', 'y', 'b']
		self.ID_counter = 0
		self.Linestyles = ['-', '--', '-.', ':', (0, (5, 10))]
		self.SymmLCDM_ids = np.array([])
		self.FofrLCDM_ids = np.array([])
		self.Get_legends(models)
		self.SymmLCDM_ids = self.SymmLCDM_ids.astype(np.int32)
		self.FofrLCDM_ids = self.FofrLCDM_ids.astype(np.int32)

	def Get_legends(self, models):
		def append_legends(name, mod=0):
			if mod == 1:
				self.Symm_legends.append(name)
			elif mod == 2:
				self.fofr_legends.append(name)
			else:
				self.Symm_legends.append(name)
				self.fofr_legends.append(name)
			self.All_legends.append(name)
		
		for modnames in models:
			if modnames == 'lcdm':
				append_legends('$\Lambda$CDM')
				self.SymmLCDM_ids = np.append(self.SymmLCDM_ids, self.ID_counter)
				self.FofrLCDM_ids = np.append(self.FofrLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'symmA':
				append_legends('Symm_A', mod=1)
				self.SymmLCDM_ids = np.append(self.SymmLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'symmB':
				append_legends('Symm_B', mod=1)
				self.SymmLCDM_ids = np.append(self.SymmLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'symmC':
				append_legends('Symm_C', mod=1)
				self.SymmLCDM_ids = np.append(self.SymmLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'symmD':
				append_legends('Symm_D', mod=1)
				self.SymmLCDM_ids = np.append(self.SymmLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'fofr4':
				append_legends('fofr4', mod=2)
				self.FofrLCDM_ids = np.append(self.FofrLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'fofr5':
				append_legends('fofr5', mod=2)
				self.FofrLCDM_ids = np.append(self.FofrLCDM_ids, self.ID_counter)
				self.ID_counter += 1
			elif modnames == 'fofr6':
				append_legends('fofr6', mod=2)
				self.FofrLCDM_ids = np.append(self.FofrLCDM_ids, self.ID_counter)
				self.ID_counter += 1

	def savefigure(self, figure, name, savedir=False):
		""" Function that calls savefig based on figure instance and filename. """
		if not savedir:
			savedir_ = self.results_dir
		else:
			savedir_ = savedir
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(savedir_ + name + self.filetype, bbox_inches='tight', dpi=100)

	def Do_plot(self, EmptyHaloMass, NonEmptyHaloMass, AllHaloMass, NumCPs_list):
		NModels = len(AllHaloMass)
		SymmLCDM = self.SymmLCDM_ids
		FofrLCDM = self.FofrLCDM_ids
		Symm_only = self.SymmLCDM_ids[1:]
		Fofr_only = self.FofrLCDM_ids[1:]
		#####
		##### Computing some data
		#####
		### Number of empty vs non empty halos as a function of mass
		EmptyHaloMass = np.asarray(EmptyHaloMass)
		NonEmptyHaloMass = np.asarray(NonEmptyHaloMass)
		AllHaloMass = np.asarray(AllHaloMass)
		Common_bin_halo_mass = OF.Get_common_bin_logX(AllHaloMass, binnum=30)
		Number_EmptyHalo = []
		Number_NonEmpty = []
		for i in range(NModels):
			bin_value, bin_std = OF.Bin_numbers_common(EmptyHaloMass[i], EmptyHaloMass[i], Common_bin_halo_mass, std='poisson')
			Number_EmptyHalo.append(bin_value)
			bin_value, bin_std = OF.Bin_numbers_common(NonEmptyHaloMass[i], NonEmptyHaloMass[i], Common_bin_halo_mass, std='poisson')
			Number_NonEmpty.append(bin_value)
		### Number of critical points per halo
		Numer_CP_Halo_mass = []
		for i in range(NModels):
			index_bin = np.digitize(NonEmptyHaloMass[i], Common_bin_halo_mass)
			binval = np.array([np.sum(np.unique(NumCPs_list[i][j == index_bin])) for j in range(len(Common_bin_halo_mass))])
			print binval
		#####
		##### Plotting data
		#####
		#### Number of empty vs nonempty halos, as a function of mass
		Halo_mass_distribution_symm = plt.figure()
		ax = plt.axes()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(Common_bin_halo_mass, Number_EmptyHalo[i], linestyle='--', color=self.Plot_colors_symm[ij], alpha=0.5)
			ax.plot(Common_bin_halo_mass, Number_NonEmpty[i], linestyle='-', color=self.Plot_colors_symm[ij], label=self.Symm_legends[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		h,l=ax.get_legend_handles_labels()
		ax.legend(h,l)

		Halo_mass_distribution_fofr = plt.figure()
		ax = plt.axes()
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Common_bin_halo_mass, Number_EmptyHalo[i], linestyle='--', color=self.Plot_colors_fofr[ij], alpha=0.5)
			ax.plot(Common_bin_halo_mass, Number_NonEmpty[i], linestyle='-', color=self.Plot_colors_fofr[ij], label=self.fofr_legends[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		h,l=ax.get_legend_handles_labels()
		ax.legend(h,l)

		#### Number of critical points per halo
		CPs_per_halo_symm = plt.figure()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], NumCPs_list[i], 'o', color=self.Plot_colors_symm[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.Symm_legends)

		CPs_per_halo_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = SymmLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], NumCPs_list[i], 'o', color=self.Plot_colors_fofr[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.fofr_legends)

		print '--- SAVING IN: ', self.results_dir, ' ---'
		self.savefigure(Halo_mass_distribution_symm, 'Halo_mass_distribution_cSymmetron')
		self.savefigure(Halo_mass_distribution_fofr, 'Halo_mass_distribution_cFofr')
		self.savefigure(CPs_per_halo_symm, 'CP_per_halo_cSymmetron')
		self.savefigure(CPs_per_halo_fofr, 'CP_per_halo_cFofr')

def Argument_parser():
	""" Parses optional argument when program is run from the command line """
	#print 'Run python code with -h argument to see extra optional arguments'
	parser = argparse.ArgumentParser()
	# Optional arguments
	parser.add_argument("-Model", "--Model", help="Determines which model to run."\
	 				+ "Models may be: lcdm, symmX (X=A,B,C,D) or fofrY (Y=4,5,6). Use argument 'all' to run all models. Runs none by default."\
	 				+ "Do not run seperate models and 'all' at once! Defaults at lcdm", type=str, default='lcdm')
	parser.add_argument("-Nparts", "--NumberParticles", help="Run with a set number of particles. Default 64.", type=int, default=64)
	parser.add_argument("-Nsigma", "--Nsigma_disperse", help="Sigma value used for DisPerSE. Default at 3.", type=int, default=3)
	parser.add_argument("-filetype", "--Filetype", help="Filetype the figures are saved in. Can be pdf or png (default)", type=str, default='png')
	# Parse arguments
	args = parser.parse_args()
	# Some checks
	Model_ok = False
	Model_disperse = False
	Model_check = ['lcdm', 'symmA', 'symmB', 'symmC', 'symmD', 'fofr4', 'fofr5', 'fofr6', 'all']
	Models_to_be_run = []
	if args.Model:
		for models in Model_check:
			if args.Model == models:
				Model_ok = True
				if args.Model == 'all':
					Models_to_be_run = [Model_check[i] for i in range(len(Model_check)-1)]
				else:
					Models_to_be_run.append(args.Model)
		if not Model_ok:
			raise ValueError('Model input name %s not correctly set -Model argument.' %args.Model)
	return args, Models_to_be_run

if __name__ == '__main__':
	savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
	# Parse the arguments
	parsed_arguments, Models_to_be_run = Argument_parser()
	p_model = parsed_arguments.Model
	N_parts = parsed_arguments.NumberParticles
	N_sigma = parsed_arguments.Nsigma_disperse
	Filetype = parsed_arguments.Filetype
	print '====== Information ======'
	print 'Running for model:', p_model
	print 'Number of particles: ', N_parts
	print 'Sigma value: ', N_sigma
	print '========================='
	Empty_masses = []
	NonEmpty_masses = []
	All_masses = []
	Models_included = []
	NumCPs = []
	def append_data(mass, emass, nemass, NumCP, modelname):
		Empty_masses.append(emass)
		NonEmpty_masses.append(nemass)
		All_masses.append(mass)
		Models_included.append(modelname)
		NumCPs.append(NumCP)

	for p_model in Models_to_be_run:
		Instance = AnalyseCritPts(p_model, N_parts, N_sigma)
		NumCP_in_halo, CP_ids_in_halo, CP3_halo, CP2_halo = Instance.Check_cps()
		HaloMasses = Instance.HaloMass
		yesfils = NumCP_in_halo > 0
		nofils = NumCP_in_halo == 0
		print "Nonempty haloes: ", len(np.where(yesfils)[0])
		print "Total haloes: ", len(NumCP_in_halo)
		append_data(HaloMasses, HaloMasses[nofils], HaloMasses[yesfils], NumCP_in_halo[yesfils], p_model)

	Plot_instance = Plot_results(Models_included, N_sigma, 'ModelComparisons/HaloAnalysis/', filetype=Filetype)
	Plot_instance.Do_plot(Empty_masses, NonEmpty_masses, All_masses, NumCPs)
	"""
	ok_ids = np.array([])
	timer = time.time()
	nonempty_ids = CP_ids_in_halo[yesfils]
	for i in range(len(nonempty_ids)):
		ok_ids = np.unique(np.concatenate((ok_ids, nonempty_ids[i])))
	print len(np.unique(ok_ids))
	print len(Instance.CP_3DPos)
	print time.time() - timer, 's'
	"""