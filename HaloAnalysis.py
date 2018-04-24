# Common modules
import numpy as np
import cPickle as pickle
import argparse
import os
import time
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# Own modules
import ReadHaloData as RHD
import OtherFunctions as OF

# Global parameters
halfbox = 256./2.0  # Half the simulation box

class AnalyseCritPts():
	def __init__(self, model, npart, sigma):
		self.npart = npart
		self.sigma = sigma
		self.model = model
		if model == 'lcdm':
			Dispersemodel = 'LCDM'
		elif model == 'symmA':
			self.model = 'symm_A'
			Dispersemodel = 'SymmA'
		elif model == 'symmB':
			self.model = 'symm_B'
			Dispersemodel = 'SymmB'
		elif model == 'symmC':
			self.model = 'symm_C'
			Dispersemodel = 'SymmC'
		elif model == 'symmD':
			self.model = 'symm_D'
			Dispersemodel = 'SymmD'
		else:
			Dispersemodel = model
		self.read_filament_data(Dispersemodel, npart, sigma)
		self.HaloID, self.HaloPos, self.HaloMass, self.VirRadius = RHD.Read_halo_data(self.model)

	def read_filament_data(self, model, npart, sigma):
		""" Reads filament data from DisPerSE """
		print "===== Reading filament data of model " + model + " ... " 
		cachedir_foldername_extra = model + 'npart'+str(npart) + 'nsig'+str(sigma)
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/FilamentData/' + cachedir_foldername_extra + '/'
		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'Conditions_check.p'
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'
		# Read critical points from pickle and creates 3D coordinates of critical point positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(Filament_coordinate_data, CP_coordinate_data)
		CP_3DPos = []
		for j in range(len(self.CritPointXpos)):
			CP_3DPos.append(np.array((self.CritPointXpos[j], self.CritPointYpos[j], self.CritPointZpos[j])))
		self.CP_3DPos = np.array(CP_3DPos)
		Relevant_cps = self.CP_type >= 2
		Maxima_cp = self.CP_type[Relevant_cps] == 3
		Saddle_cp = self.CP_type[Relevant_cps] == 2
		self.CP_3DPos_max = self.CP_3DPos[Maxima_cp]
		self.CP_3DPos_saddle = self.CP_3DPos[Saddle_cp]
		self.Number_filaments_connecting_to_CP = self.Number_filaments_connecting_to_CP[self.CP_type >= 2]

	def Unpack_filament_data(self, Fil_coord, CP_coord):
		""" Unpacks filament and critial point data from the read filament data module. """
		self.NFils = Fil_coord[0]
		self.FilamentPos = Fil_coord[1]
		self.xdimPos = Fil_coord[2]
		self.ydimPos = Fil_coord[3]
		self.zdimPos = Fil_coord[4]
		self.NFilamentPoints = Fil_coord[5]
		self.FilID = Fil_coord[6]
		self.PairIDS = np.array(Fil_coord[7])        

		self.CritPointXpos = CP_coord[0]
		self.CritPointYpos = CP_coord[1]
		self.CritPointZpos = CP_coord[2]
		self.CP_type = np.array(CP_coord[3])
		self.CP_persistent_pair = CP_coord[4] 
		self.Critpts_filamentID = np.array(CP_coord[5]) 
		self.CP_id_of_connecting_filament = CP_coord[6]
		self.Number_filaments_connecting_to_CP = np.array(CP_coord[7])

	def CP_in_halo(self, HaloCord, CPCord, radius, norm_axis=1):
		""" Checks whether critical points are within a halo with a given radius """
		Difference = HaloCord - CPCord
		diffx = Difference[:,0]
		diffy = Difference[:,1]
		diffz = Difference[:,2]
		Difference[diffx <= -halfbox,0] += 256.0
		Difference[diffx >= halfbox,0] -= 256.0
		Difference[diffy <= -halfbox,1] += 256.0
		Difference[diffy >= halfbox,1] -= 256.0
		Difference[diffz <= -halfbox,2] += 256.0
		Difference[diffz >= halfbox,2] -= 256.0
		Diff_len = np.linalg.norm(Difference, axis=norm_axis)**2
		InHaloes = Diff_len <= radius**2
		return InHaloes

	def Check_cps(self):
		""" Checking how many critical points there are within each halo """
		cachedir_NumCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/NumCPInHalo/'
		cachedir_CPIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/CPIDsInHalo/'
		cachedir_MaximaCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/MaximaCPInHalo/'
		cachedir_SaddleCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/SaddleCPInHalo/'
		cachedir_MaximaCP_CPIDS = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/MaximaCPInHalo/CPIDsInHalo/'
		cachedir_SaddleCP_CPIDS = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/SaddleCPInHalo/CPIDsInHalo/'
		cachedir_NumFil = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/NumFilInHalo/'
		if not os.path.isdir(cachedir_NumCP):
			os.makedirs(cachedir_NumCP)
		if not os.path.isdir(cachedir_CPIDs):
			os.makedirs(cachedir_CPIDs)
		if not os.path.isdir(cachedir_MaximaCP):
			os.makedirs(cachedir_MaximaCP)
		if not os.path.isdir(cachedir_SaddleCP):
			os.makedirs(cachedir_SaddleCP)
		if not os.path.isdir(cachedir_MaximaCP_CPIDS):
			os.makedirs(cachedir_MaximaCP_CPIDS)
		if not os.path.isdir(cachedir_SaddleCP_CPIDS):
			os.makedirs(cachedir_SaddleCP_CPIDS)
		if not os.path.isdir(cachedir_NumFil):
			os.makedirs(cachedir_NumFil)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_NumCP = cachedir_NumCP + Common_filename
		cachefile_CPIDs = cachedir_CPIDs + Common_filename
		cachefile_MaximaCP = cachedir_MaximaCP + Common_filename
		cachefile_SaddleCP = cachedir_SaddleCP + Common_filename
		cachefile_MaximaCP_CPIDs = cachedir_MaximaCP_CPIDS + Common_filename
		cachefile_SaddleCP_CPIDs = cachedir_SaddleCP_CPIDS + Common_filename
		cachefile_NumFil = cachedir_NumFil + Common_filename

		if os.path.isfile(cachefile_NumCP):
			print "Reading CPs in halos of model " + self.model + "  ..."
			NumCP_in_halo = np.load(cachefile_NumCP)
			CP_ids_in_halo = np.load(cachefile_CPIDs)
			NumFils_in_halo = np.load(cachefile_NumFil)
		else:
			print "Calculating CPs in halos"
			timer = time.time()
			NumCP_in_halo = []
			CP_ids_in_halo = []
			NumFils_in_halo = []
			for i in range(len(self.HaloPos)):
				InHalo = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos, self.VirRadius[i])
				OK_cps = np.where(InHalo)[0]
				NumCP_in_halo.append(len(OK_cps))
				NumFils_in_halo.append(np.sum(self.Number_filaments_connecting_to_CP[OK_cps]))
				CP_ids_in_halo.append(OK_cps.astype(np.int32))
			NumCP_in_halo = np.asarray(NumCP_in_halo)
			CP_ids_in_halo = np.asarray(CP_ids_in_halo)
			NumFils_in_halo = np.asarray(NumFils_in_halo)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_NumCP, NumCP_in_halo)
			np.save(cachefile_CPIDs, CP_ids_in_halo)
			np.save(cachefile_NumFil, NumFils_in_halo)

		if os.path.isfile(cachefile_MaximaCP):
			print "Reading maxima and saddle cps in halos ..."
			MaximaCP_in_halo = np.load(cachefile_MaximaCP)
			SaddleCP_in_halo = np.load(cachefile_SaddleCP)
			CP_ids_in_halo_m = np.load(cachefile_MaximaCP_CPIDs)
			CP_ids_in_halo_s = np.load(cachefile_SaddleCP_CPIDs)
		else:
			print "Calculating CPs (maxima and saddle seperate) in halos"
			timer = time.time()
			MaximaCP_in_halo = []
			SaddleCP_in_halo = []
			CP_ids_in_halo_m = []
			CP_ids_in_halo_s = []
			for i in range(len(self.HaloPos)):
				InHalo_max = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos_max, self.VirRadius[i])
				InHalo_saddle = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos_saddle, self.VirRadius[i])
				OK_cps_m = np.where(InHalo_max)[0]
				OK_cps_s = np.where(InHalo_saddle)[0]
				MaximaCP_in_halo.append(len(OK_cps_m))
				SaddleCP_in_halo.append(len(OK_cps_s))
				CP_ids_in_halo_m.append(OK_cps_m.astype(np.int32))
				CP_ids_in_halo_s.append(OK_cps_s.astype(np.int32))
			MaximaCP_in_halo = np.asarray(MaximaCP_in_halo)
			SaddleCP_in_halo = np.asarray(SaddleCP_in_halo)
			CP_ids_in_halo_m = np.asarray(CP_ids_in_halo_m)
			CP_ids_in_halo_s = np.asarray(CP_ids_in_halo_s)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_MaximaCP, MaximaCP_in_halo)
			np.save(cachefile_SaddleCP, SaddleCP_in_halo)
			np.save(cachefile_MaximaCP_CPIDs, CP_ids_in_halo_m)
			np.save(cachefile_SaddleCP_CPIDs, CP_ids_in_halo_s)
		return NumCP_in_halo, CP_ids_in_halo, MaximaCP_in_halo, SaddleCP_in_halo, NumFils_in_halo, CP_ids_in_halo_m, CP_ids_in_halo_s

	def Check_cps_filtered(self):
		""" 
		Checking how many critical points there are within each halo 
		This function filters out filaments/critical points based on the data from ParticleAnalysis.py
		"""
		cachedir_NumCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/NumCPInHalo/'
		cachedir_CPIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/CPIDsInHalo/'
		cachedir_MaximaCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/MaximaCPInHalo/'
		cachedir_SaddleCP = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/SaddleCPInHalo/'
		cachedir_NumFil = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/NumFilInHalo/'
		cachedir_OKfils = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/IncludedFilaments/'
		if not os.path.isdir(cachedir_NumCP):
			os.makedirs(cachedir_NumCP)
		if not os.path.isdir(cachedir_CPIDs):
			os.makedirs(cachedir_CPIDs)
		if not os.path.isdir(cachedir_MaximaCP):
			os.makedirs(cachedir_MaximaCP)
		if not os.path.isdir(cachedir_SaddleCP):
			os.makedirs(cachedir_SaddleCP)
		if not os.path.isdir(cachedir_NumFil):
			os.makedirs(cachedir_NumFil)
			
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_NumCP = cachedir_NumCP + Common_filename
		cachefile_CPIDs = cachedir_CPIDs + Common_filename
		cachefile_MaximaCP = cachedir_MaximaCP + Common_filename
		cachefile_SaddleCP = cachedir_SaddleCP + Common_filename
		cachefile_NumFil = cachedir_NumFil + Common_filename
		cachefile_okfils = cachedir_OKfils + Common_filename

		# Obtains filament length and filters out filaments smaller than 1 Mpc/h
		Included_fils = np.load(cachefile_okfils)
		Filament_3DPos = []
		for j in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[j], self.ydimPos[j], self.zdimPos[j])))
		Filament_3DPos = np.array(Filament_3DPos)
		FilamentLength = []
		for i in range(len(Filament_3DPos)):
			fillen = OF.Get_filament_length(Filament_3DPos[i])
			FilamentLength.append(fillen)
		FilamentLength = np.asarray(FilamentLength)
		Small_filaments = FilamentLength > 1.0   # Filter for filaments smaller than 1 Mpc/h
		self.NFils_filter = len(FilamentLength[Small_filaments][Included_fils])
		# Filter pair ids (CP end points of filaments) based on filament length
		PairIDS_lenfilt = self.PairIDS[Small_filaments]
		# Filter pair ids based on density computed from ParticleAnalysis.py
		PairIDS_denfilt = PairIDS_lenfilt[Included_fils]
		Accepted_CP_ids = np.concatenate((np.unique(PairIDS_denfilt[:,0]),np.unique(PairIDS_denfilt[:,1])))
		self.Critpts_filamentID_filter = self.Critpts_filamentID[Accepted_CP_ids]
		Accepted_CP_ids_o2 = np.unique(PairIDS_denfilt[:,0])
		Accepted_CP_ids_o3 = np.unique(PairIDS_denfilt[:,1])
		self.CP_3DPos_filter = self.CP_3DPos[Accepted_CP_ids]
		CP_3DPos_filter_02 = self.CP_3DPos[Accepted_CP_ids_o2]
		CP_3DPos_filter_03 = self.CP_3DPos[Accepted_CP_ids_o3]
		self.Number_filaments_connecting_to_CP_filter = self.Number_filaments_connecting_to_CP[Accepted_CP_ids]

		if os.path.isfile(cachefile_NumCP):
			print "Reading CPs in halos of model " + self.model + "  ..."
			NumCP_in_halo = np.load(cachefile_NumCP)
			CP_ids_in_halo = np.load(cachefile_CPIDs)
			NumFils_in_halo = np.load(cachefile_NumFil)
		else:
			print "Calculating CPs in halos, using filtered filaments"
			timer = time.time()
			NumCP_in_halo = []
			CP_ids_in_halo = []
			NumFils_in_halo = []
			for i in range(len(self.HaloPos)):
				InHalo = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos_filter, self.VirRadius[i])
				OK_cps = np.where(InHalo)[0]
				NumCP_in_halo.append(len(OK_cps))
				NumFils_in_halo.append(np.sum(self.Number_filaments_connecting_to_CP_filter[OK_cps]))
				CP_ids_in_halo.append(OK_cps.astype(np.int32))
			NumCP_in_halo = np.asarray(NumCP_in_halo)
			CP_ids_in_halo = np.asarray(CP_ids_in_halo)
			NumFils_in_halo = np.asarray(NumFils_in_halo)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_NumCP, NumCP_in_halo)
			np.save(cachefile_CPIDs, CP_ids_in_halo)
			np.save(cachefile_NumFil, NumFils_in_halo)

		if os.path.isfile(cachefile_MaximaCP):
			print "Reading maxima and saddle cps in halos ..."
			MaximaCP_in_halo = np.load(cachefile_MaximaCP)
			SaddleCP_in_halo = np.load(cachefile_SaddleCP)
		else:
			print "Calculating CPs (maxima and saddle seperate) in halos"
			timer = time.time()
			MaximaCP_in_halo = []
			SaddleCP_in_halo = []
			for i in range(len(self.HaloPos)):
				InHalo_max = self.CP_in_halo(self.HaloPos[i], CP_3DPos_filter_02, self.VirRadius[i])
				InHalo_saddle = self.CP_in_halo(self.HaloPos[i], CP_3DPos_filter_02, self.VirRadius[i])
				OK_cps_m = np.where(InHalo_max)[0]
				OK_cps_s = np.where(InHalo_saddle)[0]
				MaximaCP_in_halo.append(len(OK_cps_m))
				SaddleCP_in_halo.append(len(OK_cps_s))
			MaximaCP_in_halo = np.asarray(MaximaCP_in_halo)
			SaddleCP_in_halo = np.asarray(SaddleCP_in_halo)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_MaximaCP, MaximaCP_in_halo)
			np.save(cachefile_SaddleCP, SaddleCP_in_halo)
		return NumCP_in_halo, CP_ids_in_halo, MaximaCP_in_halo, SaddleCP_in_halo, NumFils_in_halo

	def Check_halo_in_cp(self):
		""" Checking how many critical points there are within each halo """
		cachedir_NumHalo = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/HaloInCP/Number_halos_per_cp/'
		cachedir_HaloIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/HaloInCP/Halo_ids_per_cp/'
		if not os.path.isdir(cachedir_NumHalo):
			os.makedirs(cachedir_NumHalo)
		if not os.path.isdir(cachedir_HaloIDs):
			os.makedirs(cachedir_HaloIDs)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_NumHalo = cachedir_NumHalo + Common_filename
		cachefile_HaloIDs = cachedir_HaloIDs + Common_filename
		if os.path.isfile(cachefile_NumHalo):
			NumHalo_in_cp = np.load(cachefile_NumHalo)
			Halo_ids_in_cp = np.load(cachefile_HaloIDs)
		else:
			print "Computing number of halos per critical point ..."
			timer = time.time()
			NumHalo_in_cp = []
			Halo_ids_in_cp = []
			for i in range(len(self.CP_3DPos)):
				InHalo = self.CP_in_halo(self.HaloPos, self.CP_3DPos[i], self.VirRadius)
				OK_cps = np.where(InHalo)[0]
				NumHalo_in_cp.append(len(OK_cps))
				Halo_ids_in_cp.append(OK_cps.astype(np.int32))
			NumHalo_in_cp = np.asarray(NumHalo_in_cp)
			Halo_ids_in_cp = np.asarray(Halo_ids_in_cp)
			print "time took: ", time.time() - timer, "s"
			np.save(cachefile_NumHalo, NumHalo_in_cp)
			np.save(cachefile_HaloIDs, Halo_ids_in_cp)
		return NumHalo_in_cp, Halo_ids_in_cp

class Plot_results():
	def __init__(self, models, Nsigma, foldername, filetype='png', filter_fils=False):
		if filetype == 'png':
			foldername += 'PNG/'
		elif filetype == 'pdf':
			foldername += 'PDF/'
		self.filter_fils = filter_fils
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
		if self.filter_fils:
			figure.savefig(savedir_ + name + 'FilteredFilament' + self.filetype, bbox_inches='tight', dpi=100)
		else:
			figure.savefig(savedir_ + name + self.filetype, bbox_inches='tight', dpi=100)

	def Do_plot(self, EmptyHaloMass, NonEmptyHaloMass, AllHaloMass, NumCPs_list, Number_filaments):
		NModels = len(NumCPs_list)
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
			#print binval
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
		plt.ylabel('$N$ halos')
		plt.xscale('log')
		plt.yscale('log')
		h,l=ax.get_legend_handles_labels()
		ax.legend(h,l)

		Halo_mass_distribution_fofr = plt.figure()
		ax = plt.axes()
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Common_bin_halo_mass, Number_EmptyHalo[i], linestyle='--', color=self.Plot_colors_fofr[ij], alpha=0.5)
			ax.plot(Common_bin_halo_mass, Number_NonEmpty[i], linestyle='-', color=self.Plot_colors_fofr[ij], label=self.fofr_legends[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$ halos')
		plt.xscale('log')
		plt.yscale('log')
		h,l=ax.get_legend_handles_labels()
		ax.legend(h,l)

		#### Number of critical points per halo
		CPs_per_halo_symm = plt.figure()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], NumCPs_list[i], 'o', color=self.Plot_colors_symm[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$ critical points')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.Symm_legends)

		CPs_per_halo_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], NumCPs_list[i], 'o', color=self.Plot_colors_fofr[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$ critical points')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.fofr_legends)
		#### Number of filaments connected to a halo as a function of the halo mass
		Filaments_per_halo_symm = plt.figure()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], Number_filaments[i], 'o', color=self.Plot_colors_symm[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$ filaments')
		plt.xscale('log')
		plt.legend(self.Symm_legends)
        
		Filaments_per_halo_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(NonEmptyHaloMass[i], Number_filaments[i], 'o', color=self.Plot_colors_fofr[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$ filaments')
		plt.xscale('log')
		plt.legend(self.fofr_legends)
                
		print '--- SAVING IN: ', self.results_dir, ' ---'
		self.savefigure(Halo_mass_distribution_symm, 'Halo_mass_distribution_cSymmetron')
		self.savefigure(Halo_mass_distribution_fofr, 'Halo_mass_distribution_cFofr')
		self.savefigure(CPs_per_halo_symm, 'CP_per_halo_cSymmetron')
		self.savefigure(CPs_per_halo_fofr, 'CP_per_halo_cFofr')
		self.savefigure(Filaments_per_halo_symm, 'Filaments_per_halo_cSymmetron')
		self.savefigure(Filaments_per_halo_fofr, 'Filaments_per_halo_cFofr')
        
	def Order2_3_plots(self, AllHaloMass, O2_Cps, O3_Cps, O2_masses, O3_masses):
		NModels = len(O2_Cps)
		SymmLCDM = self.SymmLCDM_ids
		FofrLCDM = self.FofrLCDM_ids
		Symm_only = self.SymmLCDM_ids[1:]
		Fofr_only = self.FofrLCDM_ids[1:]
		#####
		##### Computing some data
		#####
		### Number of empty vs non empty halos as a function of mass
		AllHaloMass = np.asarray(AllHaloMass)
		Common_bin_halo_mass = OF.Get_common_bin_logX(AllHaloMass, binnum=30)
		#####
		##### Plotting data
		#####
		### Plotting order 2 critical points as a function of mass 
		O2CPs_per_halo_symm = plt.figure()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(O2_masses[i], O2_Cps[i], 'o', color=self.Plot_colors_symm[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.Symm_legends)

		O2CPs_per_halo_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = SymmLCDM[ij]
			plt.plot(O2_masses[i], O2_Cps[i], 'o', color=self.Plot_colors_fofr[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.fofr_legends)

		### Plotting order 3 critical points as a function of mass 
		O3CPs_per_halo_symm = plt.figure()
		for ij in range(len(SymmLCDM)):
			i = SymmLCDM[ij]
			plt.plot(O3_masses[i], O3_Cps[i], 'o', color=self.Plot_colors_symm[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.Symm_legends)

		O3CPs_per_halo_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = SymmLCDM[ij]
			plt.plot(O3_masses[i], O3_Cps[i], 'o', color=self.Plot_colors_fofr[ij])
		plt.xlabel('Halo mass - $[M_\odot/h]$')
		plt.ylabel('$N$')
		plt.xscale('log')
		#plt.yscale('log')
		plt.legend(self.fofr_legends)
		print '--- SAVING IN: ', self.results_dir, ' ---'
		self.savefigure(O2CPs_per_halo_symm, 'O2CP_per_halo_cSymmetron')
		self.savefigure(O2CPs_per_halo_fofr, 'O2CP_per_halo_cFofr')
		self.savefigure(O3CPs_per_halo_symm, 'O3CP_per_halo_cSymmetron')
		self.savefigure(O3CPs_per_halo_fofr, 'O3CP_per_halo_cFofr')

	def Other_plots(self, OK_PC, Out_PC, OK_PC_alt, Out_PC_alt, EmptyH, NonEmptyH, Nfilaments, Nfilaments_filter, NumUniqueFils):
		def autolabel(rects, do_int=False):
			"""
			Attach a text label above each bar displaying its height
			"""
			for rect in rects:
				height = rect.get_height()
				if do_int:
					ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,'%.d' % int(height), ha='center', va='bottom')
				else:
					ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,'%.1f' % float(height), ha='center', va='bottom')
		NModels = len(OK_PC)
		width=0.35
		Percentage_CP_in_halo, ax = plt.subplots()
		ind = np.arange(NModels)
		rects1 = ax.bar(ind, OK_PC, width, color='b')
		rects2 = ax.bar(ind+width, Out_PC, width, color='r')
		ax.set_ylabel('% of critical points')
		ax.set_xticks(ind + width/2)
		ax.set_xticklabels((r'$\Lambda$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6'))
		ax.legend((rects1[0], rects2[0]), ('In halo', 'Outside halo'), loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		autolabel(rects1)
		autolabel(rects2)

		#NumberHalos_per_CP, ax = plt.subplots()
		#ind = np.arange(NModels)
		#rects1 = ax.bar(ind, OK_PC_alt, width, color='b')
		#rects2 = ax.bar(ind+width, Out_PC_alt, width, color='r')
		#ax.set_ylabel('% of critical points')
		#ax.set_xticks(ind + width/2)
		#ax.set_xticklabels((r'$\Lambda$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6'))
		#ax.legend((rects1[0], rects2[0]), ('In halo', 'Outside halo'), loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		#autolabel(rects1)
		#autolabel(rects2)

		HalosWithCPs, ax = plt.subplots()
		ind = np.arange(NModels)
		rects1 = ax.bar(ind, NonEmptyH, width, color='b')
		rects2 = ax.bar(ind+width, EmptyH, width, color='r')
		ax.set_ylabel('% of critical points')
		ax.set_xticks(ind + width/2)
		ax.set_xticklabels((r'$\Lambda$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6'))
		ax.legend((rects1[0], rects2[0]), ('Halos with critpt(s)', 'Halos without critpt'), loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		autolabel(rects1, do_int=True)
		autolabel(rects2, do_int=True)

		NumberFilaments, ax = plt.subplots()
		ind = np.arange(NModels)
		rects1 = ax.bar(ind, Nfilaments, width, color='b')
		rects2 = ax.bar(ind+width, Nfilaments_filter, width, color='r')
		ax.set_ylabel('$N$ filaments')
		ax.set_xticks(ind + width/2)
		#ax.set_xticklabels((r'$\Lambda$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6'))
		ax.set_xticklabels(self.All_legends)
		ax.legend((rects1[0], rects2[0]), ('Total filaments', 'Filaments after \n filtering CPs'), loc = 'lower left',
				 bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		autolabel(rects1, do_int=True)
		autolabel(rects2, do_int=True)

		NumberFilaments_unique, ax = plt.subplots()
		ind = np.arange(NModels)
		rects1 = ax.bar(ind, Nfilaments, width, color='b')
		rects2 = ax.bar(ind+width, NumUniqueFils, width, color='r')
		ax.set_ylabel('$N$ filaments')
		ax.set_xticks(ind + width/2)
		#ax.set_xticklabels((r'$\Lambda$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6'))
		ax.set_xticklabels(self.All_legends)
		ax.legend((rects1[0], rects2[0]), ('Total filaments', 'Filaments after \n filtering CPs'), loc = 'lower left',
				 bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		autolabel(rects1, do_int=True)
		autolabel(rects2, do_int=True)

		self.savefigure(Percentage_CP_in_halo, 'Percentage_CP_in_halo')
		#self.savefigure(NumberHalos_per_CP, 'NumberHalos_per_CP')
		self.savefigure(HalosWithCPs, 'HalosWithCps')
		self.savefigure(NumberFilaments, 'Number_filaments_after_cp_filter')
		self.savefigure(NumberFilaments_unique, 'Number_filaments_after_cp_filter_UNIQUE')
    
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
	# Global data
	All_masses = []
	Models_included = []
	# Non filtered data
	Empty_masses = []
	NonEmpty_masses = []
	NumCPs = []
	Order2_cps = []
	Order3_cps = []
	Order2_mass = []
	Order3_mass = []
	NumFilaments = []
	InHaloPC = []
	OutHaloPC = []
	InHalo_alternative = []
	OutHalo_alternative = []
	Empty_halos = []
	NonEmpty_halos = []
	Numfils_all = []
	Numfils_filter = []
	Num_unique_fils = []
	# Filtered data
	Models_included_filter = []
	Empty_masses_f = []
	NonEmpty_masses_f = []
	NumCPs_f = []
	Order2_cps_f = []
	Order3_cps_f = []
	Order2_mass_f = []
	Order3_mass_f = []
	NumFilaments_f = []
	InHaloPC_f = []
	OutHaloPC_f = []
	Empty_halos_f = []
	NonEmpty_halos_f = []
	Numfils_all_f = []
	Numfils_filter_f = []
	Num_unique_fils_f = []

	def append_data(mass, emass, nemass, NumCP, O2CP, O3CP, O2Mass, O3Mass, Nfils, modelname):
		Empty_masses.append(emass)
		NonEmpty_masses.append(nemass)
		All_masses.append(mass)
		Models_included.append(modelname)
		NumCPs.append(NumCP)
		Order2_cps.append(O2CP)
		Order3_cps.append(O3CP)
		Order2_mass.append(O2Mass)
		Order3_mass.append(O3Mass)
		NumFilaments.append(Nfils)

	def append_data_filter(emass, nemass, NumCP, O2CP, O3CP, O2Mass, O3Mass, Nfil):
		Empty_masses_f.append(emass)
		NonEmpty_masses_f.append(nemass)
		NumCPs_f.append(NumCP)
		Order2_cps_f.append(O2CP)
		Order3_cps_f.append(O3CP)
		Order2_mass_f.append(O2Mass)
		Order3_mass_f.append(O3Mass)
		NumFilaments_f.append(Nfil)
        
	def Percentage_computing(nonempty_ids_, CPPos, model, npart, sigma, do_filter=False, savefile=True):
		if do_filter:
			cachedir_PCT = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/CPInHaloPercentage/'
			cachedir_okID = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/FilteredFils/CPInHaloPercentage/ok_ids/'
		else:
			cachedir_PCT = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/CPInHaloPercentage/'
			cachedir_okID = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/HaloAnalysis/CPInHaloPercentage/ok_ids/'
		if not os.path.isdir(cachedir_PCT):
			os.makedirs(cachedir_PCT)
		if not os.path.isdir(cachedir_okID):
			os.makedirs(cachedir_okID)
		Common_filename =  model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_PCT = cachedir_PCT + Common_filename
		cachefile_okID = cachedir_okID + Common_filename
		if os.path.isfile(cachefile_PCT) and savefile:
			savedPCT = np.load(cachefile_PCT)
			savedokIDs = np.load(cachefile_okID)
			Percentage = savedPCT[0]
			Out_percentage = savedPCT[1]
			ok_ids_pc = savedokIDs
		else:
			ok_ids_pc = np.array([])
			timer = time.time()
			for i in range(len(nonempty_ids_)):
				ok_ids_pc = np.unique(np.concatenate((ok_ids_pc, nonempty_ids_[i])))
			Percentage = 100.0*float(len(ok_ids_pc))/len(CPPos)
			ok_ids_pc = ok_ids_pc.astype(np.int32)
			Out_percentage = 100.0 - Percentage
			np.save(cachefile_PCT, np.array([Percentage, Out_percentage]))
			np.save(cachefile_okID, ok_ids_pc)
			print "percentage time: ", time.time() - timer, 's'
		return Percentage, Out_percentage, ok_ids_pc

	for p_model in Models_to_be_run:
		Instance = AnalyseCritPts(p_model, N_parts, N_sigma)
		# For non filtered filaments
		NumCP_in_halo, CP_ids_in_halo, CP3_halo, CP2_halo, NumberFils, CP_ids_in_halo_O3, CP_ids_in_halo_O2 = Instance.Check_cps()
		HaloMasses = Instance.HaloMass
		yesfils = NumCP_in_halo > 0
		nofils = NumCP_in_halo == 0
		Empty_halos.append(len(NumCP_in_halo[nofils]))
		NonEmpty_halos.append(len(NumCP_in_halo[yesfils]))
		yesfils_o3 = CP3_halo > 0
		yesfils_o2 = CP2_halo > 0
		append_data(HaloMasses, HaloMasses[nofils], HaloMasses[yesfils], NumCP_in_halo[yesfils], CP2_halo[yesfils_o2],
					CP3_halo[yesfils_o3], HaloMasses[yesfils_o2], HaloMasses[yesfils_o3], NumberFils[yesfils], p_model)
		PCT, OutPCT, ok_ids = Percentage_computing(CP_ids_in_halo[yesfils], Instance.CP_3DPos, p_model, N_parts, N_sigma)
		InHaloPC.append(PCT)
		OutHaloPC.append(OutPCT)
		Numfils_all.append(np.sum(Instance.NFils))
		Numfils_filter.append(np.sum(Instance.Number_filaments_connecting_to_CP[ok_ids]))
		Unique_fil_ids = np.array([])
		for i in range(len(Instance.Critpts_filamentID[ok_ids])):
			Unique_fil_ids = np.unique(np.concatenate((Unique_fil_ids, Instance.Critpts_filamentID[ok_ids][i])))
		Num_unique_fils.append(len(Unique_fil_ids))
		# For filtered filaments
		if (p_model != 'symmC' and N_parts == 64) or N_parts == 188:
			Models_included_filter.append(p_model)
			NumCP_in_halo_f, CP_ids_in_halo_f, CP3_halo_f, CP2_halo_f, NumberFils_f = Instance.Check_cps_filtered()
			yesfils_f = NumCP_in_halo_f > 0
			nofils_f = NumCP_in_halo_f == 0
			Empty_halos_f.append(len(NumCP_in_halo_f[nofils_f]))
			NonEmpty_halos_f.append(len(NumCP_in_halo_f[yesfils_f]))
			yesfils_o3_f = CP3_halo_f > 0
			yesfils_o2_f = CP2_halo_f > 0
			append_data_filter(HaloMasses[nofils_f], HaloMasses[yesfils_f], NumCP_in_halo_f[yesfils_f], CP2_halo_f[yesfils_o2_f],
						CP3_halo_f[yesfils_o3_f], HaloMasses[yesfils_o2_f], HaloMasses[yesfils_o3_f], NumberFils_f[yesfils_f])

			PCT_f, OutPCT_f, ok_ids_f = Percentage_computing(CP_ids_in_halo_f[yesfils_f], Instance.CP_3DPos_filter, p_model, N_parts, N_sigma, do_filter=True)
			InHaloPC_f.append(PCT_f)
			OutHaloPC_f.append(OutPCT_f)
			Numfils_all_f.append(np.sum(Instance.NFils_filter))
			Numfils_filter_f.append(np.sum(Instance.Number_filaments_connecting_to_CP_filter[ok_ids_f]))
			Unique_fil_ids_f = np.array([])
			###
			### CHECK WHY THERE ARE MORE FILAMENTS COMPARED TO FILTERED ONES
			### DO A DENSITY PLOTS VS FILAMENTS IN A SLICED BOX, COPY FROM DISPERSE_ANALYSIS
			### 
			for i in range(len(Instance.Critpts_filamentID_filter[ok_ids_f])):
				Unique_fil_ids_f = np.unique(np.concatenate((Unique_fil_ids_f, Instance.Critpts_filamentID_filter[ok_ids_f][i])))
			Num_unique_fils_f.append(len(Unique_fil_ids_f))
		elif p_model == 'symmC' and N_parts == 64:
			Dummy = np.array([0,0,0,0,0])
			append_data_filter(Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy)
			InHaloPC_f.append(0)
			OutHaloPC_f.append(0)
			Empty_halos_f.append(0)
			NonEmpty_halos_f.append(0)
			Numfils_all_f.append(0)
			Numfils_filter_f.append(0)
			Num_unique_fils_f.append(0)
		print "Nonempty haloes: ", len(np.where(yesfils)[0])
		#print "Nonempty haloes after filter: ", len(np.where(yesfils_f)[0])
		print "Total haloes: ", len(NumCP_in_halo)
		#print stop
		#CPinhalo, HaloIds_inhalo = Instance.Check_halo_in_cp()
		#Percentage_inhalo = 100.0*(len(CPinhalo[CPinhalo > 0]))/len(CPinhalo)
		#Percentage_outhalo = 100.0 - Percentage_inhalo
		#InHalo_alternative.append(Percentage_inhalo)
		#OutHalo_alternative.append(Percentage_outhalo)

	Plot_instance = Plot_results(Models_included, N_sigma, 'ModelComparisons/HaloAnalysis/', filetype=Filetype)
	Plot_instance.Do_plot(Empty_masses, NonEmpty_masses, All_masses, NumCPs, NumFilaments)
	Plot_instance.Order2_3_plots(All_masses, Order2_cps, Order3_cps, Order2_mass, Order3_mass)
	Plot_instance.Other_plots(InHaloPC, OutHaloPC, InHalo_alternative, OutHalo_alternative, Empty_halos, NonEmpty_halos, Numfils_all,
							 Numfils_filter, Num_unique_fils)
	print "Plotting for filtered filaments"
	Plot_instance2 = Plot_results(Models_included, N_sigma, 'ModelComparisons/HaloAnalysis/', filetype=Filetype, filter_fils=True)
	Plot_instance2.Do_plot(Empty_masses_f, NonEmpty_masses_f, All_masses, NumCPs_f, NumFilaments_f)
	Plot_instance2.Order2_3_plots(All_masses, Order2_cps_f, Order3_cps_f, Order2_mass_f, Order3_mass_f)
	Plot_instance2.Other_plots(InHaloPC_f, OutHaloPC_f, 1, 2, Empty_halos_f, NonEmpty_halos_f, Numfils_all_f, Numfils_filter_f, Num_unique_fils_f)