import numpy as np
import cPickle as pickle
import ReadGadgetFile as RGF
import OtherFunctions as OF

class ParticleAnalysis():
	def __init__(self, model, npart, sigma):
		if model == 'lcdm':
			self.Dispersemodel = 'LCDM'
		elif model == 'symm_A':
			self.Dispersemodel = 'SymmA'
		elif model == 'symm_B':
			self.Dispersemodel = 'SymmD'
		elif model == 'symm_C':
			self.Dispersemodel = 'SymmC'
		elif model == 'symm_D':
			self.Dispersemodel = 'SymmD'
		else:
			self.Dispersemodel = model

		self.Read_basic_data

	def Read_basic_data(self, model, npart, sigma):
		""" Reads filament data from DisPerSE """
		cachedir_foldername_extra = self.Dispersemodel + 'npart'+str(npart)
		if sigma:
			cachedir_foldername_extra += 'nsig'+str(sigma)
		else:
			cachedir_foldername_extra += 'nsig3'
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/FilamentData/' + cachedir_foldername_extra + '/'
		
		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'Conditions_check.p'
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'
		# Read filament from pickle and creates 3D coordinates of filament positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(Box_boundaries, Filament_coordinate_data)
		for j in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[j], self.ydimPos[j], self.zdimPos[j])))
		self.Filament_3DPos = np.array(Filament_3DPos)

		# Reads particle data from Gadget
		MaskDM = np.array([0,0,0,0,0,0])*256.0
		GadgetInstance = RGF.Read_Gadget_file([0,0,0],MaskDM)
		PartPosX, PartPosY, PartPosZ, PartID, HistDM, BinXEDM, BinYEDM, KDTREE = GadgetInstance.Get_particles(model, includeKDTree=False)
		Part_velx, Part_vely, Part_velz = GadgetInstance.Get_velocities()
		
		self.ParticlePos, self.ParticleVel = OF.Get_3D(PartPosX, PartPosY, PartPosZ, Part_velx, Part_vely, Part_velz)

	def Unpack_filament_data(self, Box_info, Fil_coord):
		""" Unpacks filament data from the read filament data module. """
		self.xmin = Box_info[0]
		self.xmax = Box_info[1]
		self.ymin = Box_info[2]
		self.ymax = Box_info[3]
		self.zmin = Box_info[4]
		self.zmax = Box_info[5]

		self.NFils = Fil_coord[0]
		self.FilamentPos = Fil_coord[1]
		self.xdimPos_nonBC = Fil_coord[2]
		self.ydimPos_nonBC = Fil_coord[3]
		self.zdimPos_nonBC = Fil_coord[4]
		self.NFilamentPoints = Fil_coord[5]
		self.FilID = Fil_coord[6]
		self.PairIDS = Fil_coord[7]


	def Read_particle_data(self, model, npart, sigma, boxexpand):
		""" Reads particle boxes, masked IDs and computed particle distances """
		# Distances computed
		cachedir_ppf_distances = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Distances/'
		cachefile_distances = cachedir_ppf_distances + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + 'Analytic_Periodic.p'
		self.Distances = pickle.load(open(cachefile_distances, 'rb'))
		# Segment the particle belongs to
		cachedir_ppf_segIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/SegIDs/'
		cachefile_segIDs = cachedir_ppf_segIDs  + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + 'Analytic_Periodic.p'
		self.Segment_index = pickle.load(open(cachefile_segIDs, 'rb'))
		# T value that gives shortest distance, D(t) = |s(t) - p|
		cachedir_ppf_tsols = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Tsolutions/'
		cachefile_tsols = cachedir_ppf_tsols + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + 'Analytic_Periodic.p'
		self.Tsolutions = pickle.load(open(cachefile_tsols, 'rb'))
		# Particle IDs of the particle box
		cachedir_ppf_masks = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/MaskedParticlesIDs/'
		cachefile_masks = cachedir_ppf_masks + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '_Periodic.p'
		self.Masked_IDs = pickle.load(open(cachefile_masks, 'rb'))
		# Particle position of the particle box
		cachedir_ppf_pbox = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ParticleBoxes/'
		cachefile_pbox = cachedir_ppf_pbox + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '_Periodic.p'
		self.ParticleBoxes = pickle.load(open(cachefile_pbox, 'rb'))

	def accepted_distance_func(self, distance_arr, binnum=150):
		"""
		Algorithm that determines the accepted particle distance.
		See thesis for details
		Input are computed distances to a filament
		"""
		maximum_thickness = 0.7   # Units of Mpc/h
		hist, bins_list = np.histogram(distance_arr, bins=binnum)
		# Find largest peak in histogram, ensures it smaller than a given threshold
		Max_bin_idx = np.where(hist == np.max(hist))[0]
		Max_bin_idx = Max_bin_idx[0]
		if bins_list[Max_bin_idx] > maximum_thickness:
			nlargest = -2
			for i in range(len(hist)):
				Max_bin_idx_new = np.where(hist == np.partition(hist, nlargest)[nlargest])[0]
				if bins_list[Max_bin_idx_new[0]] < maximum_thickness:
					Max_bin_idx = Max_bin_idx_new[0]
					break
				else:
					nlargest -= 1
		accept_distance_threshold = bins_list[Max_bin_idx]
		current_count = hist[Max_bin_idx]
		# Loops through increasing bin and see if histogram increases or decreases
		counter = 0
		for i in range(Max_bin_idx, len(hist)):
			count = hist[i]
			if count > current_count and counter == 0:
				accept_distance_threshold = bins_list[i]
				current_count = count
				counter += 1
			elif count > current_count and counter == 1:
				accept_distance_threshold = bins_list[i-1]
				break
			elif count < current_count:
				accept_distance_threshold = bins_list[i]
				current_count = count
				counter = 0
		return accept_distance_threshold

	def Get_accepted_particles(self, Distance_threshold):
		Particles_accepted = []
		Distances_accepted = []
		Particles_not_accepted = []
		Accepted_box_particles = []
		for i in range(len(Distances)):
			Filament_thickness = np.where(self.Distances[i] <= Distance_threshold[i])[0]
			Bigger_distances = np.where(self.Distances[i] > Distance_threshold[i])[0]
			Accepted_box_particles.append(Filament_thickness)
			Particles_accepted.append(self.Masked_IDs[i][Filament_thickness])
			Particles_not_accepted.append(self.Masked_IDs[i][Bigger_distances])
			Distances_accepted.append(self.Distances[i][Filament_thickness])

		self.Accepted_box_particles = np.asarray(Accepted_box_particles)
		self.Particles_accepted = np.asarray(Particles_accepted)
		self.Distances_accepted = np.asarray(Distances_accepted)
		Particles_not_accepted = np.asarray(Particles_not_accepted)
