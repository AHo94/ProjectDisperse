import numpy as np
import cPickle as pickle
import ReadGadgetFile as RGF
import OtherFunctions as OF
import os

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

		self.model = model
		self.npart = npart
		self.sigma = sigma
		self.Read_basic_data(model, npart, sigma)
		self.Read_particle_data(model, npart, sigma, 3)
		self.Do_filter_particles()   # Filters the particles

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
		self.FilamentLength = []
		for i in range(len(Filament_3DPos)):
			fillen = OF.Get_filament_length(Filament_3DPos[i], self.xmax-self.xmin)
			self.FilamentLength.append(fillen)
		self.FilamentLength = np.asarray(self.FilamentLength)

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
		#Common_filename =  model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		cachedir_ppf_distances = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Distances/'
		cachefile_distances = cachedir_ppf_distances + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		self.Distances = pickle.load(open(cachefile_distances, 'rb'))
		# Segment the particle belongs to
		cachedir_ppf_segIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/SegIDs/'
		cachefile_segIDs = cachedir_ppf_segIDs  + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		self.Segment_index = pickle.load(open(cachefile_segIDs, 'rb'))
		# T value that gives shortest distance, D(t) = |s(t) - p|
		cachedir_ppf_tsols = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Tsolutions/'
		cachefile_tsols = cachedir_ppf_tsols + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		self.Tsolutions = pickle.load(open(cachefile_tsols, 'rb'))
		# Particle IDs of the particle box
		cachedir_ppf_masks = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/MaskedParticlesIDs/'
		cachefile_masks = cachedir_ppf_masks + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		self.Masked_IDs = pickle.load(open(cachefile_masks, 'rb'))
		# Particle position of the particle box
		cachedir_ppf_pbox = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ParticleBoxes/'
		cachefile_pbox = cachedir_ppf_pbox + model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.p'
		self.ParticleBoxes = pickle.load(open(cachefile_pbox, 'rb'))

	def accepted_distance_func(self, distance_arr, binnum=150):
		"""
		Algorithm that determines the accepted particle distance.
		See thesis for details
		Input are computed distances to a filament
		"""
		maximum_thickness = 0.7   # Units of Mpc/h
		hist, bins_list = np.histogram(distance_arr, bins='fd')
		# Find largest peak in histogram, ensures it smaller than a given threshold
		Max_bin_idx = np.where(hist == np.max(hist))[0]
		Max_bin_idx = Max_bin_idx[0]
		if bins_list[Max_bin_idx] > maximum_thickness:
			nlargest = -2
			for i in range(len(hist)-1):
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

	def Filter_halo_particles(self, FilamentPos, SegIDs, Tsols, fillen, Accepted_box_part, binnum=40):
		""" 
		Function that filters halo particles. 
		Checks the t-values, used to determine at which part the particle is with respect to the filament, to determine the particle axis.
		Particle axis normalized from 0 to number of segments.
		Creates a set number of bins for the axis, default at 40.
		Finds the average number of particles per bin. 
		Checks bins from the start and end point of the filament. Bins that are larger than the average are filtered out.
		Stops when the bins are smaller than the average. Bins in the middle that are larger than the average are not filtered.

		!! AS OF 14.03.2018: Only properly implemented for one filament at a time !!
		"""
		# Determines the particle axis based on the t-value and segment number
		FilamentAxis_ID = []
		Seglens = np.linalg.norm(FilamentPos[1:] - FilamentPos[:-1], axis=1)
		if type(FilamentPos[0][0]) == np.float32:
			FilamentAxis = np.zeros(len(Accepted_box_part))
			Percentage = Seglens/fillen
			for i in range(len(FilamentPos)-1):
				segment_id = np.where(SegIDs[Accepted_box_part] == i)[0]
				tvals = Tsols[Accepted_box_part][segment_id]
				tvals += i
				FilamentAxis[segment_id] = tvals
			FilamentAxis_ID.append(FilamentAxis*100)
		else:
			for i in range(len(FilamentPos)):
				segmentpts = FilamentPos[i][1:] - FilamentPos[i][:-1]
				seglens = np.linalg.norm(segmentpts, axis=1)
				FilamentAxis = np.array([])
				for j in range(len(FilamentPos)-1):
					segment_id = np.where(SegIDs[Accepted_box_part[i]] == j)[0]
					Point_distance = np.linalg.norm(FilamentPos[j] - Fil_coord_closest[segment_id], axis=1)
					Percentage = np.round((Point_distance/seglens[j])*100*(j+1)).astype(np.int32)
					FilamentAxis = np.concatenate([FilamentAxis, Percentage])
				FilamentAxis_ID.append(FilamentAxis)
		FilamentAxis_ID = np.array(FilamentAxis_ID)
		
		# Bins the particle axis
		#binnum = 40
		Bins = np.linspace(0, len(Seglens), binnum)
		Bin_values = []
		if len(FilamentAxis_ID) <= 1:
			hist, bins = np.histogram(FilamentAxis_ID[0], bins=binnum)
			Bin_values.append(hist)
		else:
			for i in range(len(FilamentPos)):
				hist, bins = np.histogram(FilamentAxis_ID[i], bins=binnum)
				Bin_values.append(hist)
				
		# Currently only filters for one filament at a time! Add method for multiple filaments
		if len(Bin_values) == 1:
			checker_F = 0
			checker_L = 0
			average_num = np.average(Bin_values[0])
			Filter_F = 0
			Filter_L = len(Bin_values[0])
			for i in range(len(Bin_values[0])):
				if Bin_values[0][i] > average_num:
					Filter_F += 1
				elif Bin_values[0][-i] < average_num and not checker_F:
					checker_F = 1
					Filter_F += 1
				else:
					break
			for i in range(1,len(Bin_values[0])+1):
				if Bin_values[0][-i] > average_num:
					Filter_L -= 1
				elif Bin_values[0][-i] < average_num and not checker_L:
					checker_L = 1
					Filter_L -= 1
				else:
					break
			Filter_firstPts = np.round(bins[Filter_F:Filter_L][0])
			Filter_lastPts = np.round(bins[Filter_F:Filter_L][-1])
			Check = FilamentAxis_ID[0] > Filter_firstPts
			Check2 = FilamentAxis_ID[0] < Filter_lastPts
			Filtered_parts = Check*Check2
		return Filtered_parts

	def Do_filter_particles(self):
		""" 
		Filters the particles around the filaments. This should be run first before the thickness is defined.
		All the filtered values are saved in pickle files.
		"""
		boxexpand = 3
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(boxexpand) + 'Filtered.p'
		cachedir_ppf_distances = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/Distances/'
		cachefile_distances = cachedir_ppf_distances + Common_filename
		# Segment the particle belongs to
		cachedir_ppf_segIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/SegIDs/'
		cachefile_segIDs = cachedir_ppf_segIDs + Common_filename
		# T value that gives shortest distance, D(t) = |s(t) - p|
		cachedir_ppf_tsols = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/Tsolutions/'
		cachefile_tsols = cachedir_ppf_tsols + Common_filename
		# Particle IDs of the particle box
		cachedir_ppf_masks = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/MaskedParticlesIDs/'
		cachefile_masks = cachedir_ppf_masks + Common_filename
		# Particle position of the particle box
		cachedir_ppf_pbox = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/ParticleBoxes/'
		cachefile_pbox = cachedir_ppf_pbox + Common_filename

		if not os.path.isdir(cachedir_ppf_distances):
			os.makedirs(cachedir_ppf_distances)
		if not os.path.isdir(cachedir_ppf_segIDs):
			os.makedirs(cachedir_ppf_segIDs)
		if not os.path.isdir(cachedir_ppf_tsols):
			os.makedirs(cachedir_ppf_tsols)
		if not os.path.isdir(cachedir_ppf_masks):
			os.makedirs(cachedir_ppf_masks)
		if not os.path.isdir(cachedir_ppf_pbox):
			os.makedirs(cachedir_ppf_pbox)

		if os.path.isfile(cachefile_distances):
			Filtered_distances = pickle.load(open(cachefile_distances, 'rb'))
			Filtered_masks = pickle.load(open(cachefile_masks, 'rb'))
			Filtered_tsols = pickle.load(open(cachefile_tsols, 'rb'))
			Filtered_segids = pickle.load(open(cachefile_segIDs, 'rb'))
			Filtered_partbox = pickle.load(open(cachefile_pbox, 'rb'))		
		else:
			Filtered_distances = []
			Filtered_masks = []
			Filtered_tsols = []
			Filtered_segids = []
			Filtered_partbox = []
			for i in range(len(self.Filament_3DPos)):
				All_ids = np.array(range(len(Distances[i])))
				Part_filter = self.Filter_halo_particles(self.Filament_3DPos[i], self.Segment_index[i], self.Tsolutions[i],
														 self.FilamentLength[i], All_ids)
				Filtered_distances.append(self.Distances[Part_filter])
				Filtered_masks.append(self.Masked_IDs[Part_filter])
				Filtered_tsols.append(self.Tsolutions[Part_filter])
				Filtered_segids.append(self.Segment_index[Part_filter])
				Filtered_partbox.append(self.ParticleBoxes[Part_filter])
			Filtered_distances = np.asarray(Filtered_distances)
			Filtered_masks = np.asarray(Filtered_masks)
			Filtered_tsols = np.asarray(Filtered_tsols)
			Filtered_segids = np.asarray(Filtered_segids)
			Filtered_partbox = np.asarray(Filtered_partbox)
			pickle.dump(Filtered_distances, open(cachefile_distances, 'wb'))
			pickle.dump(Filtered_masks, open(cachefile_masks, 'wb'))
			pickle.dump(Filtered_tsols, open(cachefile_tsols, 'wb'))
			pickle.dump(Filtered_segids, open(cachefile_segIDs, 'wb'))
			pickle.dump(Filtered_partbox, open(cachefile_pbox, 'wb'))

