# Common modules
import numpy as np
import cPickle as pickle
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import time
import sys

# Own modules
import ReadGadgetFile as RGF
import OtherFunctions as OF
import ParticleMasking as PMA
import ParticlesPerFilament as PPF
import PlotFunctions as pf

# Global variables as constants
Mpc = 3.08568025e22
G_grav = 6.67258e-11
H_0 = 100*1e3/Mpc   # h/s, h = some constant usually h = 0.7
Solmass = 1.98191*1e30 # kg
rho_crit = 3.0*H_0**2/(8*np.pi*G_grav)  # kg*h^2/m^3
Npart_box_total = 512.0**3
Box_volume = (256.0*Mpc)**3.0/Npart_box_total   # m^3/h^3
DM_mass = 0.23*rho_crit*Box_volume/Solmass  # Units of M_sun/h
pre_rho = DM_mass*Solmass/(Mpc**3)

class FilterParticlesAndFilaments():
	def __init__(self, model, npart, sigma):

		self.model = model
		self.npart = npart
		self.sigma = sigma
		if model == 'lcdm':
			self.Dispersemodel = 'LCDM'
		elif model == 'symmA':
			self.model = 'symm_A'
			self.Dispersemodel = 'SymmA'
		elif model == 'symmB':
			self.model = 'symm_B'
			self.Dispersemodel = 'SymmB'
		elif model == 'symmC':
			self.model = 'symm_C'
			self.Dispersemodel = 'SymmC'
		elif model == 'symmD':
			self.model = 'symm_D'
			self.Dispersemodel = 'SymmD'
		else:
			self.Dispersemodel = model

		self.do_read_data = 0
		"""
		cachedir_OKfils = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/IncludedFilaments/'
		cachedir_speed = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/Speed/'
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_okfils = cachedir_OKfils + Common_filename
		cachefile_speed = cachedir_speed + Common_filename
		
		if not os.path.isfile(cachefile_okfils):#  or not os.path: 	# FIXME: include lack of particle speed file
			print 'Reading filament and particle data for model: ', model
			self.Read_basic_data(model, npart, sigma)
			print 'Filtering particles'
			filter_time = time.time()
			self.Do_filter_particles()   # Filters the halo particles
			print 'Filering time:', time.time() - filter_time, 's'
		"""
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
		Filament_3DPos = []
		for j in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[j], self.ydimPos[j], self.zdimPos[j])))
		self.Filament_3DPos = np.array(Filament_3DPos)
		self.FilamentLength = []
		for i in range(len(Filament_3DPos)):
			fillen = OF.Get_filament_length(Filament_3DPos[i])
			self.FilamentLength.append(fillen)
		self.FilamentLength = np.asarray(self.FilamentLength)

		# Reads particle data from Gadget
		MaskDM = np.array([0,0,0,0,0,0])*256.0
		GadgetInstance = RGF.Read_Gadget_file([0,0,0],MaskDM)
		PartPosX, PartPosY, PartPosZ, PartID, HistDM, BinXEDM, BinYEDM, KDTREE = GadgetInstance.Get_particles(model, includeKDTree=False)
		Part_velx, Part_vely, Part_velz = GadgetInstance.Get_velocities()
		
		self.Small_filaments = self.FilamentLength > 1.0   # Filter for filaments smaller than 1 Mpc/h
		
		self.ParticlePos, self.ParticleVel = OF.Get_3D(PartPosX, PartPosY, PartPosZ, Part_velx, Part_vely, Part_velz)

	def Get_filament_length(self):
		""" Only reads and returns filament lengths of a given model """
		cachedir_foldername_extra = self.Dispersemodel + 'npart'+str(self.npart)
		if self.sigma:
			cachedir_foldername_extra += 'nsig'+str(self.sigma)
		else:
			cachedir_foldername_extra += 'nsig3'
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/FilamentData/' + cachedir_foldername_extra + '/'
		
		# Pickle filenames and folder directory
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'
		# Read filament from pickle and creates 3D coordinates of filament positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(Box_boundaries, Filament_coordinate_data)
		Filament_3DPos = []
		for j in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[j], self.ydimPos[j], self.zdimPos[j])))
		self.Filament_3DPos = np.array(Filament_3DPos)
		self.FilamentLength = []
		for i in range(len(Filament_3DPos)):
			fillen = OF.Get_filament_length(Filament_3DPos[i])
			self.FilamentLength.append(fillen)
		self.FilamentLength = np.asarray(self.FilamentLength)
		self.Small_filaments = self.FilamentLength > 1.0   # Filter for filaments smaller than 1 Mpc/h
		return self.FilamentLength[self.Small_filaments]
		

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
		self.xdimPos = Fil_coord[2]
		self.ydimPos = Fil_coord[3]
		self.zdimPos = Fil_coord[4]
		self.NFilamentPoints = Fil_coord[5]
		self.FilID = Fil_coord[6]
		self.PairIDS = Fil_coord[7]

	def Read_particle_data(self, model, npart, sigma, boxexpand):
		""" Reads particle boxes, masked IDs and computed particle distances. Uses pickle files """
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

	def Read_particle_data_numpy(self, model, npart, sigma, boxexpand):
		""" Reads particle boxes, masked IDs and computed particle distances. Uses numpy save and load. """
		# Distances computed
		Common_filename =  model + '_' + str(npart) + 'part_nsig' + str(sigma)+ '_BoxExpand' + str(boxexpand) + '.npy'
		cachedir_ppf_distances = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Distances/'
		cachefile_distances = cachedir_ppf_distances + Common_filename
		Distances = np.load(cachefile_distances)
		# Segment the particle belongs to
		cachedir_ppf_segIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/SegIDs/'
		cachefile_segIDs = cachedir_ppf_segIDs  + Common_filename
		Segment_index = np.load(cachefile_segIDs)
		# T value that gives shortest distance, D(t) = |s(t) - p|
		cachedir_ppf_tsols = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Tsolutions/'
		cachefile_tsols = cachedir_ppf_tsols + Common_filename
		Tsolutions = np.load(cachefile_tsols)
		# Particle IDs of the particle box
		cachedir_ppf_masks = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/MaskedParticlesIDs/'
		cachefile_masks = cachedir_ppf_masks + Common_filename
		Masked_IDs = np.load(cachefile_masks)
		return Distances, Segment_index, Tsolutions, Masked_IDs

	def accepted_distance_func_histogram(self, distance_arr, binnum=150):
		"""
		Algorithm that determines the accepted particle distance.
		Input are computed distances to a filament
		Old method, use density method
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

	def Accepted_distance_density(self, particle_distances, fillen):
		""" 
		Determines the thickness of the filament, based on the filament density. 
		Filament will start with a large density and decrease as the radius increases.
		The thickness is at the distance at which the average density is smaller than the threshold.
		Threshold = 10*rho_crit
		"""
		distances_sorted = np.sort(particle_distances)

		Nparts2 = len(distances_sorted)
		N_smallR = np.linspace(1, Nparts2, Nparts2)
		Volumes = np.pi*fillen*distances_sorted**2
		average_densities = pre_rho*N_smallR/Volumes
		first_smaller_crit = (np.where(average_densities[20:]/rho_crit < 10.0)[0])[0]
		OK_dist = distances_sorted[first_smaller_crit]
		return OK_dist

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
		Filters the halo particles around the filaments. This should be run first before the thickness is defined.
		All the filtered values are saved in pickle files.
		"""
		boxexpand = 6
		#Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(boxexpand) + '.p'  # Pickle
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(boxexpand) + '.npy' # Numpy
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

		if not os.path.isdir(cachedir_ppf_distances):
			os.makedirs(cachedir_ppf_distances)
		if not os.path.isdir(cachedir_ppf_segIDs):
			os.makedirs(cachedir_ppf_segIDs)
		if not os.path.isdir(cachedir_ppf_tsols):
			os.makedirs(cachedir_ppf_tsols)
		if not os.path.isdir(cachedir_ppf_masks):
			os.makedirs(cachedir_ppf_masks)

		if os.path.isfile(cachefile_distances):
			# Read as Pickle file
			#self.Filtered_distances = pickle.load(open(cachefile_distances, 'rb'))
			#self.Filtered_masks = pickle.load(open(cachefile_masks, 'rb'))
			#self.Filtered_tsols = pickle.load(open(cachefile_tsols, 'rb'))
			#self.Filtered_segids = pickle.load(open(cachefile_segIDs, 'rb'))	
			# Read as numpy file
			self.Filtered_distances = np.load(cachefile_distances)
			self.Filtered_masks = np.load(cachefile_masks)
			self.Filtered_tsols = np.load(cachefile_tsols)
			self.Filtered_segids = np.load(cachefile_segIDs)
		else:
			Distances, Segment_index, Tsolutions, Masked_IDs = self.Read_particle_data_numpy(self.model, self.npart, self.sigma, 6)
			self.Filtered_distances = []
			self.Filtered_masks = []
			self.Filtered_tsols = []
			self.Filtered_segids = []
			self.Filtered_partbox = []
			filter_time = time.time()
			for i in range(len(self.Filament_3DPos)):
				All_ids = np.array(range(len(Distances[i])))
				Part_filter = self.Filter_halo_particles(self.Filament_3DPos[i], Segment_index[i], Tsolutions[i],
														 self.FilamentLength[i], All_ids)
				self.Filtered_distances.append(Distances[i][Part_filter])
				self.Filtered_masks.append(Masked_IDs[i][Part_filter])
				self.Filtered_tsols.append(Tsolutions[i][Part_filter])
				self.Filtered_segids.append(Segment_index[i][Part_filter])
			self.Filtered_distances = np.asarray(self.Filtered_distances)
			self.Filtered_masks = np.asarray(self.Filtered_masks)
			self.Filtered_tsols = np.asarray(self.Filtered_tsols)
			self.Filtered_segids = np.asarray(self.Filtered_segids)
			print 'Filering time:', time.time() - filter_time, 's'
			# Save as pickle
			#pickle.dump(self.Filtered_distances, open(cachefile_distances, 'wb'))
			#pickle.dump(self.Filtered_masks, open(cachefile_masks, 'wb'))
			#pickle.dump(self.Filtered_tsols, open(cachefile_tsols, 'wb'))
			#pickle.dump(self.Filtered_segids, open(cachefile_segIDs, 'wb'))
			#pickle.dump(self.Filtered_partbox, open(cachefile_pbox, 'wb'))
			# Save as numpy files
			np.save(cachefile_distances, self.Filtered_distances)
			np.save(cachefile_masks, self.Filtered_masks)
			np.save(cachefile_tsols, self.Filtered_tsols)
			np.save(cachefile_segIDs, self.Filtered_segids)
	
	def Filter_filament_density_threshold(self, Fpos, particle_distances, fillen):
		""" 
		Filter filaments that does not reach a given threshold density.
		Loop through the distances, in a sorted array, as a radius R. Computes the density of the filament as
		rho = M_p*N_p/V
		V = pi*L*R**2, 
		L = filament Length
		M_p, N_p = Mass of particle and number of particles respectively
		1) If density/rho_crit < 2*threshold for the first 20 points, the filament is considered noise.
		2) If max(density/rho_crit) < threshold, the filament is considered noise
		3) If min(density/rho_crit) > threshold, the filament must be recomputed with a larger box expand
		Current threshold = 10
		"""
		smallest_dist = np.array([np.min(particle_distances[j]) for j in range(len(particle_distances))])
		largest_dist = np.array([np.max(particle_distances[j]) for j in range(len(particle_distances))])
		
		larger_threshold = []
		Filtered_filaments = []
		Included_filaments = []
		
		for i in range(len(particle_distances)):
			distances_sorted = np.sort(particle_distances[i])
			average_densities = []
			pre_volume = np.pi*fillen[i]
			Noise_filament = 0
			NumParts = len(distances_sorted)
			N_smallR = np.linspace(1, NumParts, NumParts)
			Volumes = pre_volume*distances_sorted**2
			average_densities = pre_rho*N_smallR/Volumes
			First_few_densities = average_densities[0:20]/rho_crit < 20 	# Check if first few points are below a double threshold
			if First_few_densities.any():
				Noise_filament = 1
				Filtered_filaments.append(i)

			if not Noise_filament:
				if np.max(np.array(average_densities)/rho_crit) < 10:
					Filtered_filaments.append(i)
				if np.min(average_densities[1:])/rho_crit > 10:
					larger_threshold.append(i)
				else:
					Included_filaments.append(i)

		return np.array(Included_filaments), np.array(Filtered_filaments), np.array(larger_threshold)

	def Recompute_densities(self, Over_IDs, box_exp_multiplier):
		""" 
		Recomputes particle distances and masks for an increased box expand. 
		Only used for filaments whose density is always larger than the given threshold.
		"""
		SlicedParts, SlicedIDs, SliceRanges = OF.Get_particle_box_slices(self.ParticlePos)
		Nf_Distances = []
		Nf_Segids = []
		Nf_Tsolutions = []
		Nf_ParticleBoxes = []
		Nf_Masked_IDs = []
		for number in Over_IDs:
			#print self.Filament_3DPos[number]
			Dist_new, SegID_new, Tsols_new, Pbox_new, MaskIDs_new = self.Mask_and_compute_distances_again(self.Filament_3DPos[number], box_exp_multiplier*6, 
																										SlicedParts, SlicedIDs, SliceRanges)
			Nf_Distances.append(Dist_new)
			Nf_Segids.append(SegID_new)
			Nf_Tsolutions.append(Tsols_new)
			Nf_ParticleBoxes.append(Pbox_new)
			Nf_Masked_IDs.append(MaskIDs_new)
		Nf_Distances = np.array(Nf_Distances)
		Nf_Segids = np.array(Nf_Segids)
		Nf_Tsolutions = np.array(Nf_Tsolutions)
		Nf_ParticleBoxes = np.array(Nf_ParticleBoxes)
		Nf_Masked_IDs = np.array(Nf_Masked_IDs)

		Included_fils, Filtered_fils, OverThreshold = self.Filter_filament_density_threshold(self.Filament_3DPos[Over_IDs], Nf_Distances,
																							 self.FilamentLength[Over_IDs])
		Return_included = Over_IDs[Included_fils] if Included_fils.any() else np.array([])
		Return_excluded = Over_IDs[Filtered_fils] if Filtered_fils.any() else np.array([])
		Return_overIDs = Over_IDs[OverThreshold] if OverThreshold.any() else np.array([])
		return_dist = Nf_Distances[Included_fils] if Included_fils.any() else Nf_Distances
		return_masks = Nf_Masked_IDs[Included_fils] if Included_fils.any() else Nf_Masked_IDs
		return_segids = Nf_Segids[Included_fils] if Included_fils.any() else Nf_Segids
		return Return_included, Return_excluded, Return_overIDs, return_dist, return_masks, return_segids

	def Mask_and_compute_distances_again(self, FilPos, boxexp, SlicedParts, SlicedIDs, SliceRanges):
		""" 
		Recomputes particle distances and masks for an increased box expand. 
		Only used for filaments whose density is always larger than the given threshold.
		"""
		# Get indices of slices and particles to be sent
		idY, idZ = OF.get_indices_slicing(FilPos, SliceRanges, boxexp)
		Sent_particles = []
		IDs_sent = []
		for idys in idY:
			for idzs in idZ:
				Sent_particles.append(SlicedParts[idys][idzs])
				IDs_sent.append(SlicedIDs[idys][idzs])
		Sent_particles = np.array(Sent_particles)
		IDs_sent = np.array(IDs_sent)
		if len(Sent_particles) >= 2:
			Particle_indices = (np.concatenate((IDs_sent[0], IDs_sent[1])))
			Sent_particles_real = np.concatenate((Sent_particles[0], Sent_particles[1]))
			for k in range(2,len(Sent_particles)):
				Particle_indices = (np.concatenate((Particle_indices, IDs_sent[k])))
				Sent_particles_real = np.concatenate((Sent_particles_real, Sent_particles[k]))
		else:
			Particle_indices = IDs_sent[0]
			Sent_particles_real = Sent_particles[0]	

		filbox = PMA.filament_box(FilPos)
		masked_ids = PMA.masked_particle_indices(filbox, Sent_particles_real, boxexp)
		partbox_a = Sent_particles_real[masked_ids]
		masked_ids_true = Particle_indices[masked_ids]
		true_dist = []
		filaxis_temp = []
		t_solutions_all = []
		Nsegs = len(FilPos)
		for i in range(Nsegs-1):
			distances, t_sol = PPF.get_distance_analytic(FilPos[i], FilPos[i+1], partbox_a)
			true_dist.append(distances)
			t_solutions_all.append(t_sol)
		true_dist = np.array(true_dist).swapaxes(0,1)
		t_solutions_all = np.array(t_solutions_all).swapaxes(0,1)
		Shortest_dist = np.min(true_dist, axis=1)
		
		Segment_ids = true_dist.argmin(axis=1)
		Shortest_t_sol = []
		for i in range(len(Segment_ids)):
			Shortest_t_sol.append(t_solutions_all[i][Segment_ids[i]])
		return Shortest_dist, Segment_ids.astype(np.int32), np.array(Shortest_t_sol), partbox_a, masked_ids_true

	def Compute_threshold_and_noise(self):
		"""
		Determines which filament is filtered out, based on the density thresholds.
		Recomputes filaments whose density is always larger than the density thresholds.
		Stops after 3 iterations, where box_expandNew = 3*box_expand, currently box_expand = 6 
		"""
		# Filter filaments of lengths less or equal to 1 Mpc/h
		FilamentPos = self.Filament_3DPos[self.Small_filaments]
		Distances = self.Filtered_distances[self.Small_filaments]
		FilLengths = self.FilamentLength[self.Small_filaments]
		Tsols = self.Filtered_tsols[self.Small_filaments]
		SegIDs = self.Filtered_segids[self.Small_filaments]
		Masks = self.Filtered_masks[self.Small_filaments]
        
		Included_fils, Filtered_fils, OverThreshold = self.Filter_filament_density_threshold(FilamentPos, Distances, FilLengths)
		# Recomputes filaments where densities are always larger than threshold
		multiplier = 1
		while OverThreshold.any():
			multiplier += 1
			New_incFils, New_filtFils, OverThreshold, Nf_distances, Nf_masks, Nf_Segids = self.Recompute_densities(OverThreshold, multiplier)
			Included_fils = np.concatenate((Included_fils, New_incFils)) if New_incFils.any() else Included_fils
			Filtered_fils = np.concatenate((Filtered_fils, New_filtFils)) if New_filtFils.any() else Filtered_fils
			if New_incFils.any():
				Distances[New_incFils] = Nf_distances
				Masks[New_incFils] = Nf_masks
				SegIDs[New_incFils] = Nf_Segids
			if multiplier > 3:
				print 'Multiplier threshold over 3, stopping the box expand computation'
				break
		Distance_thresholds = []
		for index in Included_fils:
			OK_distance = self.Accepted_distance_density(Distances[index], FilLengths[index])
			Distance_thresholds.append(OK_distance)
		Distance_thresholds = np.asarray(Distance_thresholds)
		# ADD PART TO RECOMPUTE THICKNESS IF 2*THICKNESS < THRESHOLD 
		Particles_accepted = []
		Distances_accepted = []
		SegIDs_accepted = []
		for i in range(len(Distance_thresholds)):
			index = Included_fils[i]
			Filament_thickness = np.where(Distances[index] <= Distance_thresholds[i])[0]
			Particles_accepted.append(Masks[index][Filament_thickness])
			Distances_accepted.append(Distances[index][Filament_thickness])
			SegIDs_accepted.append(SegIDs[index][Filament_thickness])
		Particles_accepted = np.asarray(Particles_accepted)
		Distances_accepted = np.asarray(Distances_accepted)
		SegIDs_accepted = np.asarray(SegIDs_accepted)
        
		# Filter filaments with less than 100 particles
		Number_particles = np.array([len(Particles_accepted[i]) for i in range(len(Particles_accepted))])
		Few_particles = np.where(Number_particles > 100)[0]
		Included_fils = Included_fils[Few_particles]
		Distance_thresholds = Distance_thresholds[Few_particles]
		Particles_accepted = Particles_accepted[Few_particles]
		Distances_accepted = Distances_accepted[Few_particles]
		SegIDs_accepted = SegIDs_accepted[Few_particles]
		return Included_fils, Distance_thresholds, Particles_accepted, Distances_accepted, SegIDs_accepted
		
	def Get_threshold_and_noise(self):
		"""
		Reads data of distance thresholds etc. if it exist. Computes it if otherwise.
		"""
		cachedir_OKfils = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/IncludedFilaments/'
		cachedir_thresholds = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/DistanceThresholds/'
		cachedir_OKParts = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/AcceptedParticleIDs/'
		cachedir_OKDists = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/AcceptedDistances/'
		cachedir_SegIDs = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/ExpandedSegIDs/'
		if not os.path.isdir(cachedir_OKfils):
			os.makedirs(cachedir_OKfils)
		if not os.path.isdir(cachedir_thresholds):
			os.makedirs(cachedir_thresholds)
		if not os.path.isdir(cachedir_OKParts):
			os.makedirs(cachedir_OKParts)
		if not os.path.isdir(cachedir_OKDists):
			os.makedirs(cachedir_OKDists)
		if not os.path.isdir(cachedir_SegIDs):
			os.makedirs(cachedir_SegIDs)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_okfils = cachedir_OKfils + Common_filename
		cachefile_thresholds = cachedir_thresholds + Common_filename
		cachefile_okparts = cachedir_OKParts + Common_filename
		cachefile_okdists = cachedir_OKDists + Common_filename
		cachefile_segIDs = cachedir_SegIDs + Common_filename
		if os.path.isfile(cachefile_okfils):
			print 'Reading thresholds and accepted particles from numpy files...'
			Included_fils = np.load(cachefile_okfils)
			Distance_thresholds = np.load(cachefile_thresholds)
			Particles_accepted = np.load(cachefile_okparts)
			Distances_accepted = np.load(cachefile_okdists)
			SegIDs = np.load(cachefile_segIDs)
		else:
			if not self.do_read_data:
				print 'Reading filament and particle data for model: ', self.model
				self.Read_basic_data(self.model, self.npart, self.sigma)
				print 'Filtering particles'
				filter_time = time.time()
				self.Do_filter_particles()   # Filters the halo particles
				print 'Filering time:', time.time() - filter_time, 's'
				self.do_read_data = 1
			print 'Computing thresholds and accepted particles ...'
			computing_time = time.time()
			Included_fils, Distance_thresholds, Particles_accepted, Distances_accepted, SegIDs = self.Compute_threshold_and_noise()
			print 'Threshold computing time: ', time.time() - computing_time, 's. Now saving to numpy files'
			np.save(cachefile_okfils, Included_fils)
			np.save(cachefile_thresholds, Distance_thresholds)
			np.save(cachefile_okparts, Particles_accepted)
			np.save(cachefile_okdists, Distances_accepted)
			np.save(cachefile_segIDs, SegIDs)
		
		return Included_fils, Distance_thresholds, Particles_accepted, Distances_accepted, SegIDs

	def Velocity_components(self, FilamentPos, segids, PartVel_3D):
		""" 
		Computes the speeds of the different velocity components. 
		"""
		Speed = np.linalg.norm(PartVel_3D, axis=1)
		s_vectors = FilamentPos[1:] - FilamentPos[:-1]
		s_vectors_normed = np.array([stuff/np.linalg.norm(stuff) for stuff in s_vectors])
		s_vector_part = np.array([s_vectors_normed[segids[i]] for i in range(len(segids))])
		
		# Projects normal vector on the filament to the plane to the particle
		Plane_vectors = np.cross(PartVel_3D, s_vector_part)
		PlanevsSegment = np.cross(Plane_vectors, s_vector_part)
		
		Normal_vectors = PlanevsSegment/np.linalg.norm(PlanevsSegment, axis=1).reshape(len(Plane_vectors),1)
		# Computes the speeds of the parallel and orthogonal components
		Orthogonal_velocity = np.array([np.dot(PartVel_3D[i], Normal_vectors[i]) for i in range(len(PartVel_3D))])
		Parallel_velocity = np.array([np.dot(PartVel_3D[i], s_vector_part[i]) for i in range(len(PartVel_3D))])
		
		return Parallel_velocity, Orthogonal_velocity, Speed

	def Compute_speed_components(self, filpos, Part_accept, segids):
		Parallel_speeds = []
		Orthogonal_speeds = []
		All_speeds = []
		for i in range(len(Part_accept)):
			Pvel3D = self.ParticleVel[Part_accept[i]]
			Para_speed, Orth_speed, Speed = self.Velocity_components(filpos[i], segids[i], Pvel3D)
			Parallel_speeds.append(Para_speed), Orthogonal_speeds.append(Orth_speed), All_speeds.append(Speed)
		return np.array(Parallel_speeds), np.array(Orthogonal_speeds), np.array(All_speeds)

	def Get_speed_components(self, Part_accept, segids, Fils_accept):
		cachedir_speed = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/Speed/'
		cachedir_Ospeed = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/OrthogonalComp/'
		cachedir_Pspeed = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/ParallelComp/'
		if not os.path.isdir(cachedir_speed):
			os.makedirs(cachedir_speed)
		if not os.path.isdir(cachedir_Ospeed):
			os.makedirs(cachedir_Ospeed)
		if not os.path.isdir(cachedir_Pspeed):
			os.makedirs(cachedir_Pspeed)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_speed = cachedir_speed + Common_filename
		cachefile_Ospeed = cachedir_Ospeed + Common_filename
		cachefile_Pspeed = cachedir_Pspeed + Common_filename
		if os.path.isfile(cachefile_speed):
			print 'Reading from speed component numpy files...'
			All_speeds = np.load(cachefile_speed)
			Orthogonal_speeds = np.load(cachefile_Ospeed)
			Parallel_speeds = np.load(cachefile_Pspeed)
		else:
			if not self.do_read_data:
				print 'Reading filament and particle data for model: ', self.model
				self.Read_basic_data(self.model, self.npart, self.sigma)
				#print 'Filtering particles'
				#filter_time = time.time()
				#self.Do_filter_particles()
				self.do_read_data = 1
			print 'Computing speed components'
			Speed_timer = time.time()
			Parallel_speeds, Orthogonal_speeds, All_speeds = self.Compute_speed_components(self.Filament_3DPos[self.Small_filaments][Fils_accept], Part_accept, segids)
			print 'Speed computing time: ', time.time() - Speed_timer, 's'
			np.save(cachefile_speed, All_speeds)
			np.save(cachefile_Ospeed, Orthogonal_speeds)
			np.save(cachefile_Pspeed, Parallel_speeds)
		return All_speeds, Orthogonal_speeds, Parallel_speeds


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
		self.Get_legends(models)

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
			elif modnames == 'symmA':
				append_legends('Symm_A', mod=1)
			elif modnames == 'symmB':
				append_legends('Symm_B', mod=1)
			elif modnames == 'symmC':
				append_legends('Symm_C', mod=1)
			elif modnames == 'symmD':
				append_legends('Symm_D', mod=1)
			elif modnames == 'fofr4':
				append_legends('fofr4', mod=2)
			elif modnames == 'fofr5':
				append_legends('fofr5', mod=2)
			elif modnames == 'fofr6':
				append_legends('fofr6', mod=2)

	def savefigure(self, figure, name, savedir=False):
		""" Function that calls savefig based on figure instance and filename. """
		if not savedir:
			savedir_ = self.results_dir
		else:
			savedir_ = savedir
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(savedir_ + name + self.filetype, bbox_inches='tight')
        
        
	def Compute_means(self, bins_common, Part_distances, speed):
		Common_filename = 'SpeedBins_' + str(N_parts) + 'part_nsig' + str(self.Nsigma)+ '_BoxExpand' + str(6) + 'SpeedBins.npy'
		Common_filename_std = 'SpeedBins_std_' + str(N_parts) + 'part_nsig' + str(self.Nsigma)+ '_BoxExpand' + str(6) + 'SpeedBins_std.npy'
		cachedir = '/mn/stornext/u3/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/SpeedBins/'
		if not os.path.isdir(cachedir):
			os.makedirs(cachedir)
		cachefile = cachedir + Common_filename
		cachefile_std = cachedir + Common_filename_std
		if os.path.isfile(cachefile):
			Mean_speed_normalized = np.load(cachefile)
			Mean_speed_normalized_std = np.load(cachefile_std)
		else:
			Mean_speed_normalized = []
			Mean_speed_normalized_std = []
			for i in range(len(Part_distances)):
				Temp_arr = []
				Temp_arr_std = []
				for j in range(len(Part_distances[i])):
					bin_value, bin_std = OF.Bin_mean_common(Part_distances[i][j]/self.Thresholds[i][j], speed[i][j], bins_common)
					Temp_arr.append(bin_value)
					Temp_arr_std.append(bin_std)
				Mean_speed_normalized.append(np.array(Temp_arr))
				Mean_speed_normalized_std.append(np.array(Temp_arr_std))
			Mean_speed_normalized = np.asarray(Mean_speed_normalized)
			Mean_speed_normalized_std = np.asarray(Mean_speed_normalized_std)
			np.save(cachefile, Mean_speed_normalized)
			np.save(cachefile_std, Mean_speed_normalized_std)
		return Mean_speed_normalized, Mean_speed_normalized_std

	def Compute_similar_profiles(self, prop, speeds, distances, thickness, bins, botrange, uprange, prop_name, models):
		Folder = '/mn/stornext/u3/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/SimilarProperties/'
		Filename = 'Similar' + prop_name + 'Range' + str(botrange) + '-' + str(uprange) + '_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
		Filename_std = 'Similar' + prop_name + 'Range' + str(botrange) + '-' + str(uprange) + 'STD_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
		SavedFile = Folder + Filename
		SavedFile_std = Folder + Filename_std
		if not os.path.isdir(Folder):
			os.makedirs(Folder)

		if os.path.isfile(SavedFile):
			Mean_profile = np.load(SavedFile)
			Mean_profile_std = np.load(SavedFile_std)
		else:
			Similar_prop = (prop >= botrange) & (prop <= uprange)
			Speeds_included = speeds[Similar_prop]
			Parts_included = distances[Similar_prop]/thickness[Similar_prop]
			Temp_arr = []
			Temp_arr_std = []
			for k in range(len(Speeds_included)):
				bin_value, bin_std = OF.Bin_mean_common(Parts_included[k], Speeds_included[k], bins)
				Temp_arr.append(bin_value)
				Temp_arr_std.append(bin_std)
			Mean_profile, Mean_profile_std = OF.Histogram_average(bins, np.array(Temp_arr))
			np.save(SavedFile, Mean_profile)
			np.save(SavedFile_std, Mean_profile_std)
		return Mean_profile, Mean_profile_std

	def Particle_profiles(self, Thresholds, Accepted_parts, FilLengths):
		""" Plots data related to the particles """
		if len(Thresholds) != len(Accepted_parts):
			raise ValueError("Threshold data not same length as Accepted particles!")
		######## Computing relevant data ########
		s_variable = 0.7
		binnum = 40
		NModels = len(Thresholds)
		Mass_label = 'Filament mass - [$M_\odot / h$]'
		Number_label = '$N$ filaments'
		Thickness_label = 'Filament thickness - [Mpc/h]'
		Length_label = 'Filament length - [Mpc/h]'
		Mean_Mass_label = 'Mean filament mass - [$M_\odot / h$]'
		Mean_Thickness_label = 'Mean filament thickness - [Mpc/h]'
		SymmLCDM = np.array([0,1,2,3])
		FofrLCDM = np.array([0,4,5,6])
		Symm_only = np.array([1,2,3])
		Fofr_only = np.array([4,5,6])
		self.Thresholds = np.array(Thresholds)
		self.FilLengths = np.array(FilLengths)
		#### Computes masses of the filament. Done for all models
		self.Filament_masses = []
		for i in range(NModels):
			Mass_array_temp = np.array([len(Accepted_parts[i][j]) for j in range(len(Accepted_parts[i]))]).astype(np.float32)
			self.Filament_masses.append(Mass_array_temp*DM_mass)
		Common_bin_mass = OF.Get_common_bin_logX(self.Filament_masses, binnum=binnum)
		Number_mass = []
		Error_mass = []
		for i in range(NModels):
			bin_value, bin_err = OF.Bin_numbers_common(self.Filament_masses[i], self.Filament_masses[i], Common_bin_mass, std='poisson')
			Number_mass.append(bin_value)
			Error_mass.append(bin_err)
		Number_mass = np.asarray(Number_mass)
		Error_mass = np.asarray(Error_mass)
		# Relative difference of the masses. Basemodel lcdm
		RelativeDiff_mass = [OF.relative_deviation(Number_mass, i) for i in range(1, NModels)]
		Prop_error_mass = [OF.Propagate_error_reldiff(Number_mass[0], Number_mass[i], Error_mass[0], Error_mass[i]) for i in range(1, NModels)]
		# Thickness of filaments as a binned histogram, using np.digitize
		Common_bin_thickness = OF.Get_common_bin_logX(Thresholds, binnum=binnum)
		Number_thickness = []
		for i in range(NModels):
			bin_value, bin_std = OF.Bin_numbers_common(Thresholds[i], Thresholds[i], Common_bin_thickness)
			Number_thickness.append(bin_value)
		Number_thickness = np.asarray(Number_thickness)
		# Comparing different properties. E.g. length vs mass or length vs thickness etc.
		Common_bin_length = OF.Get_common_bin_logX(FilLengths, binnum=binnum)

		# Mean thickness as a function of length + relative difference
		Mean_thickness = []
		Mean_thickness_std = []
		for i in range(NModels):
			binval, bin_std = OF.Bin_mean_common(FilLengths[i], Thresholds[i], Common_bin_length)
			Mean_thickness.append(binval)
			Mean_thickness_std.append(bin_std)
		RelDiff_mean_thickness = np.array([OF.relative_deviation(Mean_thickness, i) for i in range(1, NModels)])
		Mean_thickness = np.asarray(Mean_thickness)
		# Mean mass as a function of length
		Mean_mass = []
		Mean_mass_std = []
		for i in range(NModels):
			binval, bin_std = OF.Bin_mean_common(FilLengths[i], self.Filament_masses[i], Common_bin_length)
			Mean_mass.append(binval)
			Mean_mass_std.append(bin_std)
		RelDiff_mean_mass = np.array([OF.relative_deviation(Mean_mass, i) for i in range(1, NModels)])
		Mean_mass = np.asarray(Mean_mass)
		# Mean mass as a function of thickness
		Mean_mass_vT = []
		Mean_mass_vT_std = []
		for i in range(NModels):
			binval, bin_std = OF.Bin_mean_common(Thresholds[i], self.Filament_masses[i], Common_bin_thickness)
			Mean_mass_vT.append(binval)
			Mean_mass_vT_std.append(bin_std)
		RelDiff_mean_mass_vT = np.array([OF.relative_deviation(Mean_mass_vT, i) for i in range(1, NModels)])
		Mean_mass_vT = np.asarray(Mean_mass_vT)
		Prop_err_mean_thickness = np.array([OF.Propagate_error_reldiff(Mean_thickness[0], Mean_thickness[i], Mean_thickness_std[0], Mean_thickness_std[i])
									for i in range(1,NModels)])
		Prop_err_mean_mass = np.array([OF.Propagate_error_reldiff(Mean_mass[0], Mean_mass[i], Mean_mass_std[0], Mean_mass_std[i]) for i in range(1,NModels)])
		Prop_err_mean_mass_vT = np.array([OF.Propagate_error_reldiff(Mean_mass_vT[0], Mean_mass_vT[i], Mean_mass_vT_std[0], Mean_mass_vT_std[i]) 
									for i in range(1,NModels)])
		return
		""" 
		!!!!!!!!!
		AS OF 21.03.2018, RANGES DOES NOT INCLUDE SYMMETRON C MODEL. FIX RANGES WHEN SYMMETRON C MODEL IS FIXED
		When symm C is included, ranges should be (not including LCDM):
		Symmetron: range(1, 5)  --- Currently range(0,4)
		Fofr: range(5,NModels) --- Currently range(4, NModels) 
		!!!!!!!!!
		"""
		######## Plotting ########
		######## Mass histograms 
		### Mass histogram of all filaments
		NumMass_all = pf.Call_plot_sameX(Common_bin_mass, Number_mass, Mass_label, Number_label, self.All_legends, logscale='loglog')
		### Mass histogram of lcdm + symmetron filaments
		NumMass_Symm = pf.Call_plot_sameX(Common_bin_mass, Number_mass[SymmLCDM], Mass_label, Number_label, self.Symm_legends, logscale='loglog')
		### Mass histogram of lcdm + f(R) filaments
		NumMass_fofr = pf.Call_plot_sameX(Common_bin_mass, Number_mass[FofrLCDM], Mass_label, Number_label, self.fofr_legends, logscale='loglog')
		### Mass histogram of lcdm + symmetron filaments - Semilog x scale
		NumMass_Symm_logx = pf.Call_plot_sameX(Common_bin_mass, Number_mass[SymmLCDM], Mass_label, Number_label, self.Symm_legends, logscale='logx')
		### Mass histogram of lcdm + f(R) filaments - Semilog x scale
		NumMass_fofr_logx = pf.Call_plot_sameX(Common_bin_mass, Number_mass[FofrLCDM], Mass_label, Number_label, self.fofr_legends, logscale='logx')
		### Mass histograms with errors, lcdm + symmetron
		NumMass_error_symm = plt.figure()
		for i in SymmLCDM:
			plt.plot(Common_bin_mass, Number_mass[i], alpha=0.7)
			plt.fill_between(Common_bin_mass, Number_mass[i]-Error_mass[i], Number_mass[i]+Error_mass[i], alpha=0.3)
		plt.legend(self.Symm_legends)
		plt.xlabel('Filament mass - $M_\odot / h$')
		plt.ylabel('$N$ filaments')
		plt.xscale('log')
		### Mass histograms with errors, lcdm + f(R)
		NumMass_error_fofr = plt.figure()
		for i in FofrLCDM:
			plt.plot(Common_bin_mass, Number_mass[i], alpha=0.7)
			plt.fill_between(Common_bin_mass, Number_mass[i]-Error_mass[i], Number_mass[i]+Error_mass[i], alpha=0.3)
		plt.legend(self.fofr_legends)
		plt.xlabel('Filament mass - $M_\odot/h$')
		plt.ylabel('$N$ filaments')
		plt.xscale('log')

		### Relative differene of all models
		RelDiff_mass_all = plt.figure()
		plt.semilogx(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		for i in range(NModels-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
		plt.legend(self.All_legends)
		plt.xlabel('Filament mass - $M_\odot / h$')
		plt.ylabel('Relative difference of $N$ filament')
		### Relative difference of lcdm + symmetron
		RelDiff_mass_Symm = plt.figure()
		#plt.semilogx(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		for i in range(3):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
		plt.legend(self.Symm_legends[1:])
		plt.xlabel('Filament mass - $M_\odot / h$')
		plt.ylabel('Relative difference of $N$ filament')
		### Relative difference of lcdm + f(R)
		RelDiff_mass_fofr = plt.figure()
		#plt.semilogx(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		for i in range(3, NModels-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
		plt.legend(self.fofr_legends[1:])
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('Relative difference of $N$ filament')
		### Relative difference of lcdm + symmetron, with error
		RelDiff_mass_Symm_err = plt.figure()
		#plt.plot(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		#plt.fill_between(Common_bin_mass, np.zeros(len(Common_bin_mass)), np.zeros(len(Common_bin_mass)))
		for i in (Symm_only-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
			plt.fill_between(Common_bin_mass, RelativeDiff_mass[i]-Prop_error_mass[i], RelativeDiff_mass[i]+Prop_error_mass[i], alpha=0.3)
		plt.legend(self.Symm_legends[1:])
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('Relative difference of filament masses')
		plt.xscale('log')
		RelDiff_mass_fofr_err = plt.figure()
		#plt.plot(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		#plt.fill_between(Common_bin_mass, np.zeros(len(Common_bin_mass)), np.zeros(len(Common_bin_mass)))
		for i in (Fofr_only-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
			plt.fill_between(Common_bin_mass, RelativeDiff_mass[i]-Prop_error_mass[i], RelativeDiff_mass[i]+Prop_error_mass[i], alpha=0.3)
		plt.legend(self.fofr_legends[1:])
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('Relative difference of filament masses')
		plt.xscale('log')
		#plt.yscale('log')

		######## Thickness histograms
		### Thickness histogram of all filaments
		NumThickness_all = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness, Thickness_label, Number_label, self.All_legends, logscale='loglog')
		### Thickness hisotgram of lcdm + symmetron filaments, including logX scale
		NumThickness_Symm = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[SymmLCDM], Thickness_label, Number_label, self.Symm_legends, logscale='loglog')
		NumThickness_Symm_logX = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[SymmLCDM], Thickness_label, Number_label, self.Symm_legends, logscale='logx')
		### Thickness histogram of lcdm + f(R) filaments, including LogX scale
		NumThickness_fofr = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[FofrLCDM], Thickness_label, Number_label, self.fofr_legends, logscale='loglog')
		NumThickness_fofr_logX = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[FofrLCDM], Thickness_label, Number_label, self.fofr_legends, logscale='logx')
		
		######## Compare different properties
		### Thickness as a function of length, LCDM + Symmetron and LCDM + f(R)
		ThickVsLen_Symm = pf.Call_plot_sameX(Common_bin_length, Mean_thickness[SymmLCDM], Length_label, Mean_Thickness_label, self.Symm_legends, logscale='loglog')
		ThickVsLen_fofr = pf.Call_plot_sameX(Common_bin_length, Mean_thickness[FofrLCDM], Length_label, Mean_Thickness_label, self.fofr_legends, logscale='loglog')
		### Relative difference with errobar, Symmetron and f(R) seperate, base model = LCDM
		ThickVsLen_RelErr_Symm = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Symm_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_thickness[i], alpha=0.7)
			plt.fill_between(Common_bin_length, RelDiff_mean_thickness[i]-Prop_err_mean_thickness[i], RelDiff_mean_thickness[i]+Prop_err_mean_thickness[i], alpha=0.4)
		plt.legend(self.Symm_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('Relative difference of thickness')
		plt.xscale('log')

		ThickVsLen_RelErr_fofr = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Fofr_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_thickness[i], alpha=0.7)
			plt.fill_between(Common_bin_length, RelDiff_mean_thickness[i]-Prop_err_mean_thickness[i], RelDiff_mean_thickness[i]+Prop_err_mean_thickness[i], alpha=0.4)
		plt.legend(self.fofr_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('Relative difference of thickness')
		plt.xscale('log')

		### Mass as a function of length, LCDM + Symmetron and LCDM + f(R)
		MassVsLen_Symm = pf.Call_plot_sameX(Common_bin_length, Mean_mass[SymmLCDM], Length_label, Mean_Mass_label, self.Symm_legends, logscale='loglog')
		MassVsLen_fofr = pf.Call_plot_sameX(Common_bin_length, Mean_mass[FofrLCDM], Length_label, Mean_Mass_label, self.fofr_legends, logscale='loglog')
		### Relative difference with errorbar, Symmetron and f(R) seperate, base model = LCDM
		MassVsLen_RelErr_Symm = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Symm_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_mass[i], alpha=0.7)
			plt.fill_between(Common_bin_length, RelDiff_mean_mass[i]-Prop_err_mean_mass[i], RelDiff_mean_mass[i]+Prop_err_mean_mass[i], alpha=0.4)
		plt.legend(self.Symm_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('Relative difference of mass')
		plt.xscale('log')
		MassVsLen_RelErr_fofr = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Fofr_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_mass[i], alpha=0.7)
			plt.fill_between(Common_bin_length, RelDiff_mean_mass[i]-Prop_err_mean_mass[i], RelDiff_mean_mass[i]+Prop_err_mean_mass[i], alpha=0.4)
		plt.legend(self.fofr_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('Relative difference of mass')
		plt.xscale('log')

		### Mass as a function of thickness, LCDM + Symmetron and LCDM + f(R)
		MassVsThick = plt.figure(figsize=(8,4))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		plt.subplot(1,2,1)
		for i in SymmLCDM:
			plt.loglog(Common_bin_thickness, Mean_mass_vT[i])
		plt.legend(self.Symm_legends)
		plt.subplot(1,2,2)
		for i in FofrLCDM:
			plt.loglog(Common_bin_thickness, Mean_mass_vT[i])
		plt.legend(self.fofr_legends)
		MassVsThick.text(0.5, 0.01, Thickness_label, ha='center', fontsize=10)
		MassVsThick.text(0.04, 0.6, Mass_label, ha='center', rotation='vertical', fontsize=10)
		#plt.tight_layout()
						
		print '--- SAVING IN: ', self.results_dir, ' ---'
		######## Mass histograms 
		self.savefigure(NumMass_all, 'Filament_mass_distribution')
		self.savefigure(NumMass_Symm, 'Filament_mass_distribution_cSymmetron')
		self.savefigure(NumMass_fofr, 'Filament_mass_distribution_cFofr')
		self.savefigure(NumMass_Symm_logx, 'Filament_mass_distribution_logX_cSymmetron')
		self.savefigure(NumMass_fofr_logx, 'Filament_mass_distribution_logX_cFofr')
		self.savefigure(NumMass_error_symm, 'Filament_mass_distribution_error_cSymmetron')
		self.savefigure(NumMass_error_fofr , 'Filament_mass_distribution_error_cFofr')
		self.savefigure(RelDiff_mass_all, 'Relative_difference_mass')
		self.savefigure(RelDiff_mass_Symm, 'Relative_difference_mass_cSymmetron')
		self.savefigure(RelDiff_mass_fofr, 'Relative_difference_mass_cFofr')
		self.savefigure(RelDiff_mass_Symm_err, 'Relative_difference_mass_error_cSymmetron')
		######## Thickness histograms
		self.savefigure(NumThickness_all, 'Filament_Thickness_distribution')
		self.savefigure(NumThickness_Symm, 'Filament_Thickness_distribution_cSymmetron')
		self.savefigure(NumThickness_fofr, 'Filament_Thickness_distribution_cFofr')
		self.savefigure(NumThickness_Symm_logX, 'Filament_Thickness_distribution_logX_cSymmetron')
		self.savefigure(NumThickness_fofr_logX, 'Filament_Thickness_distribution_logX_cFofr')
		######## Compare different properties
		self.savefigure(ThickVsLen_Symm, 'ThicknessVsLength_cSymmetron')
		self.savefigure(ThickVsLen_fofr, 'ThicknessVsLength_cFofr')
		self.savefigure(ThickVsLen_RelErr_Symm, 'ThicknessVsLength_Reldiff_cSymmetron')
		self.savefigure(ThickVsLen_RelErr_fofr, 'ThicknessVsLength_Reldiff_cFofr')
		self.savefigure(MassVsLen_Symm, 'MassVsLength_cSymmetron')
		self.savefigure(MassVsLen_fofr, 'MassVsLength_cFofr')
		self.savefigure(MassVsLen_RelErr_Symm, 'MassVsLength_Reldiff_cSymmetron')
		self.savefigure(MassVsLen_RelErr_fofr, 'MassVsLength_Reldiff_cFofr')		
		self.savefigure(MassVsThick, 'MassVsThickness')
		plt.close('all')	# Clear all current windows to free memory

	def Velocity_profiles(self, All_speeds, Orthogonal_speeds, Parallel_speeds, Part_distances):
		""" Plots data related to the velocity profiles """
		velocity_savefile_dir = 'ModelComparisons/VelocityAnalysis/'
		if self.raw_filetype == 'png':
			velocity_savefile_dir += 'PNG/'
		elif self.raw_filetype == 'pdf':
			velocity_savefile_dir += 'PDF/'
		#self.filetype = '.' + self.raw_filetype
		sigma_name_folder = 'Sigma'+str(self.Nsigma) + '/'
		velocity_savefile_dir += sigma_name_folder

		velocity_results_dir = os.path.join(savefile_directory, velocity_savefile_dir)
		if not os.path.isdir(velocity_results_dir):
			os.makedirs(velocity_results_dir)

		######## Compute relevant stuff ########
		binnum = 40
		NModels = len(All_speeds)
		s_variable = 0.7
		Mass_label = 'Filament mass - [$M_\odot / h$]'
		Number_label = '$N$ filaments'
		Thickness_label = 'Filament thickness - [Mpc/h]'
		Length_label = 'Filament length - [Mpc/h]'
		Mean_Mass_label = 'Mean filament mass - [$M_\odot / h$]'
		Mean_Thickness_label = 'Mean filament thickness - [Mpc/h]'
		Average_speed_label =  r'$\langle v \rangle$ - [km/s]'
		SymmLCDM = np.array([0,1,2,3])
		FofrLCDM = np.array([0,4,5,6])
		Symm_only = np.array([1,2,3])
		Fofr_only = np.array([4,5,6])
		All_speeds = np.asarray(All_speeds)
		Orthogonal_speeds = np.asarray(Orthogonal_speeds)
		Parallel_speeds = np.asarray(Parallel_speeds)
		Symm_filenames = ['LCDM', 'Symm_A', 'Symm_B', 'Symm_D']
		Fofr_filenames = ['LCDM', 'fofr4', 'fofr5', 'fofr6']
		#### Average speed of all filaments as a function of 
		Normalized_part_distances = []
		Distance_limits = []
		Distance_limits_normalized = []
		for i in range(NModels):
			Normalized_part_distances.append(Part_distances[i]/self.Thresholds[i])
			New_min = np.array([np.min(Part_distances[i][j]) for j in range(len(Part_distances[i]))])
			New_max = np.array([np.max(Part_distances[i][j]) for j in range(len(Part_distances[i]))])
			New_min_norm = np.array([np.min(Normalized_part_distances[-1][j]) for j in range(len(Normalized_part_distances[-1]))])
			New_max_norm = np.array([np.max(Normalized_part_distances[-1][j]) for j in range(len(Normalized_part_distances[-1]))])
			Distance_limits.append(np.concatenate((New_min, New_max)))
			Distance_limits_normalized.append(np.concatenate((New_min_norm, New_max_norm)))
		Normalized_part_distances = np.asarray(Normalized_part_distances)
		Common_bin_thickness = OF.Get_common_bin_logX(self.Thresholds, binnum=binnum)
		Common_bin_distances = OF.Get_common_bin(Distance_limits, binnum=binnum)
		Common_bin_distances_normalized = OF.Get_common_bin(Distance_limits_normalized, binnum=binnum)        
		Common_bin_mass = OF.Get_common_bin_logX(self.Filament_masses, binnum=binnum)
		print 'Got common bins'
		Mean_speed_normalized, Mean_speed_normalized_std = self.Compute_means(Common_bin_distances_normalized, Part_distances, All_speeds)
		#### Compute the average of the average speed of each filament. That is, average profile for --all-- filaments
		Mean_speed_allFils = []
		Mean_speed_allFils_std = []
		for i in range(NModels):
			#Mean, std = OF.Mean_per_bin(Common_bin_distances_normalized, Mean_speed_normalized[i])
			Mean, std = OF.Histogram_average(Common_bin_distances_normalized, Mean_speed_normalized[i])
			Mean_speed_allFils.append(Mean)
			Mean_speed_allFils_std.append(std)
		print 'Mean speeds done computing'
		### Compute the average of particle speeds for filament with similar masses

		######## Plotting data ########
		#### Average particle speed for all filaments
		Mass_titles = ['$M \in [10^{12}, 10^{13}]M_\odot$', '$M \in [10^{13}, 10^{14}]M_\odot$', '$M \in [10^{14}, 10^{15}]M_\odot$']
		AverageSpeed_AllFils = plt.figure(figsize=(12,5))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		plt.subplot(1,2,1)
		for i in SymmLCDM:
			plt.plot(Common_bin_distances_normalized, Mean_speed_allFils[i])
		plt.legend(self.Symm_legends)
		plt.subplot(1,2,2)
		for i in FofrLCDM:
			plt.plot(Common_bin_distances_normalized, Mean_speed_allFils[i])
		plt.legend(self.fofr_legends)
		AverageSpeed_AllFils.text(0.5, 0.01, 'Particle distance to filament center (normalized)', ha='center', fontsize=10)
		AverageSpeed_AllFils.text(0.04, 0.6, Average_speed_label, ha='center', rotation='vertical', fontsize=10)

		### Average speed of filaments with similar masses, comparing LCDM + Symmetron
		AverageSpeed_SimilarMass_Symm = plt.figure(figsize=(30,8))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		Mass_ranges = [1e12, 1e13, 1e14, 1e15]   # Units of M_sun/h, maybe use min and max of mass bin?
		for j in range(len(Mass_ranges)-1):
			ax = plt.subplot(1,3,j+1)
			print 'iteration ',j 
			for i in SymmLCDM:
				Mean_profile, Mean_profile_std = self.Compute_similar_profiles(self.Filament_masses[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
																			Common_bin_distances_normalized, Mass_ranges[j], Mass_ranges[j+1], 'Mass', Symm_filenames[i])
				plt.plot(Common_bin_distances_normalized, Mean_profile, label=self.Symm_legends[i])
				#plt.fill_between(Common_bin_distances_normalized, Mean_profile-Mean_profile_std, Mean_profile+Mean_profile_std, alpha=0.3)
			plt.title(Mass_titles[j], fontsize=10)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		AverageSpeed_SimilarMass_Symm.text(0.5, 0.01, 'Particle distance to filament center (normalized)', ha='center', fontsize=10)
		AverageSpeed_SimilarMass_Symm.text(0.04, 0.7, Average_speed_label, ha='center', rotation='vertical', fontsize=10)
		### Average speed of filaments with similar masses, comparing LCDM + f(R)
		AverageSpeed_SimilarMass_fofr = plt.figure(figsize=(20,5))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for j in range(len(Mass_ranges)-1):
			ax = plt.subplot(1,3,j+1)
			print 'iteration ',j 
			for ij in range(len(FofrLCDM)):
				i = FofrLCDM[ij]
				Mean_profile, Mean_profile_std = self.Compute_similar_profiles(self.Filament_masses[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
																			Common_bin_distances_normalized, Mass_ranges[j], Mass_ranges[j+1], 'Mass', Fofr_filenames[ij])
				plt.plot(Common_bin_distances_normalized, Mean_profile, label=self.fofr_legends[ij])
				#plt.fill_between(Common_bin_distances_normalized, Mean_profile-Mean_profile_std, Mean_profile+Mean_profile_std, alpha=0.3)
			plt.title(Mass_titles[j], fontsize=10)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		AverageSpeed_SimilarMass_fofr.text(0.5, 0.01, 'Particle distance to filament center (normalized)', ha='center', fontsize=10)
		AverageSpeed_SimilarMass_fofr.text(0.04, 0.7, Average_speed_label, ha='center', rotation='vertical', fontsize=10)
		### Average speed of filament of similar length. Symmetron + LCDM comparison
		Length_titles = ['$L \in [1,5]$ Mpc/h', '$L \in [5,10]$ Mpc/h', '$L \in [10,20]$ Mpc/h']
		AverageSpeed_SimilarLength_Symm = plt.figure(figsize=(20,5))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		Length_ranges = [1, 5, 10, 20]   # Units of Mpc/h, maybe use min and max of mass bin?
		for j in range(len(Length_ranges)-1):
			ax = plt.subplot(1,3,j+1)
			print 'iteration ',j
			for i in SymmLCDM:
				Mean_profile, Mean_profile_std = self.Compute_similar_profiles(self.FilLengths[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
																			Common_bin_distances_normalized, Length_ranges[j], Length_ranges[j+1], 'Length', Symm_filenames[i])
				plt.plot(Common_bin_distances_normalized, Mean_profile, label=self.Symm_legends[i])
				#plt.fill_between(Common_bin_distances_normalized, Mean_profile-Mean_profile_std, Mean_profile+Mean_profile_std, alpha=0.3)
			plt.title(Length_titles[j], fontsize=10)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		AverageSpeed_SimilarLength_Symm.text(0.5, 0.01, 'Particle distance to filament center (normalized)', ha='center', fontsize=10)
		AverageSpeed_SimilarLength_Symm.text(0.04, 0.7, Average_speed_label, ha='center', rotation='vertical', fontsize=10)
		#### Average speed of filament of similar length. f(R) + LCDM comparison
		AverageSpeed_SimilarLength_fofr = plt.figure(figsize=(20,5))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		Length_ranges = [1, 5, 10, 20]   # Units of Mpc/h, maybe use min and max of mass bin?
		for j in range(len(Length_ranges)-1):
			ax = plt.subplot(1,3,j+1)
			print 'iteration ',j
			for ij in range(len(FofrLCDM)):
				i = FofrLCDM[ij]
				Mean_profile, Mean_profile_std = self.Compute_similar_profiles(self.FilLengths[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
																			Common_bin_distances_normalized, Length_ranges[j], Length_ranges[j+1], 'Length', Fofr_filenames[ij])
				plt.plot(Common_bin_distances_normalized, Mean_profile, label=self.fofr_legends[ij])
				#plt.fill_between(Common_bin_distances_normalized, Mean_profile-Mean_profile_std, Mean_profile+Mean_profile_std, alpha=0.3)
			plt.title(Length_titles[j], fontsize=10)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		AverageSpeed_SimilarLength_fofr.text(0.5, 0.01, 'Particle distance to filament center (normalized)', ha='center', fontsize=10)
		AverageSpeed_SimilarLength_fofr.text(0.04, 0.7, Average_speed_label, ha='center', rotation='vertical', fontsize=10)
		#### Average speed for different mass bins, with LCDM and Symmetron models
		AverageSpeed_MassBins = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		plt.subplot(1,2,1)
		for k in SymmLCDM:
			Average_speed = []
			Error_speed = []
			for i in range(len(Common_bin_mass)-1):
				Similar_masses = (self.Filament_masses[k] >= Common_bin_mass[i]) & (self.Filament_masses[k] <= Common_bin_mass[i+1])
				Speeds_included = All_speeds[k][Similar_masses]
				Average_per_fil = np.array([np.average(Speeds_included[j]) for j in range(len(Speeds_included))])
				Standard_deviation = np.array([np.std(Average_per_fil)])/np.sqrt(len(Speeds_included))
				Error_speed.append(Standard_deviation)
				Average_speed.append(np.average(Average_per_fil))
			#plt.plot(Common_bin_mass[1:], Average_speed, 'x-')
			plt.errorbar(Common_bin_mass[1:], Average_speed, Error_speed)
		plt.legend(self.Symm_legends)
		plt.subplot(1,2,2)
		for k in FofrLCDM:
			Average_speed = []
			Error_speed = []
			for i in range(len(Common_bin_mass)-1):
				Similar_masses = (self.Filament_masses[k] >= Common_bin_mass[i]) & (self.Filament_masses[k] <= Common_bin_mass[i+1])
				Speeds_included = All_speeds[k][Similar_masses]
				Average_per_fil = np.array([np.average(Speeds_included[j]) for j in range(len(Speeds_included))])
				Standard_deviation = np.array([np.std(Average_per_fil)])/np.sqrt(len(Speeds_included))
				Error_speed.append(Standard_deviation)
				Average_speed.append(np.average(Average_per_fil))
			plt.errorbar(Common_bin_mass[1:], Average_speed, Error_speed)
			#plt.plot(Common_bin_mass[1:], Average_speed, 'x-')
		plt.legend(self.fofr_legends)
		AverageSpeed_MassBins.text(0.5, 0.01, Mass_label, ha='center', fontsize=10)
		AverageSpeed_MassBins.text(0.04, 0.7, Average_speed_label, ha='center', rotation='vertical', fontsize=10)
		
		####### Save figures #######
		print '--- SAVING IN: ', velocity_results_dir, ' ---'
		self.savefigure(AverageSpeed_AllFils, 'Average_speed_all_filaments', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarMass_Symm, 'Average_speed_similar_mass_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarMass_fofr, 'Average_speed_similar_mass_cFofr', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarLength_Symm, 'Average_speed_similar_length_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarLength_fofr, 'Average_speed_similar_length_cFofr', velocity_results_dir)
		self.savefigure(AverageSpeed_MassBins, 'Average_speed_massbins', velocity_results_dir)
		plt.close('all')	# Clear all current windows to free memory

	#def Similar_profiles(self, All_speeds, ):


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
	# Common directories
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
	Filament_ids = []
	Dist_thresholds = []
	Part_accepted = []
	Dist_accepted = []
	Models_included = []

	Filament_lengths = []
	All_speed_list = []
	Par_speed_list = []
	Orth_speed_list = []

	def Append_data(Fid, D_thres,  OK_part, OK_dist, modelname):
		""" Appends common data to a common list """
		Filament_ids.append(Fid)
		Dist_thresholds.append(D_thres)
		Part_accepted.append(OK_part)
		Dist_accepted.append(OK_dist)
		Models_included.append(modelname)

	def Append_data_speeds(AllS, OrthS, ParaS):
		All_speed_list.append(AllS)
		Par_speed_list.append(ParaS)
		Orth_speed_list.append(OrthS)

	for p_model in Models_to_be_run:
		if p_model == 'symmC':
			continue
		else:
			Instance = FilterParticlesAndFilaments(p_model, N_parts, N_sigma)
			OK_fils, thresholds, OK_particles, OK_distances, SegIDs = Instance.Get_threshold_and_noise()
			Filament_lengths.append(Instance.Get_filament_length()[OK_fils])
			Append_data(OK_fils, thresholds, OK_particles, OK_distances, p_model)
			Speeds, Ospeed, Pspeed = Instance.Get_speed_components(OK_particles, SegIDs, OK_fils)
			Append_data_speeds(Speeds, Ospeed, Pspeed)
	Plot_instance = Plot_results(Models_included, N_sigma, 'ModelComparisons/ParticleAnalysis/', filetype=Filetype)
	Plot_instance.Particle_profiles(Dist_thresholds, Part_accepted, Filament_lengths)
	Plot_instance.Velocity_profiles(All_speed_list, Orth_speed_list, Par_speed_list, Dist_accepted)