# Common modules
import numpy as np
import cPickle as pickle
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import time

# Own modules
import ReadGadgetFile as RGF
import OtherFunctions as OF
import ParticleMasking as PMA
import ParticlesPerFilament as PPF

# Global variables as constants
Mpc = 3.08568025e22
G_grav = 6.67258e-11
H_0 = 100*1e3/Mpc   # h/s, h = some constant usually h = 0.7
Solmass = 1.98191*1e30 # kg
rho_crit = 3.0*H_0**2/(8*np.pi*G_grav)  # kg*h^2/m^3
Npart_box_total = 512.0**3
Box_volume = (256.0*Mpc)**3.0/Npart_box_total   # m^3 per particle
DM_mass = 0.23*rho_crit*Box_volume/Solmass  # Units of solar masses * h^2 per particle
pre_rho = DM_mass*Solmass/(Mpc**3)	

class FilterParticlesAndFilaments():
	def __init__(self, model, npart, sigma):
		if model == 'lcdm':
			self.Dispersemodel = 'LCDM'
		elif model == 'symm_A':
			self.Dispersemodel = 'SymmA'
		elif model == 'symm_B':
			self.Dispersemodel = 'SymmB'
		elif model == 'symm_C':
			self.Dispersemodel = 'SymmC'
		elif model == 'symm_D':
			self.Dispersemodel = 'SymmD'
		else:
			self.Dispersemodel = model

		self.model = model
		self.npart = npart
		self.sigma = sigma
		print 'Reading filament and particle data for model: ', model
		self.Read_basic_data(model, npart, sigma)
		print 'Filtering particles'
		filter_time = time.time()
		self.Do_filter_particles()   # Filters the halo particles
		print 'Filering time:', time.time() - filter_time, 's'

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

	def Accepted_distance_density(self, particle_distances, fillen, filpos, tsols, segids):
		""" 
		Determines the thickness of the filament, based on the filament density. 
		Filament will start with a large density and decrease as the radius increases.
		The thickness is at the distance at which the average density is smaller than the threshold.
		Threshold = 10*rho_crit
		"""
		#Nparts = len(particle_distances)
		#Filtered_part = self.Filter_halo_particles(filpos, segids, tsols, fillen, np.array(range(Nparts)))
		#p_dist = particle_distances[Filtered_part]
		#distances_sorted = np.sort(p_dist)
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
	
	def Filter_filament_density_threshold(self, Fpos, particle_distances, fillen, number, tsols, segids):
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
			# First filter halo particles
			#Filtered_part = Get_filament_axis_filter_analytic(Fpos[number[i]], segids[number[i]], tsols[number[i]], fillen[i], 
			#													np.array(range(len(segids[number[i]]))))
			#p_dist = particle_distances[i][Filtered_part]
			
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
																							 self.FilamentLength[Over_IDs], range(0, len(Over_IDs)),
																							 Nf_Tsolutions, Nf_Segids)
		Return_included = Over_IDs[Included_fils] if Included_fils.any() else np.array([])
		Return_excluded = Over_IDs[Filtered_fils] if Filtered_fils.any() else np.array([])
		Return_overIDs = Over_IDs[OverThreshold] if OverThreshold.any() else np.array([])
		return Return_included, Return_excluded, Return_overIDs

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
		Small_filaments = self.FilamentLength > 1.0
		FilamentPos = self.Filament_3DPos[Small_filaments]
		Distances = self.Filtered_distances[Small_filaments]
		FilLengths = self.FilamentLength[Small_filaments]
		Tsols = self.Filtered_tsols[Small_filaments]
		SegIDs = self.Filtered_segids[Small_filaments]
		Masks = self.Filtered_masks[Small_filaments]

		Included_fils, Filtered_fils, OverThreshold = self.Filter_filament_density_threshold(FilamentPos, Distances, FilLengths, range(0, len(Distances)),
																							 Tsols, SegIDs)
		# Recomputes filaments where densities are always larger than threshold
		multiplier = 1
		while OverThreshold.any():
			multiplier += 1
			New_incFils, New_filtFils, OverThreshold = self.Recompute_densities(OverThreshold, multiplier)
			Included_fils = np.concatenate((Included_fils, New_incFils)) if New_incFils.any() else Included_fils
			Filtered_fils = np.concatenate((Filtered_fils, New_filtFils)) if New_filtFils.any() else Filtered_fils
			if multiplier > 3:
				print 'Multiplier threshold over 3, stopping the box expand computation'
				break

		Distance_thresholds = []
		for index in Included_fils:
			OK_distance = self.Accepted_distance_density(Distances[index], FilLengths[index], FilPos[index], Tsols[index], SegIDs[index])
			Distance_thresholds.append(OK_distance)
		Distance_thresholds = np.asarray(Distance_thresholds)
		# ADD PART TO RECOMPUTE THICKNESS IF 2*THICKNESS < THRESHOLD
		
		Particles_accepted = []
		Distances_accepted = []
		Accepted_box_particles = []
		for i in range(len(Distance_thresholds)):
			index = Included_fils[i]
			Filament_thickness = np.where(Distances[index] <= Distance_thresholds[i])[0]
			Accepted_box_particles.append(Filament_thickness)
			Particles_accepted.append(Masks[index][Filament_thickness])
			Distances_accepted.append(Distances[index][Filament_thickness])
		Accepted_box_particles = np.asarray(Accepted_box_particles)
		Particles_accepted = np.asarray(Particles_accepted)
		Distances_accepted = np.asarray(Distances_accepted)
		return Included_fils, Distance_thresholds, Accepted_box_particles, Particles_accepted, Distances_accepted
		
	def Get_threshold_and_noise(self):
		"""
		Reads data of distance thresholds etc. if it exist. Computes it if otherwise.
		"""
		cachedir_OKfils = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/IncludedFilaments/'
		cachedir_thresholds = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/DistanceThresholds/'
		cachedir_BoxParts = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/AcceptedBoxPartIDs/'
		cachedir_OKParts = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/AcceptedParticleIDs/'
		cachedir_OKDists = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/AcceptedDistances/'
		if not os.path.isdir(cachedir_OKfils):
			os.makedirs(cachedir_OKfils)
		if not os.path.isdir(cachedir_thresholds):
			os.makedirs(cachedir_thresholds)
		if not os.path.isdir(cachedir_BoxParts):
			os.makedirs(cachedir_BoxParts)
		if not os.path.isdir(cachedir_OKParts):
			os.makedirs(cachedir_OKParts)
		if not os.path.isdir(cachedir_OKDists):
			os.makedirs(cachedir_OKDists)

		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_okfils = cachedir_OKfils + Common_filename
		cachefile_thresholds = cachedir_thresholds + Common_filename
		cachefile_boxparts = cachedir_BoxParts + Common_filename
		cachefile_okparts = cachedir_OKParts + Common_filename
		cachefile_okdists = cachedir_OKDists + Common_filename
		if os.path.isfile(cachefile_okfils):
			print 'Reading thresholds and accepted particles from numpy files...'
			Included_fils = np.load(cachefile_okfils)
			Distance_thresholds = np.load(cachefile_thresholds)
			Accepted_box_particles = np.load(cachefile_boxparts)
			Particles_accepted = np.load(cachefile_okparts)
			Distances_accepted = np.load(cachefile_okdists)
		else:
			print 'Computing thresholds and accepted particles ...'
			computing_time = time.time()
			Included_fils, Distance_thresholds, Accepted_box_particles, Particles_accepted, Distances_accepted = self.Compute_threshold_and_noise()
			print 'Threshold computing time: ', time.time() - computing_time, 's. Now saving to numpy files'
			np.save(cachefile_okfils, Included_fils)
			np.save(cachefile_thresholds, Distance_thresholds)
			np.save(cachefile_boxparts, Accepted_box_particles)
			np.save(cachefile_okparts, Particles_accepted)
			np.save(cachefile_okdists, Distances_accepted)
		
		return Included_fils, Distance_thresholds, Accepted_box_particles, Particles_accepted, Distances_accepted

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

	def Compute_speed_components(self, filpos, Part_accept):
		segids = self.Filtered_segids
		Parallel_speeds = []
		Orthogonal_speeds = []
		All_speeds = []
		#for accepted_ids in Part_accept:
		for i in range(len(Part_accept)):
			Pvel3D = self.ParticleVel[Part_accept[i]]
			Para_speed, Orth_speed, Speed = Velocity_components_A(filpos[i], segids[Part_accept[i]], PVel3D)
			Parallel_speeds.append(Para_speed), Orthogonal_speeds.append(Orth_speed), All_speeds.append(Speed)
		return np.array(Parallel_speeds), np.array(Orthogonal_speeds), np.array(All_speeds)

	def Get_speed_components(self, filpos, Part_accept):
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
			print 'Computing speed components'
			Speed_timer = time.time()
			Parallel_speeds, Orthogonal_speeds, All_speeds = self.Compute_speed_components(filpos, Part_accept)
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
		self.filetype = filetype
		sigma_name_folder = 'Sigma'+str(Nsigma) + '/'
		foldername += sigma_name_folder

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

	def savefigure(self, figure, name):
		""" Function that calls savefig based on figure instance and filename. """
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(self.results_dir + name + self.filetype, bbox_inches='tight')

	def Particle_profiles(self, Thresholds, Accepted_parts):
		""" Plots data related to the particles """
		######## Computing relevant data
		# Computes masses of the filament. Done for all models
		self.Filament_masses = []
		for i in range(len(Accepted_parts)):
			Mass_array_temp = np.array([len(Accepted_parts[i][j]) for j in range(len(Accepted_parts[i]))])
			self.Filament_masses.append(Mass_array_temp)
		Common_bin_mass = OF.Get_common_bin_logX(self.Filament_masses, binnum=40)
		Number_mass = []
		Error_mass = []
		for i in range(len(self.Filament_masses)):
			bin_value, bin_err = OF.Get_numbers_common(self.Filament_masses, self.Filament_masses, Common_bin_mass, std='Poisson')
			Number_mass.append(bin_value)
			Error_mass.append(bin_err)

		# Relative difference of the masses. Basemodel lcdm
		RelativeDiff_mass = [OF.relative_deviation(Number_mass, i) for i in range(1, len(self.Filament_masses))]

		# Thickness of filaments as a binned histogram, using np.digitize
		Common_bin_thickness = OF.Get_common_bin_logX(Thresholds, binnum=40)
		Number_thickness = []
		for i in range(len(Thresholds)):
			#index_bin = np.digitize(Thresholds[i], Common_bin_thickness)
			#bin_value = np.array([len(Thresholds[i][index_bin == j]) for j in range(len(Common_bin_thickness))])
			bin_value, bin_std = OF.Get_numbers_common(Thresholds, Thresholds, Common_bin_thickness)
			Number_thickness.append(bin_value)

		# Com
		
		####### Plotting #######
		######## Mass histograms 
		# Mass histogram of all filaments
		NumMass_all = plt.figure()
		for i in range(len(self.Filament_masses)):
			plt.loglog(Common_bin_mass, self.Filament_masses[i], 'o-')
		plt.legend(self.All_legends)
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('$N$ filaments')
		# Mass histogram of lcdm + symmetron filaments
		NumMass_Symm = plt.figure()
		for i in range(0,5):
			plt.loglog(Common_bin_mass, self.Filament_masses[i], 'o-')
		plt.legend(self.Symm_legends)
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('$N$ filaments')
		# Mass histogram of all filaments
		NumMass_fofr = plt.figure()
		for i in [0,5,6,7]:
			plt.loglog(Common_bin_mass, self.Filament_masses[i], 'o-')
		plt.legend(self.fofr_legends)
		plt.xlabel('Filament mass - $M_\odot h^2$')
		plt.ylabel('$N$ filaments')

		######## Thickness histograms
		# Thickness histogram of all filaments
		NumThickness_all = plt.figure()
		for i in range(len(Thresholds)):
			plt.loglog(Common_bin_thickness, Number_thickness[i], 'o-')
		plt.legend(self.All_legends)
		plt.xlabel('Filament thickness - [Mpc/h]')
		plt.ylabel('$N$ filaments')
		# Thickness hisotgram of lcdm + symmetron filaments
		NumThickness_Symm = plt.figure()
		for i in range(0,5):
			plt.loglog(Common_bin_thickness, Number_thickness[i], 'o-')
		plt.legend(self.Symm_legends)
		plt.xlabel('Filament thickness - [Mpc/h]')
		plt.ylabel('$N$ filaments')
		# Thickness histogram of lcdm + f(R) filaments
		NumThickness_fofr = plt.figure()
		for i in [0,5,6,7]:
			plt.loglog(Common_bin_thickness, Number_thickness[i], 'o-')
		plt.legend(self.fofr_legends)
		plt.xlabel('Filament thickness - [Mpc/h]')
		plt.ylabel('$N$ filaments')
		
		print '--- SAVING IN: ', self.results_dir, ' ---'
		self.savefigure(NumThickness_all, 'Filament_Thickness_distribution')
		self.savefigure(NumThickness_Symm, 'Filament_Thickness_distribution_cSymmetron')
		self.savefigure(NumThickness_fofr, 'Filament_Thickness_distribution_cFofr')

	def Velocity_profiles(self, All_speeds, Orthogonal_speeds, Parallel_speeds):
		""" Plots data related to the velocity profiles """
		

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
	if args.Model:
		for models in Model_check:
			if args.Model == models:
				Model_ok = True
		if not Model_ok:
			raise ValueError('Model input name %s not correctly set -Model argument.' %args.Model)
	return args


if __name__ == '__main__':
	# Common directories
	savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
	# Parse the arguments
	parsed_arguments = Argument_parser()
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
	OK_box_particles = []
	Part_accepted = []
	Dist_accepted = []
	Models_included = []
	All_speed_list = []
	Par_speed_list = []
	Orth_speed_list = []

	def Append_data(Fid, D_thres, OK_pbox, OK_part, OK_dist, modelname):
		""" Appends common data to a common list """
		Filament_ids.append(Fid)
		Dist_thresholds.append(D_thres)
		OK_box_particles.append(OK_pbox)
		Part_accepted.append(OK_part)
		Dist_accepted.append(OK_dist)
		Models_included.append(modelname)

	def Append_data_speeds(AllS, OrthS, ParaS):
		All_speed_list.append(AllS)
		Par_speed_list.append(ParaS)
		Orth_speed_list.append(OrthS)

	if p_model == 'lcdm' or p_model == 'all':
		LCDM_instance = FilterParticlesAndFilaments('lcdm', N_parts, N_sigma)
		OK_fils_LCDM, LCDM_thresholds, OK_pbox_LDCM, OK_particles_LCDM, OK_distances_LCDM = LCDM_instance.Get_threshold_and_noise()
		Append_data(OK_fils_LCDM, LCDM_thresholds, OK_pbox_LDCM, OK_particles_LCDM, OK_distances_LCDM, 'lcdm')

	if p_model == 'symmA' or p_model == 'all':
		SA_instance = FilterParticlesAndFilaments('symm_A', N_parts, N_sigma)
		OK_fils_SA, SA_thresholds, OK_pbox_SA, OK_particles_SA, OK_distances_SA = SA_instance.Get_threshold_and_noise()
		Append_data(OK_fils_SA, SA_thresholds, OK_pbox_SA, OK_particles_SA, OK_distances_SA, 'symmA')

	if p_model == 'symmB' or p_model == 'all':
		SB_instance = FilterParticlesAndFilaments('symm_B', N_parts, N_sigma)
		OK_fils_SB, SB_thresholds, OK_pbox_SB, OK_particles_SB, OK_distances_SB = SB_instance.Get_threshold_and_noise()
		Append_data(OK_fils_SB, SB_thresholds, OK_pbox_SB, OK_particles_SB, OK_distances_SB, 'symmB')

	#if p_model == 'symmC' or p_model == 'all':
	#	SC_instance = FilterParticlesAndFilaments('symm_C', N_parts, N_sigma)
	#	OK_fils_SC, SC_thresholds, OK_pbox_SC, OK_particles_SC, OK_distances_SC = SC_instance.Get_threshold_and_noise()
	#	Append_data(OK_fils_SC, SC_thresholds, OK_pbox_SC, OK_particles_SC, OK_distances_SC, 'symmC')

	if p_model == 'symmD' or p_model == 'all':
		SD_instance = FilterParticlesAndFilaments('symm_D', N_parts, N_sigma)
		OK_fils_SD, SD_thresholds, OK_pbox_SD, OK_particles_SD, OK_distances_SD = SD_instance.Get_threshold_and_noise()
		Append_data(OK_fils_SD, SD_thresholds, OK_pbox_SD, OK_particles_SD, OK_distances_SD, 'symmD')

	if p_model == 'fofr4' or p_model == 'all':
		F4_instance = FilterParticlesAndFilaments('fofr4', N_parts, N_sigma)
		OK_fils_F4, F4_thresholds, OK_pbox_F4, OK_particles_F4, OK_distances_F4 = F4_instance.Get_threshold_and_noise()
		Append_data(OK_fils_F4, F4_thresholds, OK_pbox_F4, OK_particles_F4, OK_distances_F4, 'fofr4')

	if p_model == 'fofr5' or p_model == 'all':
		F5_instance = FilterParticlesAndFilaments('fofr5', N_parts, N_sigma)
		OK_fils_F5, F5_thresholds, OK_pbox_F5, OK_particles_F5, OK_distances_F5 = F5_instance.Get_threshold_and_noise()
		Append_data(OK_fils_F5, F5_thresholds, OK_pbox_F5, OK_particles_F5, OK_distances_F5, 'fofr5')

	if p_model == 'fofr6' or p_model == 'all':
		F6_instance = FilterParticlesAndFilaments('fofr6', N_parts, N_sigma)
		OK_fils_F6, F6_thresholds, OK_pbox_F6, OK_particles_F6, OK_distances_F6 = F6_instance.Get_threshold_and_noise()
		Append_data(OK_fils_F6, F6_thresholds, OK_pbox_F6, OK_particles_F6, OK_distances_F6, 'fofr6')

	Plot_instance = Plot_results(Models_included, N_sigma, 'ModelComparisons/ParticleAnalysis/', filetype=Filetype)
	Plot_instance.Particle_profiles(Dist_thresholds)