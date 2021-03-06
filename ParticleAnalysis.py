# Common modules
import numpy as np
import cPickle as pickle
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.switch_backend('agg')
import argparse
import time
import sys
import subprocess

# ZeroMQ
import zmq
import ZMQArraySending as ZMQAS

# Own modules
import ReadGadgetFile as RGF
import OtherFunctions as OF
import ParticleMasking as PMA
import ParticlesPerFilament as PPF
import PlotFunctions as pf
import Histogram_comparison as HComp

# Global variables as constants
s_variable = 0.7  # FOR PLOTTING
Omega_m0 = 0.267
Mpc = 3.08568025e22 		# m
G_grav = 6.67258e-11
h_param = 1   # 0.72
H_0 = h_param*100*1e3/Mpc   # 1/s, h = some constant usually h = 0.7
H_tilde = 100.0*1e3/Mpc 	# Units of h/s
Solmass = 1.98191*1e30 		# kg
rho_crit = 3.0*H_tilde**2/(8*np.pi*G_grav)  # kg*h^2/m^3
Npart_box_total = 512.0**3
Box_volume = (256.0*Mpc)**3 		# Units of m^3/h^3
DM_mass = (3.0*H_tilde**2*Omega_m0*Box_volume/(8*np.pi*G_grav*Npart_box_total))/Solmass 	# Units of M_sun/h
pre_rho = DM_mass*Solmass/(Mpc**3) 		# kg h^2/m^3

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
		self.Unpack_filament_data(Box_boundaries, Filament_coordinate_data, CP_coordinate_data)
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
		
		# Reads particle data from Gadget
		if model == 'symm_C':
			self.ParticlePos, self.ParticleVel = RGF.read_file_ascii_inclVel('/mn/stornext/d5/aleh/SymmC_data/SymmC_particles.solve')
		else:
			MaskDM = np.array([0,0,0,0,0,0])*256.0
			GadgetInstance = RGF.Read_Gadget_file([0,0,0],MaskDM)
			PartPosX, PartPosY, PartPosZ, PartID, HistDM, BinXEDM, BinYEDM, KDTREE = GadgetInstance.Get_particles(model, includeKDTree=False)
			Part_velx, Part_vely, Part_velz = GadgetInstance.Get_velocities()
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
		cachedir_nanvalues = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/NanValues/'
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_nanvalues = cachedir_nanvalues + Common_filename
		if os.path.isfile(cachefile_nanvalues):
			RemoveNanFils = np.load(cachefile_nanvalues)
		else:
			RemoveNanFils = np.array([])
				
		# Read filament from pickle and creates 3D coordinates of filament positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(Box_boundaries, Filament_coordinate_data, CP_coordinate_data)
		Filament_3DPos = []
		for j in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[j], self.ydimPos[j], self.zdimPos[j])))
		Filament_3DPos = np.array(Filament_3DPos)
		print "Number of filaments for " + self.Dispersemodel + ":", len(Filament_3DPos)
		FilamentLength = []
		for i in range(len(Filament_3DPos)):
			fillen = OF.Get_filament_length(Filament_3DPos[i])
			FilamentLength.append(fillen)
		FilamentLength = np.asarray(FilamentLength)
		self.FilamentLength = FilamentLength
		if RemoveNanFils.any():
			FilamentLength = np.delete(FilamentLength, RemoveNanFils)
		Small_filaments = FilamentLength > 1.0   # Filter for filaments smaller than 1 Mpc/h
		Only_small_filaments = FilamentLength <= 1.0
		print 'All filaments: ', len(FilamentLength)
		print 'Number of short filaments: ', len(FilamentLength[Only_small_filaments])
		return FilamentLength[Small_filaments]
		

	def Unpack_filament_data(self, Box_info, Fil_coord, CP_coord):
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

		self.CritPointXpos = CP_coord[0]
		self.CritPointYpos = CP_coord[1]
		self.CritPointZpos = CP_coord[2]
		self.CP_type = CP_coord[3]
		self.CP_persistent_pair = CP_coord[4] 
		self.Critpts_filamentID = CP_coord[5] 
		self.CP_id_of_connecting_filament = CP_coord[6] 
		self.Number_filaments_connecting_to_CP = CP_coord[7]

	def Number_filament_connections(self):
		""" Computes the number of filament connections one filament has"""
		#print 'Computing number connections'
		NumFilamentConnections = []
		for Connected_CP in self.PairIDS:
			counter = 0
			for CP_ids in Connected_CP:
				counter += self.Number_filaments_connecting_to_CP[CP_ids]
			NumFilamentConnections.append(counter)

		NumFilamentConnections = np.asarray(NumFilamentConnections)
		return NumFilamentConnections

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
		# Nanvalues encountered, indices of filament
		cachedir_nanvalues = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/NanValues/'
		cachefile_nanvalues = cachedir_nanvalues + Common_filename
		if not os.path.isdir(cachedir_ppf_distances):
			os.makedirs(cachedir_ppf_distances)
		if not os.path.isdir(cachedir_ppf_segIDs):
			os.makedirs(cachedir_ppf_segIDs)
		if not os.path.isdir(cachedir_ppf_tsols):
			os.makedirs(cachedir_ppf_tsols)
		if not os.path.isdir(cachedir_ppf_masks):
			os.makedirs(cachedir_ppf_masks)
		if not os.path.isdir(cachedir_nanvalues):
			os.makedirs(cachedir_nanvalues)
		

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
			if os.path.isfile(cachefile_nanvalues):
				self.RemoveNanFils = np.load(cachefile_nanvalues)
			else:
				self.RemoveNanFils = np.array([])
		else:
			Distances, Segment_index, Tsolutions, Masked_IDs = self.Read_particle_data_numpy(self.model, self.npart, self.sigma, boxexpand)
			self.Filtered_distances = []
			self.Filtered_masks = []
			self.Filtered_tsols = []
			self.Filtered_segids = []
			self.Filtered_partbox = []
			self.RemoveNanFils = []
			filter_time = time.time()
			for i in range(len(self.Filament_3DPos)):
				All_ids = np.array(range(len(Distances[i])))
				if not np.isnan(Tsolutions[i]).any():
					Part_filter = self.Filter_halo_particles(self.Filament_3DPos[i], Segment_index[i], Tsolutions[i],
															 self.FilamentLength[i], All_ids)
					self.Filtered_distances.append(Distances[i][Part_filter])
					self.Filtered_masks.append(Masked_IDs[i][Part_filter])
					self.Filtered_tsols.append(Tsolutions[i][Part_filter])
					self.Filtered_segids.append(Segment_index[i][Part_filter])
				else:
					print "Nan values encountered for filament, ", i
					self.RemoveNanFils.append(i)
			self.Filtered_distances = np.asarray(self.Filtered_distances)
			self.Filtered_masks = np.asarray(self.Filtered_masks)
			self.Filtered_tsols = np.asarray(self.Filtered_tsols)
			self.Filtered_segids = np.asarray(self.Filtered_segids)
			self.RemoveNanFils = np.asarray(self.RemoveNanFils)
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
			np.save(cachefile_nanvalues, self.RemoveNanFils)
	
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
		#smallest_dist = np.array([np.min(particle_distances[j]) for j in range(len(particle_distances))])
		#largest_dist = np.array([np.max(particle_distances[j]) for j in range(len(particle_distances))])
		
		larger_threshold = []
		Filtered_filaments = []
		Included_filaments = []
		
		for i in range(len(particle_distances)):
			distances_sorted = np.sort(particle_distances[i])
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
			Dist_new, SegID_new, Tsols_new, Pbox_new, MaskIDs_new = self.Mask_and_compute_distances_again(self.Filament_3DPos[self.Small_filaments][number],
																									box_exp_multiplier*6, SlicedParts, SlicedIDs, SliceRanges)
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

		Included_fils, Filtered_fils, OverThreshold = self.Filter_filament_density_threshold(self.Filament_3DPos[self.Small_filaments][Over_IDs], Nf_Distances,
																							 self.FilamentLength[self.Small_filaments][Over_IDs])
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
		if self.RemoveNanFils.any():
			self.Filament_3DPos = np.delete(self.Filament_3DPos, self.RemoveNanFils)
			self.FilamentLength = np.delete(self.FilamentLength, self.RemoveNanFils)
			self.Small_filaments = self.FilamentLength > 1.0
			FilamentPos = self.Filament_3DPos[self.Small_filaments]
			Distances = self.Filtered_distances[self.Small_filaments]
			FilLengths = self.FilamentLength[self.Small_filaments]
			Tsols = self.Filtered_tsols[self.Small_filaments]
			SegIDs = self.Filtered_segids[self.Small_filaments]
			Masks = self.Filtered_masks[self.Small_filaments]	
		else:
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
			print 'Some filaments over the threshold. Recomputing with increased boxsize...'
			multiplier += 1
			New_incFils, New_filtFils, OverThreshold, Nf_distances, Nf_masks, Nf_Segids = self.Recompute_densities(OverThreshold, multiplier)
			Included_fils = np.concatenate((Included_fils, New_incFils)) if New_incFils.any() else Included_fils
			Filtered_fils = np.concatenate((Filtered_fils, New_filtFils)) if New_filtFils.any() else Filtered_fils
			if New_incFils.any():
				for j in range(len(New_incFils)):
					index = New_incFils[j]
					Distances[index] = Nf_distances[j]
					Masks[index] = Nf_masks[j]
					SegIDs[index] = Nf_Segids[j]
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
		cachedir_nanvalues = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/FilteredParts/NanValues/'
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
		cachefile_nanvalues = cachedir_nanvalues + Common_filename
		read_secondtime = 0
		if os.path.isfile(cachefile_speed):
			print 'Reading from speed component numpy files...'
			All_speeds = np.load(cachefile_speed)
			Orthogonal_speeds = np.load(cachefile_Ospeed)
			Parallel_speeds = np.load(cachefile_Pspeed)
		else:
			if not self.do_read_data:
				print 'Reading filament and particle data for model: ', self.model
				self.Read_basic_data(self.model, self.npart, self.sigma)
				if os.path.isfile(cachefile_nanvalues):
					self.RemoveNanFils = np.load(cachefile_nanvalues)
				else:
					print 'Warning: Remove nan file does not exist for this model!'
					self.RemoveNanFils = np.array([])
				#print 'Filtering particles'
				#filter_time = time.time()
				#self.Do_filter_particles()
				read_secondtime = 1
				self.do_read_data = 1
			if self.RemoveNanFils.any() and read_secondtime:
				self.Filament_3DPos = np.delete(self.Filament_3DPos, self.RemoveNanFils)
				self.FilamentLength = np.delete(self.FilamentLength, self.RemoveNanFils)
				self.Small_filaments = self.FilamentLength > 1.0
			print 'Computing speed components'
			Speed_timer = time.time()
			Parallel_speeds, Orthogonal_speeds, All_speeds = self.Compute_speed_components(self.Filament_3DPos[self.Small_filaments][Fils_accept], Part_accept, segids)
			print 'Speed computing time: ', time.time() - Speed_timer, 's'
			np.save(cachefile_speed, All_speeds)
			np.save(cachefile_Ospeed, Orthogonal_speeds)
			np.save(cachefile_Pspeed, Parallel_speeds)
		return All_speeds, Orthogonal_speeds, Parallel_speeds

	def Compute_density_profile(self, distances, fillen):
		""" Computes density of the filament based on the distance input """
		cachedir_density = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/DensityProfiles/'
		if not os.path.isdir(cachedir_density):
			os.makedirs(cachedir_density)
		
		Common_filename =  self.model + '_' + str(self.npart) + 'part_nsig' + str(self.sigma)+ '_BoxExpand' + str(6) + '.npy'
		cachefile_density = cachedir_density + Common_filename
		if os.path.isfile(cachefile_density):
			Density_profiles = np.load(cachefile_density)
		else:		
			Density_profiles = []
			for i in range(len(distances)):
				distances_sorted = np.sort(distances[i])
				pre_volume = np.pi*fillen[i]
				NumParts = len(distances_sorted)
				N_smallR = np.linspace(1, NumParts, NumParts)
				Volumes = pre_volume*distances_sorted**2 		# Mpc^3/h^3
				average_densities = pre_rho*N_smallR/Volumes    #(kg h^2 / m^3 )/ (Mpc^3 /h^3) = kg/(m^3 h)
				Density_profiles.append(average_densities)
			Density_profiles = np.asarray(Density_profiles)
			np.save(cachefile_density, Density_profiles)
		return Density_profiles

	def Secondary_filter(self, Accepted_ids):
		""" 
		Does a secondary filter on the filtered filaments.
		Can be used to, for instance, filter out long filaments.
		"""		
		#self.Small_filaments = self.FilamentLength > 1.0   # Filter for filaments smaller than 1 Mpc/h
		FilamentLength = self.Get_filament_length()
		#Accepted_lengths = self.FilamentLength[self.Small_filaments][Accepted_ids]
		Accepted_lengths = FilamentLength[Accepted_ids]
		Large_filaments = Accepted_lengths <= 20 	# Units of Mpc/h. Change accordingly here
		return Large_filaments

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
		self.Linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]
		self.SymmLCDM_ids = np.array([])
		self.FofrLCDM_ids = np.array([])
		self.Get_legends(models)
		self.SymmLCDM_ids = self.SymmLCDM_ids.astype(np.int32)
		self.FofrLCDM_ids = self.FofrLCDM_ids.astype(np.int32)
		self.Dist_mass_filter = 0
        
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

	def savefigure(self, figure, name, savedir=False, dpi_mult=1):
		""" Function that calls savefig based on figure instance and filename. """
		if not savedir:
			savedir_ = self.results_dir
		else:
			savedir_ = savedir
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(savedir_ + name + self.filetype, bbox_inches='tight', dpi=100*dpi_mult)
        
        
	def Compute_means(self, bins_common, Part_distances, speed):
		Common_filename = self.Speed_filename + 'SpeedBins_' + str(N_parts) + 'part_nsig' + str(self.Nsigma)+ '_BoxExpand' + str(6) + 'SpeedBins.npy'
		Common_filename_std =self.Speed_filename + 'SpeedBins_std_' + str(N_parts) + 'part_nsig' + str(self.Nsigma)+ '_BoxExpand' + str(6) + 'SpeedBins_std.npy'
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

	def Compute_similar_profiles(self, prop, speeds, distances, thickness, bins, botrange, uprange, prop_name, models, newbinning=40):
		if DoFilter:
			Folder = '/mn/stornext/u3/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/SimilarProperties/FurtherFilter/'
		else:
			Folder = '/mn/stornext/u3/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/ProcessedData/SpeedComponents/SimilarProperties/'
		# Further filter folder used when filtering masses, thickness and long filaments
		if newbinning == 40:
			Filename = self.Speed_filename + 'Similar' + prop_name + 'Range' + str(botrange) + '-' + str(uprange) + '_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
			Filename_std = self.Speed_filename + 'Similar' + prop_name + 'Range' + str(botrange) + '-' + str(uprange) + 'STD_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
		else:
			Filename = self.Speed_filename + 'Similar' + prop_name + 'bins' + str(newbinning) + 'Range' + str(botrange) + '-' + str(uprange) \
										+ '_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
			Filename_std = self.Speed_filename + 'Similar' + prop_name + 'bins' + str(newbinning) + 'Range' + str(botrange) + '-' + str(uprange) \
										+ 'STD_Npart' + str(N_parts) + '_Nsigma' + str(self.Nsigma) + models + '.npy'
			
		SavedFile = Folder + Filename
		SavedFile_std = Folder + Filename_std
		if not os.path.isdir(Folder):
			os.makedirs(Folder)

		if os.path.isfile(SavedFile):
			Mean_profile = np.load(SavedFile)
			Mean_profile_std = np.load(SavedFile_std)
		else:
			print 'Computing similar property of: ', prop_name
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
			#Mean_profile, Mean_profile_std = self.Compute_similar_profiles_ZMQ(bins, Parts_included, Speeds_included)
			np.save(SavedFile, Mean_profile)
			np.save(SavedFile_std, Mean_profile_std)
		return Mean_profile, Mean_profile_std

	def Compute_similar_profiles_ZMQ(self, bins, Parts_included, Speeds_included):
		context = zmq.Context()
		context.linger = 0
		# Socket to send messages on
		sender = context.socket(zmq.PUSH)
		#sender.bind("tcp://127.0.0.1:5050")
		sender.bind("tcp://*:4070")

		# Socket where the data is received from
		data_receive = context.socket(zmq.PULL)
		#data_receive.bind("tcp://127.0.0.1:5052")
		data_receive.bind("tcp://*:4072")

		# Socket where end message is sent to. Used to tell workers the jobs are finished
		control_sender = context.socket(zmq.PUSH)
		#control_sender.bind("tcp://127.0.0.1:5054")
		control_sender.bind("tcp://*:4074")

		# Poller, used to check whether stuff is done or not
		poller = zmq.Poller()
		poller.register(data_receive, zmq.POLLIN)

		# Calls the script that starts up a set amount of workers.
		# Stops program a little bit to let the workers start up
		print "Starting processes"
		subprocess.call("./SpawnWorkers_speedprofile.sh 40", shell=True)
		
		time.sleep(5)
		time_dist = time.time()
		# Sends data
		print "Done, sending data"
		NumFils = len(Speeds_included)
		for i in range(NumFils):
			ZMQAS.send_zipped_pickle(sender, [Parts_included[i], Speeds_included[i], bins])	# Analytical method
			if i == NumFils-1:
				print("all data sent")
		# Give time to send data
		time.sleep(5)
		Temp_arr = []
		#print 'Looping through data receiving'
		for j in range(NumFils):
			socks = dict(poller.poll())
			if socks.get(data_receive) == zmq.POLLIN:
				data = ZMQAS.recv_zipped_pickle(data_receive)
				Temp_arr.append(data[0])
			if not socks:
				print("All data received?")
				break
		for i in range(41):
			control_sender.send_string("FINISHED")
		Mean_profile, Mean_profile_std = OF.Histogram_average(bins, np.array(Temp_arr))
		sender.close()
		data_receive.close()
		control_sender.close()
		context.term()
		return Mean_profile, Mean_profile_std

	def Compute_average_speeds_Propertybinned(self, speeds, prop, bins):
		Average_speed = []
		Error_speed = []
		for i in range(len(bins)-1):
			Similar_prop = (prop >= bins[i]) & (prop <= bins[i+1])
			Speeds_included = speeds[Similar_prop]
			Average_per_fil = np.array([np.average(Speeds_included[j]) for j in range(len(Speeds_included))])
			Standard_deviation = np.nanstd(Average_per_fil)/np.sqrt(len(Speeds_included))
			Average_speed.append(np.average(Average_per_fil))
			Error_speed.append(Standard_deviation)
		return np.array(Average_speed), np.array(Error_speed)

	def Particle_profiles(self, Thresholds, Accepted_parts, FilLengths):
		""" Plots data related to the particles """
		if len(Thresholds) != len(Accepted_parts):
			raise ValueError("Threshold data not same length as Accepted particles!")
		######## Computing relevant data ########
		binnum = 20
		NModels = len(Thresholds)
		do_fb = True
		#Mass_label = 'Filament mass - [$M_\odot / h$]'
		#Number_label = '$N$ filaments'
		#Thickness_label = 'Filament thickness - [Mpc/h]'
		#Length_label = 'Filament length - [Mpc/h]'
		#Mean_Mass_label = 'Mean filament mass - [$M_\odot / h$]'
		#Mean_Thickness_label = 'Mean filament thickness - [Mpc/h]'
		Mass_label = '$M$ - [$M_\odot / h$]'
		Number_label = '$N$'
		Number_label_reldiff = '$N_i/N_{\Lambda \mathrm{CDM}} - 1$'#'$(N_i - N_{\Lambda \mathrm{CDM}})/N_{\Lambda \mathrm{CDM}}$'
		Thickness_label = '$R_T$ - [Mpc/$h$]'
		Thickness_label_reldiff = '$R_{T,i}/R_{T,{\Lambda \mathrm{CDM}}} - 1$'#'$(T_i - T_{\Lambda\mathrm{CDM}})/T_{\Lambda \mathrm{CDM}}$'
		Length_label = '$L$ - [Mpc/$h$]'
		Mean_Mass_label = r'$\langle M \rangle - [M_\odot/h]$'
		Mean_Mass_label_reldiff = r'$\langle M_i \rangle / \langle M_{\Lambda\mathrm{CDM}} \rangle - 1$'
		Mean_Thickness_label = r'$\langle R_T \rangle - [\mathrm{Mpc}/h]$'
		Mean_Thickness_label_reldiff = r'$\langle R_{T,i} \rangle/\langle R_{T,{\Lambda\mathrm{CDM}}} \rangle - 1$'
		Density_label = r'$\langle \rho \rangle - [10^{12} M_\odot h^2/\mathrm{Mpc}^3]$'
		#Density_label_reldiff = r'$(\langle \rho_i \rangle - \langle \rho_{\Lambda \mathrm{CDM}} \rangle)/\langle \rho_{\Lambda \mathrm{CDM}} \rangle$'
		Density_label_reldiff = r'$\langle \rho_i \rangle/\langle \rho_{\Lambda \mathrm{CDM}} - 1$'
		#SymmLCDM = np.array([0,1,2,3,4])
		#FofrLCDM = np.array([0,5,6,7])
		#Symm_only = np.array([1,2,3,4])
		#Fofr_only = np.array([5,6,7])
		SymmLCDM = self.SymmLCDM_ids
		FofrLCDM = self.FofrLCDM_ids
		Symm_only = self.SymmLCDM_ids[1:]
		Fofr_only = self.FofrLCDM_ids[1:]
		
		#### Computes masses of the filament. Done for all models
		self.Filament_masses = []
		for i in range(NModels):
			Mass_array_temp = np.array([len(Accepted_parts[i][j]) for j in range(len(Accepted_parts[i]))]).astype(np.float32)
			self.Filament_masses.append(Mass_array_temp*DM_mass)
		# Further filter out based on masses
		self.Mass_filters = []
		if DoFilter:
			for i in range(NModels):
				Mass_filter = self.Filament_masses[i] < 1e15
				self.Mass_filters.append(Mass_filter)
				self.Filament_masses[i] = self.Filament_masses[i][Mass_filter]
				Thresholds[i] = Thresholds[i][Mass_filter]
				FilLengths[i] = FilLengths[i][Mass_filter]
		self.Thresholds = np.array(Thresholds)
		self.FilLengths = np.array(FilLengths)
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
		RelativeDiff_mass = np.array([OF.relative_deviation(Number_mass, i) for i in range(1, NModels)])
		Prop_error_mass = np.array([OF.Propagate_error_reldiff(Number_mass[0], Number_mass[i], Error_mass[0], Error_mass[i]) for i in range(1, NModels)])

		### Thickness of filaments as a binned histogram, using np.digitize
		Common_bin_thickness = OF.Get_common_bin_logX(Thresholds, binnum=binnum)
		Number_thickness = []
		Error_thickness = []
		for i in range(NModels):
			bin_value, bin_std = OF.Bin_numbers_common(Thresholds[i], Thresholds[i], Common_bin_thickness, std='poisson')
			bin_value[1:] = np.nan_to_num(bin_value[1:])
			Number_thickness.append(bin_value)
			Error_thickness.append(bin_std)
		Number_thickness = np.asarray(Number_thickness)
		# Relative difference of the thickness. Basemodel lcdm
		RelativeDiff_thickness = np.array([OF.relative_deviation(Number_thickness, i) for i in range(1, NModels)])
		for i in range(len(RelativeDiff_thickness)):
			RelativeDiff_thickness[i][RelativeDiff_thickness[i] == np.inf] = 0.0
			RelativeDiff_thickness[i][1:] = np.nan_to_num(RelativeDiff_thickness[i][1:])
		Prop_error_thickness = np.array([OF.Propagate_error_reldiff(Number_thickness[0], Number_thickness[i], Error_thickness[0], Error_thickness[i]) 
										for i in range(1, NModels)])

		### The average density of the filaments at radius = filament thickness
		self.Filament_density = []
		for i in range(NModels):
			Volumes = np.pi*FilLengths[i]*Thresholds[i]**2
			rho = self.Filament_masses[i]/Volumes
			self.Filament_density.append(rho/1.0e12)   # Divide by 10^12 to scale better logarithmically
		self.Filament_density = np.asarray(self.Filament_density)
		Common_bin_density = OF.Get_common_bin_logX(self.Filament_density, binnum=binnum)
		# Compute number of filament within a density bins
		Number_density = []
		Error_density = []
		for i in range(NModels):
			bin_value, bin_std = OF.Bin_numbers_common(self.Filament_density[i], self.Filament_density[i], Common_bin_density, std='poisson')
			Number_density.append(bin_value)
			Error_density.append(bin_std)
		Number_density = np.asarray(Number_density)
		Error_density = np.asarray(Error_density)
		# Relative difference of the densities. Bademodel lcdm
		RelativeDiff_density = np.array([OF.relative_deviation(Number_density, i) for i in range(1, NModels)])
		Prop_error_density = np.array([OF.Propagate_error_reldiff(Number_density[0], Number_density[i], Error_density[0], Error_density[i])
									for i in range(1, NModels)])

		### Comparing different properties. E.g. length vs mass or length vs thickness etc.
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
		Mean_thickness_std = np.asarray(Mean_thickness_std)
		# Mean mass as a function of length
		Mean_mass = []
		Mean_mass_std = []
		for i in range(NModels):
			binval, bin_std = OF.Bin_mean_common(FilLengths[i], self.Filament_masses[i], Common_bin_length)
			Mean_mass.append(binval)
			Mean_mass_std.append(bin_std)
		RelDiff_mean_mass = np.array([OF.relative_deviation(Mean_mass, i) for i in range(1, NModels)])
		Mean_mass = np.asarray(Mean_mass)
		Mean_mass_std = np.asarray(Mean_mass_std)
		# Mean mass as a function of thickness
		Mean_mass_vT = []
		Mean_mass_vT_std = []
		for i in range(NModels):
			binval, bin_std = OF.Bin_mean_common(Thresholds[i], self.Filament_masses[i], Common_bin_thickness)
			Mean_mass_vT.append(binval)
			Mean_mass_vT_std.append(bin_std)
		RelDiff_mean_mass_vT = np.array([OF.relative_deviation(Mean_mass_vT, i) for i in range(1, NModels)])
		Mean_mass_vT = np.asarray(Mean_mass_vT)
		Mean_mass_vT_std = np.asarray(Mean_mass_vT_std)

		# Computing propagated errors of the relative differences for the three properties above
		Prop_err_mean_thickness = np.array([OF.Propagate_error_reldiff(Mean_thickness[0], Mean_thickness[i], Mean_thickness_std[0], Mean_thickness_std[i])
									for i in range(1,NModels)])
		Prop_err_mean_mass = np.array([OF.Propagate_error_reldiff(Mean_mass[0], Mean_mass[i], Mean_mass_std[0], Mean_mass_std[i]) for i in range(1,NModels)])
		Prop_err_mean_mass_vT = np.array([OF.Propagate_error_reldiff(Mean_mass_vT[0], Mean_mass_vT[i], Mean_mass_vT_std[0], Mean_mass_vT_std[i]) 
									for i in range(1,NModels)])
		#sys.exit(1)
		#return
		""" 
		!!!!!!!!!
		AS OF 21.03.2018, RANGES DOES NOT INCLUDE SYMMETRON C MODEL. FIX RANGES WHEN SYMMETRON C MODEL IS FIXED
		When symm C is included, ranges should be (not including LCDM):
		Symmetron: range(1, 5)  --- Currently range(0,4)
		Fofr: range(5,NModels) --- Currently range(4, NModels) 
		AS OF 08.04.2018 --- Above problem should be fixed
		!!!!!!!!!
		"""
		######## Plotting ########
		######## Mass histograms 
		### Mass histogram of all filaments
		xlim_mass = (Common_bin_mass[0], Common_bin_mass[-1])
		NumMass_all = pf.Call_plot_sameX_OLD(Common_bin_mass, Number_mass, Mass_label, Number_label, self.All_legends, logscale='loglog')
		### Mass histogram of lcdm + symmetron filaments
		NumMass_Symm = pf.Call_plot_sameX(Common_bin_mass, Number_mass[SymmLCDM], Mass_label, Number_label, self.Symm_legends, 
										 self.Plot_colors_symm, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		### Mass histogram of lcdm + f(R) filaments
		NumMass_fofr = pf.Call_plot_sameX(Common_bin_mass, Number_mass[FofrLCDM], Mass_label, Number_label, self.fofr_legends,
										 self.Plot_colors_fofr, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		### Mass histogram of lcdm + symmetron filaments - Semilog x scale
		NumMass_Symm_logx = pf.Call_plot_sameX(Common_bin_mass, Number_mass[SymmLCDM]/1000.0, Mass_label, r'$N \times 1000$', self.Symm_legends,
										 self.Plot_colors_symm, xscale='log', linestyles=self.Linestyles, legend_anchor=False, xlim=xlim_mass)
		### Mass histogram of lcdm + f(R) filaments - Semilog x scale
		NumMass_fofr_logx = pf.Call_plot_sameX(Common_bin_mass, Number_mass[FofrLCDM]/1000.0, Mass_label, r'$N \times 1000$', self.fofr_legends,
										 self.Plot_colors_fofr, xscale='log', linestyles=self.Linestyles, legend_anchor=False, xlim=xlim_mass)


		### Mass histograms with errors, lcdm + symmetron
		NumMass_error_symm = plt.figure()
		for i in SymmLCDM:
			plt.plot(Common_bin_mass, Number_mass[i], alpha=0.7, color=self.Plot_colors_symm[i])
			plt.fill_between(Common_bin_mass, Number_mass[i]-Error_mass[i], Number_mass[i]+Error_mass[i], alpha=0.3, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xlabel(Mass_label)
		plt.ylabel(Number_label)
		plt.xscale('log')
		### Mass histograms with errors, lcdm + f(R)
		NumMass_error_fofr = plt.figure()
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Common_bin_mass, Number_mass[i], alpha=0.7, color=self.Plot_colors_fofr[ij])
			plt.fill_between(Common_bin_mass, Number_mass[i]-Error_mass[i], Number_mass[i]+Error_mass[i], alpha=0.3, facecolor=self.Plot_colors_fofr[ij])
		plt.legend(self.fofr_legends)
		plt.xlabel(Mass_label)
		plt.ylabel(Number_label)
		plt.xscale('log')

		### Relative differene of all models
		RelDiff_mass_all = plt.figure()
		plt.semilogx(Common_bin_mass, np.zeros(len(Common_bin_mass)))
		for i in range(NModels-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i])
		plt.legend(self.All_legends)
		plt.xlabel(Mass_label)
		plt.ylabel('$(N_i - N_{\Lambda CDM})/N_{\Lambda CDM}$')
		### Relative difference of lcdm + symmetron
		RelDiff_mass_Symm = plt.figure()
		for i in Symm_only-1:
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i], color=self.Plot_colors_symm[i+1])
		plt.legend(self.Symm_legends[1:])
		plt.xlabel(Mass_label)
		plt.ylabel('$(N_i - N_{\Lambda CDM})/N_{\Lambda CDM}$')
		### Relative difference of lcdm + f(R)
		RelDiff_mass_fofr = plt.figure()
		for ij in range(len(Fofr_only)):
			i = Fofr_only[ij]-1
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i], color=self.Plot_colors_fofr[ij+1])
		plt.legend(self.fofr_legends[1:])
		plt.xlabel(Mass_label)
		plt.ylabel('$(N_i - N_{\Lambda CDM})/N_{\Lambda CDM}$')

		### Relative difference of lcdm + symmetron, with error
		RelDiff_mass_Symm_err = plt.figure()
		plt.plot(Common_bin_mass[1:], np.zeros(len(Common_bin_mass[1:])), color='k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in (Symm_only-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i], color=self.Plot_colors_symm[i+1], linestyle=self.Linestyles[i])
			plt.fill_between(Common_bin_mass, RelativeDiff_mass[i]-Prop_error_mass[i], RelativeDiff_mass[i]+Prop_error_mass[i],
							 alpha=0.3, facecolor=self.Plot_colors_symm[i+1])
		plt.legend(self.Symm_legends)
		plt.xlabel(Mass_label)
		plt.ylabel('$(M_i - M_{\Lambda CDM})/M_{\Lambda CDM}$')
		plt.xscale('log')
		plt.ylim(-1,1)
		plt.xlim(xlim_mass)
		### Relative difference of lcdm + symmetron, with error
		RelDiff_mass_fofr_err = plt.figure()
		plt.plot(Common_bin_mass[1:], np.zeros(len(Common_bin_mass[1:])), color='k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in (Fofr_only-1):
			plt.semilogx(Common_bin_mass, RelativeDiff_mass[i], color=self.Plot_colors_fofr[i-3], linestyle=self.Linestyles[i-4])
			plt.fill_between(Common_bin_mass, RelativeDiff_mass[i]-Prop_error_mass[i], RelativeDiff_mass[i]+Prop_error_mass[i],
							 alpha=0.3, facecolor=self.Plot_colors_fofr[i-3])
		plt.legend(self.fofr_legends)
		plt.xlabel(Mass_label)
		plt.ylabel('$(M_i - M_{\Lambda CDM})/M_{\Lambda CDM}$')
		plt.xscale('log')
		plt.ylim(-1,1)
		plt.xlim(xlim_mass)
		#plt.yscale('log')

		######## Thickness histograms
		### Thickness histogram of all filaments
		xlim_thickness = (Common_bin_thickness[0], Common_bin_thickness[-1])
		NumThickness_all = pf.Call_plot_sameX_OLD(Common_bin_thickness, Number_thickness, Thickness_label, Number_label, self.All_legends, logscale='loglog')
		### Thickness hisotgram of lcdm + symmetron filaments, including logX scale
		NumThickness_Symm = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[SymmLCDM], Thickness_label, Number_label, self.Symm_legends,
												self.Plot_colors_symm, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		NumThickness_Symm_logX = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[SymmLCDM], Thickness_label, Number_label, self.Symm_legends,
												self.Plot_colors_symm, xscale='log', linestyles=self.Linestyles, legend_anchor=False,
												xlim=xlim_thickness)
		### Thickness histogram of lcdm + f(R) filaments, including LogX scale
		NumThickness_fofr = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[FofrLCDM], Thickness_label, Number_label, self.fofr_legends,
												self.Plot_colors_fofr, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		NumThickness_fofr_logX = pf.Call_plot_sameX(Common_bin_thickness, Number_thickness[FofrLCDM], Thickness_label, Number_label, self.fofr_legends,
												self.Plot_colors_fofr, xscale='log', linestyles=self.Linestyles, legend_anchor=False,
												xlim=xlim_thickness)
		
		### Relative differences of thickness
		Reldiff_num_thick_Symm = pf.Call_plot_sameX(Common_bin_thickness, RelativeDiff_thickness[Symm_only-1], Thickness_label, Number_label_reldiff,
											self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Prop_error_thickness[Symm_only-1], xscale='log', yscale='log',
											linestyles=self.Linestyles, reldiff=True, xlim=xlim_thickness,
											legend_anchor=False)
		Reldiff_num_thick_Symm_logX = pf.Call_plot_sameX(Common_bin_thickness, RelativeDiff_thickness[Symm_only-1], Thickness_label, Number_label_reldiff,
											self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Prop_error_thickness[Symm_only-1], xscale='log',
											linestyles=self.Linestyles, fillbetween=do_fb, reldiff=True, legend_anchor=False)
		
		Reldiff_num_thick_fofr = pf.Call_plot_sameX(Common_bin_thickness, RelativeDiff_thickness[Fofr_only-1], Thickness_label, Number_label_reldiff,
											self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Prop_error_thickness[Fofr_only-1], xscale='log', yscale='log',
											linestyles=self.Linestyles, reldiff=True, xlim=xlim_thickness,
											legend_anchor=False)
		Reldiff_num_thick_fofr_logX = pf.Call_plot_sameX(Common_bin_thickness, RelativeDiff_thickness[Fofr_only-1], Thickness_label, Number_label_reldiff,
											self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Prop_error_thickness[Fofr_only-1],  xscale='log',
											linestyles=self.Linestyles, fillbetween=do_fb, reldiff=True, legend_anchor=False)
		
		######## Density histograms
		### Density histograms with number of filaments at a given density bin
		NumDensity_Symm = pf.Call_plot_sameX(Common_bin_density, Number_density[SymmLCDM], Density_label, Number_label, self.Symm_legends,
											self.Plot_colors_symm, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		NumDensity_Symm_logX = pf.Call_plot_sameX(Common_bin_density, Number_density[SymmLCDM], Density_label, Number_label, self.Symm_legends,
											self.Plot_colors_symm, xscale='log', linestyles=self.Linestyles, legend_anchor=False)
		NumDensity_fofr = pf.Call_plot_sameX(Common_bin_density, Number_density[FofrLCDM], Density_label, Number_label, self.fofr_legends,
											self.Plot_colors_fofr, xscale='log', yscale='log', linestyles=self.Linestyles, legend_anchor=False)
		NumDensity_fofr_logX = pf.Call_plot_sameX(Common_bin_density, Number_density[FofrLCDM], Density_label, Number_label, self.fofr_legends,
											self.Plot_colors_fofr, xscale='log', linestyles=self.Linestyles, legend_anchor=False)
		### Relative difference of the above
		RelDiff_num_density_Symm = pf.Call_plot_sameX(Common_bin_density, RelativeDiff_density[Symm_only-1], Density_label, Number_label_reldiff,
											self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Prop_error_density[Symm_only-1], xscale='log', yscale='log',
											linestyles=self.Linestyles, reldiff=True, legend_anchor=False)
		RelDiff_num_density_Symm_logX = pf.Call_plot_sameX(Common_bin_density, RelativeDiff_density[Symm_only-1], Density_label, Number_label_reldiff,
											self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Prop_error_density[Symm_only-1], xscale='log',
											linestyles=self.Linestyles, reldiff=True, ylim=(-1,1), xlim=(Common_bin_density[0], Common_bin_density[-1]), 
											legend_anchor=False)
		RelDiff_num_density_fofr = pf.Call_plot_sameX(Common_bin_density, RelativeDiff_density[Fofr_only-1], Density_label, Number_label_reldiff,
											self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Prop_error_density[Fofr_only-1], xscale='log', yscale='log',
											linestyles=self.Linestyles, reldiff=True, legend_anchor=False)
		RelDiff_num_density_fofr_logX = pf.Call_plot_sameX(Common_bin_density, RelativeDiff_density[Fofr_only-1], Density_label, Number_label_reldiff,
											self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Prop_error_density[Fofr_only-1], xscale='log',
											linestyles=self.Linestyles, reldiff=True, ylim=(-1,1), xlim=(Common_bin_density[0], Common_bin_density[-1]), 
											legend_anchor=False)


		######## Compare different properties
		def Do_bin_speed_gridspec(xdata, ydata1, ydata2, err1, err2, colors, legend1, linestyle1, xlabel, ylabel1, ylabel2, ylim1=(0,0), ylim2=(0,0), rowcol=[2,2]):
			figure = plt.figure()
			gs = gridspec.GridSpec(rowcol[0],rowcol[1])
			for i in range(len(ydata1)):
				ax0 = plt.subplot(gs[0,i]) if i == 0 else plt.subplot(gs[0,i], sharey=ax0)
				plt.setp(ax0.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax0.get_yticklabels(), visible=True)
				plt.setp(ax0.get_xticklabels(), visible=False) 
				#plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				for j in range(len(ydata1[i])):
					plt.plot(xdata, ydata1[i][j], color=colors[i][j], label=legend1[i][j], linestyle=linestyle1[j])
					plt.fill_between(xdata, ydata1[i][j]-err1[i][j], ydata1[i][j]+err1[i][j], alpha=0.3, facecolor=colors[i][j])
				plt.ylabel(ylabel1) if i == 0 else plt.ylabel('')
				plt.yscale('log')
				ax1 = plt.subplot(gs[1,i], sharex=ax0) if i == 0 else plt.subplot(gs[1,i], sharex=ax0, sharey=ax1)
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				plt.setp(ax1.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax1.get_xticklabels(), visible=True)
				for j in range(len(ydata2[i])):
					plt.plot(xdata, ydata2[i][j], color=colors[i][j+1], label=legend1[i][j+1], linestyle=linestyle1[j+1])
					plt.fill_between(xdata, ydata2[i][j]-err2[i][j], ydata2[i][j]+err2[i][j], alpha=0.3, facecolor=colors[i][j+1])
				plt.ylabel(ylabel2) if i == 0 else plt.ylabel('')
				h,l=ax0.get_legend_handles_labels()
				ax0.legend(h,l, fontsize=9)
				plt.xscale('log')
				if i == 0:
					ax0.yaxis.set_tick_params(which='minor', label1On=True, label2On=False)
				else:
					ax0.yaxis.set_tick_params(which='minor', label1On=False, label2On=False)
			if not (ylim1[0] == 0 and ylim1[1] == 0):
				ax0.set_ylim(ylim1)
			if not (ylim2[0] == 0 and ylim2[1] == 0):
				ax1.set_ylim(ylim2)
			figure.text(0.5, 0, xlabel, ha='center', fontsize=10)
			#plt.tight_layout()
			plt.subplots_adjust(wspace=0.01, hspace=0.01)
			return figure
		### Thickness as a function of length, LCDM + Symmetron and LCDM + f(R)
		ThickVsLen_Symm = pf.Call_plot_sameX_OLD(Common_bin_length, Mean_thickness[SymmLCDM], Length_label, Mean_Thickness_label, self.Symm_legends,
											color=self.Plot_colors_symm, logscale='loglog', linestyles=self.Linestyles)
		ThickVsLen_fofr = pf.Call_plot_sameX_OLD(Common_bin_length, Mean_thickness[FofrLCDM], Length_label, Mean_Thickness_label, self.fofr_legends,
											color=self.Plot_colors_fofr, logscale='loglog', linestyles=self.Linestyles)
		### Relative difference with errobar, Symmetron and f(R) seperate, base model = LCDM
		ThickVsLen_RelErr_Symm = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Symm_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_thickness[i], alpha=0.7, color=self.Plot_colors_symm[i+1])
			plt.fill_between(Common_bin_length, RelDiff_mean_thickness[i]-Prop_err_mean_thickness[i], RelDiff_mean_thickness[i]+Prop_err_mean_thickness[i],
							 alpha=0.4, facecolor=self.Plot_colors_symm[i+1])
		plt.legend(self.Symm_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('$(R_{T,i} - R_{T,{\Lambda CDM}})/R_{T,{\Lambda CDM}}$')
		plt.xscale('log')

		ThickVsLen_RelErr_fofr = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Fofr_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_thickness[i], alpha=0.7, color=self.Plot_colors_fofr[i-3])
			plt.fill_between(Common_bin_length, RelDiff_mean_thickness[i]-Prop_err_mean_thickness[i], RelDiff_mean_thickness[i]+Prop_err_mean_thickness[i],
							 alpha=0.4, facecolor=self.Plot_colors_fofr[i-3])
		plt.legend(self.fofr_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('$(R_{T,i} - R_{T,{\Lambda CDM}})/R_{T,{\Lambda CDM}}$')
		plt.xscale('log')

		### Thickness vs length, using gridspec call
		ThickVsLen_Symm_gridspec = pf.Do_gridspec_sameX(Common_bin_length, [Mean_thickness[SymmLCDM]], [RelDiff_mean_thickness[Symm_only-1]], 
														Length_label, Mean_Thickness_label, Mean_Thickness_label_reldiff, self.Symm_legends, 
														self.Plot_colors_symm, Primerror=[Mean_thickness_std[SymmLCDM]], 
														Secerror=[Prop_err_mean_thickness[Symm_only-1]], linestyles=self.Linestyles, reldiff=True, 
														fillbetween=True, rowcol=[2,1], xscale='log', xscale_diff='log', yscale='log', 
														legend_anchor=False, nullform=False, figsize=(6,4))
		ThickVsLen_fofr_gridspec = pf.Do_gridspec_sameX(Common_bin_length, [Mean_thickness[FofrLCDM]], [RelDiff_mean_thickness[Fofr_only-1]], 
														Length_label, Mean_Thickness_label, Mean_Thickness_label_reldiff, self.fofr_legends,
														self.Plot_colors_fofr, Primerror=[Mean_thickness_std[FofrLCDM]], 
														Secerror=[Prop_err_mean_thickness[Fofr_only-1]], linestyles=self.Linestyles, reldiff=True, 
														fillbetween=True, rowcol=[2,1], xscale='log', xscale_diff='log', yscale='log', 
														legend_anchor=False, nullform=False, figsize=(6,4))
		ThickVsLen_both_gridspec = Do_bin_speed_gridspec(Common_bin_length, [Mean_thickness[SymmLCDM], Mean_thickness[FofrLCDM]],
														[RelDiff_mean_thickness[Symm_only-1], RelDiff_mean_thickness[Fofr_only-1]],
														[Mean_thickness_std[SymmLCDM], Mean_thickness_std[FofrLCDM]],
														[Prop_err_mean_thickness[Symm_only-1], Prop_err_mean_thickness[Fofr_only-1]],
														[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
														 self.Linestyles, Length_label, Mean_Thickness_label, Mean_Thickness_label_reldiff, ylim2=(-0.25, 0.05))

		### Mass as a function of length, LCDM + Symmetron and LCDM + f(R)
		MassVsLen_Symm = pf.Call_plot_sameX_OLD(Common_bin_length, Mean_mass[SymmLCDM], Length_label, Mean_Mass_label, self.Symm_legends,
											color=self.Plot_colors_symm, logscale='loglog')
		MassVsLen_fofr = pf.Call_plot_sameX_OLD(Common_bin_length, Mean_mass[FofrLCDM], Length_label, Mean_Mass_label, self.fofr_legends,
											color=self.Plot_colors_fofr, logscale='loglog')
		### Relative difference with errorbar, Symmetron and f(R) seperate, base model = LCDM
		MassVsLen_RelErr_Symm = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Symm_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_mass[i], alpha=0.7, color=self.Plot_colors_symm[i+1])
			plt.fill_between(Common_bin_length, RelDiff_mean_mass[i]-Prop_err_mean_mass[i], RelDiff_mean_mass[i]+Prop_err_mean_mass[i],
							 alpha=0.4, facecolor=self.Plot_colors_symm[i+1])
		plt.legend(self.Symm_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('$(M_i - M_{\Lambda CDM})/M_{\Lambda CDM}$')
		plt.xscale('log')
		MassVsLen_RelErr_fofr = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in (Fofr_only-1):
			plt.plot(Common_bin_length, RelDiff_mean_mass[i], alpha=0.7, color=self.Plot_colors_fofr[i-3])
			plt.fill_between(Common_bin_length, RelDiff_mean_mass[i]-Prop_err_mean_mass[i], RelDiff_mean_mass[i]+Prop_err_mean_mass[i],
							 alpha=0.4, facecolor=self.Plot_colors_fofr[i-3])
		plt.legend(self.fofr_legends[1:])
		plt.xlabel(Length_label)
		plt.ylabel('$(M_i - M_{\Lambda CDM})/M_{\Lambda CDM}$')
		plt.xscale('log')
		### Mass vs length, using gridspec call
		MassVsLen_Symm_gridspec = pf.Do_gridspec_sameX(Common_bin_length, [Mean_mass[SymmLCDM]], [RelDiff_mean_mass[Symm_only-1]], 
														Length_label, Mean_Mass_label, Mean_Mass_label_reldiff, self.Symm_legends, self.Plot_colors_symm,
														Primerror=[Mean_mass_std[SymmLCDM]], Secerror=[Prop_err_mean_mass[Symm_only-1]],
														linestyles=self.Linestyles, reldiff=True, fillbetween=True, rowcol=[2,1], xscale='log',
														xscale_diff='log', yscale='log', legend_anchor=False, nullform=False, figsize=(6,4))
		MassVsLen_fofr_gridspec = pf.Do_gridspec_sameX(Common_bin_length, [Mean_mass[FofrLCDM]], [RelDiff_mean_mass[Fofr_only-1]], 
														Length_label, Mean_Mass_label, Mean_Mass_label_reldiff, self.fofr_legends, self.Plot_colors_fofr,
														Primerror=[Mean_mass_std[FofrLCDM]], Secerror=[Prop_err_mean_mass[Fofr_only-1]],
														linestyles=self.Linestyles, reldiff=True, fillbetween=True, rowcol=[2,1], xscale='log',
														xscale_diff='log', yscale='log', legend_anchor=False, nullform=False, ylim_diff=(-0.4, 0.5),
														figsize=(6,4))
		MassVsLen_both_gridspec = Do_bin_speed_gridspec(Common_bin_length, [Mean_mass[SymmLCDM], Mean_mass[FofrLCDM]],
														[RelDiff_mean_mass[Symm_only-1], RelDiff_mean_mass[Fofr_only-1]],
														[Mean_mass_std[SymmLCDM], Mean_mass_std[FofrLCDM]],
														[Prop_err_mean_mass[Symm_only-1], Prop_err_mean_mass[Fofr_only-1]],
														[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
														 self.Linestyles, Length_label, Mean_Mass_label, Mean_Mass_label_reldiff, ylim2=(-0.5,0.5))
		### Mass as a function of thickness, LCDM + Symmetron and LCDM + f(R)
		MassVsThick = plt.figure(figsize=(8,4))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		plt.subplot(1,2,1)
		for i in SymmLCDM:
			plt.loglog(Common_bin_thickness, Mean_mass_vT[i], color=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.subplot(1,2,2)
		for i in range(len(FofrLCDM)):
			j = FofrLCDM[i]
			plt.loglog(Common_bin_thickness, Mean_mass_vT[j], color=self.Plot_colors_fofr[i])
		plt.legend(self.fofr_legends)
		MassVsThick.text(0.5, 0, Thickness_label, ha='center', fontsize=10)
		MassVsThick.text(0, 0.5, Mean_Mass_label, va='center', ha='center', rotation='vertical', fontsize=10)
		#plt.tight_layout()
		### Mass vs thickness, using gridspec call
		MassVsThick_Symm_gridspec = pf.Do_gridspec_sameX(Common_bin_thickness, [Mean_mass_vT[SymmLCDM]], [RelDiff_mean_mass_vT[Symm_only-1]], 
														Thickness_label, Mean_Mass_label, Mean_Mass_label_reldiff, self.Symm_legends, self.Plot_colors_symm,
														Primerror=[Mean_mass_vT_std[SymmLCDM]], Secerror=[Prop_err_mean_mass_vT[Symm_only-1]],
														linestyles=self.Linestyles, reldiff=True, fillbetween=True, rowcol=[2,1], xscale='log',
														xscale_diff='log', yscale='log', legend_anchor=False, figsize=(6,4))
		MassVsThick_fofr_gridspec = pf.Do_gridspec_sameX(Common_bin_thickness, [Mean_mass_vT[FofrLCDM]], [RelDiff_mean_mass_vT[Fofr_only-1]], 
														Thickness_label, Mean_Mass_label, Mean_Mass_label_reldiff, self.fofr_legends, self.Plot_colors_fofr,
														Primerror=[Mean_mass_vT_std[FofrLCDM]], Secerror=[Prop_err_mean_mass_vT[Fofr_only-1]],
														linestyles=self.Linestyles, reldiff=True, fillbetween=True, rowcol=[2,1], xscale='log',
														xscale_diff='log', yscale='log', legend_anchor=False, figsize=(6,4))
		MassVsThick_both_gridspec = Do_bin_speed_gridspec(Common_bin_thickness, [Mean_mass_vT[SymmLCDM], Mean_mass_vT[FofrLCDM]],
														[RelDiff_mean_mass_vT[Symm_only-1], RelDiff_mean_mass_vT[Fofr_only-1]],
														[Mean_mass_vT_std[SymmLCDM], Mean_mass_vT_std[FofrLCDM]],
														[Prop_err_mean_mass_vT[Symm_only-1], Prop_err_mean_mass_vT[Fofr_only-1]],
														[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
														 self.Linestyles, Thickness_label, Mean_Mass_label, Mean_Mass_label_reldiff, ylim2=(-0.1,0.25))

		### Mass vs length, different binning method. Length as x-label
		MassVsLength_Symm_V2 = plt.figure(figsize=(8,4))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		for i in SymmLCDM:
			Average_mass = []
			for j in range(len(Common_bin_length)-1):
				Similar_length = (self.FilLengths[i] >= Common_bin_length[j]) & (self.FilLengths[i] <= Common_bin_length[j+1])
				Masses_included = self.Filament_masses[i][Similar_length]
				Average_mass.append(np.average(Masses_included))
			plt.plot(Common_bin_length[1:], Average_mass, color=self.Plot_colors_symm[i])
		plt.xlabel('$L$ - [Mpc/h]')
		plt.ylabel(r'$\langle M \rangle$ - $[M_\odot /h]$')
		plt.legend(self.Symm_legends)


		############
		############ DO GRIDSPEC PLOTS
		############
		### Number filament mass of the models along with relative differences
		NumMass_symm_logx_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_mass, [Number_mass[SymmLCDM]/1000.0], [RelativeDiff_mass[Symm_only-1]], Mass_label, 
													r'$N \times 1000$', 'Relative difference', self.Symm_legends, self.Plot_colors_symm, 
													Secerror=[Prop_error_mass[Symm_only-1]], xscale='log', linestyles=self.Linestyles, reldiff=True, 
													fillbetween=True, xlim=xlim_mass, ylim_diff=(-0.25, 0.6), rowcol=[2,1], legend_anchor=False, figsize=(6,4),
													ylabel2_font=8)

		NumMass_fofr_logx_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_mass, [Number_mass[FofrLCDM]/1000.0], [RelativeDiff_mass[Fofr_only-1]], Mass_label, 
													r'$N \times 1000$', 'Relative difference', self.fofr_legends, self.Plot_colors_fofr, 
													Secerror=[Prop_error_mass[Fofr_only-1]], xscale='log', linestyles=self.Linestyles, reldiff=True, 
													fillbetween=True, xlim=xlim_mass, ylim_diff=(-0.25, 0.6), rowcol=[2,1], legend_anchor=False, figsize=(6,4),
													ylabel2_font=8)
		### Number filament thickness of the models along with relative differences
		NumThickness_symm_logx_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_thickness, [Number_thickness[SymmLCDM]/1000.0], [RelativeDiff_thickness[Symm_only-1]],
													Thickness_label, r'$N \times 1000$', 'Relative difference', self.Symm_legends, self.Plot_colors_symm,
													Secerror=[Prop_error_thickness[Symm_only-1]], xscale='log', linestyles=self.Linestyles, reldiff=True,
													fillbetween=True, xlim=xlim_thickness, ylim_diff=(-0.6, 0.9), rowcol=[2,1], legend_anchor=False,
													 figsize=(6,4), ylabel2_font=8)
		NumThickness_fofr_logx_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_thickness, [Number_thickness[FofrLCDM]/1000.0], [RelativeDiff_thickness[Fofr_only-1]],
													Thickness_label, r'$N \times 1000$', 'Relative difference', self.fofr_legends, self.Plot_colors_fofr,
													Secerror=[Prop_error_thickness[Fofr_only-1]], xscale='log', linestyles=self.Linestyles, reldiff=True,
													fillbetween=True, xlim=xlim_thickness, ylim_diff=(-0.6, 0.9), rowcol=[2,1], legend_anchor=False,
													figsize=(6,4), ylabel2_font=8)

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
		self.savefigure(RelDiff_mass_fofr_err, 'Relative_difference_mass_error_cFofr')
		######## Thickness histograms
		self.savefigure(NumThickness_all, 'Filament_Thickness_distribution')
		self.savefigure(NumThickness_Symm, 'Filament_Thickness_distribution_cSymmetron')
		self.savefigure(NumThickness_fofr, 'Filament_Thickness_distribution_cFofr')
		self.savefigure(NumThickness_Symm_logX, 'Filament_Thickness_distribution_logX_cSymmetron')
		self.savefigure(NumThickness_fofr_logX, 'Filament_Thickness_distribution_logX_cFofr')
		self.savefigure(Reldiff_num_thick_Symm, 'Relative_difference_thickness_cSymmetron')
		self.savefigure(Reldiff_num_thick_Symm_logX, 'Relative_difference_thickness_logX_cSymmetron')
		self.savefigure(Reldiff_num_thick_fofr, 'Relative_difference_thickness_cFofr')
		self.savefigure(Reldiff_num_thick_fofr_logX, 'Relative_difference_thickness_logX_cFofr')
		######## Density histograms
		self.savefigure(NumDensity_Symm, 'Filament_density_distribution_cSymmetron')
		self.savefigure(NumDensity_fofr, 'Filament_density_distribution_cFofr')
		self.savefigure(NumDensity_Symm_logX, 'Filament_density_distribution_logX_cSymmetron')
		self.savefigure(NumDensity_fofr_logX, 'Filament_density_distribution_logX_cFofr')
		self.savefigure(RelDiff_num_density_Symm, 'Relative_difference_density_cSymmetron')
		self.savefigure(RelDiff_num_density_fofr, 'Relative_difference_density_cFofr')
		self.savefigure(RelDiff_num_density_Symm_logX, 'Relative_difference_density_logX_cSymmetron')
		self.savefigure(RelDiff_num_density_fofr_logX, 'Relative_difference_density_logX_cSymmetron')
		######## Compare different properties
		self.savefigure(ThickVsLen_Symm, 'ThicknessVsLength_cSymmetron')
		self.savefigure(ThickVsLen_fofr, 'ThicknessVsLength_cFofr')
		self.savefigure(ThickVsLen_RelErr_Symm, 'ThicknessVsLength_Reldiff_cSymmetron')
		self.savefigure(ThickVsLen_RelErr_fofr, 'ThicknessVsLength_Reldiff_cFofr')
		self.savefigure(ThickVsLen_Symm_gridspec, 'ThicknessVsLength_cSymmetron_gridspec')
		self.savefigure(ThickVsLen_fofr_gridspec, 'ThicknessVsLength_cFofr_gridspec')
		self.savefigure(ThickVsLen_both_gridspec, 'ThicknessVsLength_both_gridspec')
		self.savefigure(MassVsLen_Symm, 'MassVsLength_cSymmetron')
		self.savefigure(MassVsLen_fofr, 'MassVsLength_cFofr')
		self.savefigure(MassVsLen_RelErr_Symm, 'MassVsLength_Reldiff_cSymmetron')
		self.savefigure(MassVsLen_RelErr_fofr, 'MassVsLength_Reldiff_cFofr')
		self.savefigure(MassVsLen_Symm_gridspec, 'MassVsLength_cSymmetron_gridspec')
		self.savefigure(MassVsLen_fofr_gridspec, 'MassVsLength_cFofr_gridspec')
		self.savefigure(MassVsLen_both_gridspec, 'MassVsLength_both_gridspec')
		self.savefigure(MassVsThick, 'MassVsThickness')
		self.savefigure(MassVsThick_Symm_gridspec, 'MassVsThickness_cSymmetron_gridspec')
		self.savefigure(MassVsThick_fofr_gridspec, 'MassVsThickness_cFofr_gridspec')
		self.savefigure(MassVsThick_both_gridspec, 'MassVsThickness_both_gridspec')
		self.savefigure(MassVsLength_Symm_V2, 'MassVsLength_Symm_V2')
		### Gridspec plots
		self.savefigure(NumMass_symm_logx_GRIDSPEC, 'NumberMass_logX_cSymmetron_gridspec', dpi_mult=2)
		self.savefigure(NumMass_fofr_logx_GRIDSPEC, 'NumberMass_logX_cFofr_gridspec', dpi_mult=2)
		self.savefigure(NumThickness_symm_logx_GRIDSPEC, 'NumberThickness_logX_cSymmetron_gridspec', dpi_mult=2)
		self.savefigure(NumThickness_fofr_logx_GRIDSPEC, 'NumberThickness_logX_cFofr_gridspec', dpi_mult=2)
		plt.close('all')	# Clear all current windows to free memory

	def Velocity_profiles(self, All_speeds, Part_distances, speedtype='Speed'):
		""" 
		Plots data related to the velocity profiles 
		Can either be the total speed, orthogonal speed or parallel speed. 
		The speedtype name must be the same as the input speed array, else shit happens.
		"""
		def get_data(Models_run, filenames, xbins, p_ranges, prop, prop_name, binnums):
			Mean_profiles = []
			Mean_stds = []
			for j in range(len(p_ranges)-1):
				Temp_prof = []
				Temp_std = []
				for ij in range(len(Models_run)):
					i = Models_run[ij]
					Mean_prof, Mean_prof_std = self.Compute_similar_profiles(prop[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
													xbins, p_ranges[j], p_ranges[j+1], prop_name, filenames[ij], newbinning=binnums)
					Temp_prof.append(Mean_prof)
					Temp_std.append(Mean_prof_std)
				Mean_profiles.append(Temp_prof)
				Mean_stds.append(Temp_std)
			return Mean_profiles, Mean_stds

		def get_data_reldiffs(Models_run, p_ranges, reldiffs, Properrs):
			Yplot = []
			Yerror = []
			for j in range(len(Mass_ranges)-1):
				yplot_temp = []
				yerror_temp =[]
				for i in Models_run:
					yplot_temp.append(reldiffs[j][i])
					yerror_temp.append(Properrs[j][i])
				Yplot.append(yplot_temp)
				Yerror.append(yerror_temp)
			return Yplot, Yerror

		def store_reldiff_data(prop_array, prop_name, p_ranges, xbins, binnums):
			RelDiff_AvgSpeed = []
			PropErr_AvgSpeed = []
			for j in range(len(p_ranges)-1):
				Reldif_sim_ranges = []
				Properr_sim_ranges = []
				Mean_profile_LCDM, Mean_profile_std_LCDM = self.Compute_similar_profiles(prop_array[0], All_speeds[0], Part_distances[0],
				 																		 self.Thresholds[0], xbins, p_ranges[j],
				 																		 p_ranges[j+1], prop_name, Symm_filenames[0], newbinning=binnums)
				for i in range(1, NModels):
					Mean_profile, Mean_profile_std = self.Compute_similar_profiles(prop_array[i], All_speeds[i], Part_distances[i], self.Thresholds[i],
																	xbins, p_ranges[j], p_ranges[j+1], prop_name, All_filenames[i], newbinning=binnums)
					RelDiff = OF.relative_deviation_singular(Mean_profile_LCDM, Mean_profile)
					PropErr = OF.Propagate_error_reldiff(Mean_profile_LCDM, Mean_profile, Mean_profile_std_LCDM, Mean_profile_std)
					Reldif_sim_ranges.append(RelDiff)
					Properr_sim_ranges.append(PropErr)
				RelDiff_AvgSpeed.append(np.array(Reldif_sim_ranges))
				PropErr_AvgSpeed.append(np.array(Properr_sim_ranges))
			return RelDiff_AvgSpeed, PropErr_AvgSpeed

		velocity_savefile_dir = 'ModelComparisons/VelocityAnalysis/'
		if self.raw_filetype == 'png':
			velocity_savefile_dir += 'PNG/'
		elif self.raw_filetype == 'pdf':
			velocity_savefile_dir += 'PDF/'
		#self.filetype = '.' + self.raw_filetype
		sigma_name_folder = 'Sigma'+str(self.Nsigma) + '/'
		velocity_savefile_dir += sigma_name_folder

		OK_speedtypes = ['Speed', 'Orthogonal', 'Parallel', 'Density']
		Check = 0
		for spt in OK_speedtypes:
			if speedtype == spt:
				Check += 1
		if not Check:
			raise ValueError("Argument speedtype not set properly! Currently speedtype = ", speedtype, ". Try speedtype=Speed, Orthogonal or Parallel")
		self.Speed_filename = speedtype
		velocity_savefile_dir += speedtype + '/'

		velocity_results_dir = os.path.join(savefile_directory, velocity_savefile_dir)
		if not os.path.isdir(velocity_results_dir):
			os.makedirs(velocity_results_dir)


		print 'Computing for ', speedtype 
		######## Compute relevant stuff ########
		do_fb = True
		binnum = 20
		NModels = len(All_speeds)
		ylimits = (0,0)
		#Mass_label = 'Filament mass - [$M_\odot / h$]'
		#Number_label = '$N$ filaments'
		#Thickness_label = 'Filament thickness - [Mpc/h]'
		#Length_label = 'Filament length - [Mpc/h]'
		#Mean_Mass_label = 'Mean filament mass - [$M_\odot / h$]'
		#Mean_Thickness_label = 'Mean filament thickness - [Mpc/h]'
		Mass_label = '$M$ - [$M_\odot / h$]'
		Number_label = '$N$'
		Thickness_label = '$R_T$ - [Mpc/$h$]'
		Length_label = '$L$ - [Mpc/$h$]'
		Mean_Mass_label = r'$\bar{M} - [\mathrm{Mpc}/h]$'
		Mean_Thickness_label = r'$\bar{R_T} - [\mathrm{Mpc}/h]$'
		Distance_normalized_label = '$r/R_T$'
		Special_y_scale_density = 'no'
		if speedtype == 'Speed':
			Average_speed_label =  r'$\langle v \rangle - [\mathrm{km}/\mathrm{s}]$'
			Average_speed_label_nounit = r'$\langle v \rangle$'
			#Reldiff_label_avgspeed = r'$(\langle v_i \rangle - \langle v_{\Lambda \mathrm{CDM}} \rangle)/\langle v_{\Lambda \mathrm{CDM}} \rangle$'
			Reldiff_label_avgspeed = r'$\langle v_i \rangle/\langle v_{\Lambda \mathrm{CDM}} \rangle - 1$'
			ylim_mass = (400, 1400)
			ylim_mass_rld = (-0.05, 0.25)
			ylim_len = (500, 800)
			ylim_len_rld = (-0.05, 0.3)
			ylim_thick = (400, 1400)
			ylim_thick_rld = (-0.05, 0.3)
		elif speedtype == 'Orthogonal':
			#Average_speed_label =  r'$\langle v_\perp \rangle - [\mathrm{km}/\mathrm{s}]$'
			#Average_speed_label_nounit = r'$\langle v_\perp \rangle$'
			#Reldiff_label_avgspeed = r'$\langle v_{\perp,i}  \rangle/\langle v_{\perp, \Lambda \mathrm{CDM}} \rangle - 1$'
			#Reldiff_label_avgspeed = r'$(\langle v_{\perp,i} \rangle - \langle v_{\perp,\Lambda \mathrm{CDM}} \rangle)/\langle v_{\perp, \Lambda \mathrm{CDM}} \rangle$'
			Average_speed_label =  r'$\langle v_r \rangle - [\mathrm{km}/\mathrm{s}]$'
			Average_speed_label_nounit = r'$\langle v_r \rangle$'
			Reldiff_label_avgspeed = r'$\langle v_{r,i}  \rangle/\langle v_{r, \Lambda \mathrm{CDM}} \rangle - 1$'
			ylim_mass = (-1200, -300)
			ylim_mass_rld = (-0.05, 0.25)
			ylim_len = (-600, -350)
			ylim_len_rld = (-0.05, 0.25)
			ylim_thick = (-1200, -300)
			ylim_thick_rld = (-0.05, 0.3)
		elif speedtype == 'Parallel':
			#Average_speed_label =  r'$\langle v_\parallel \rangle - [\mathrm{km}/s]$'
			#Average_speed_label_nounit = r'$\langle v_{\parallel} \rangle$'
			#Reldiff_label_avgspeed = r'$\langle v_{\parallel,i}  \rangle/\langle v_{\parallel, \Lambda \mathrm{CDM}} \rangle - 1$'
			#Reldiff_label_avgspeed = r'$(\langle v_{\parallel,i} \rangle - \langle v_{\parallel,\Lambda \mathrm{CDM}} \rangle)/\langle v_{\parallel,\Lambda \mathrm{CDM}} \rangle$'
			Average_speed_label =  r'$\langle v_t \rangle - [\mathrm{km}/s]$'
			Average_speed_label_nounit = r'$\langle v_{t} \rangle$'
			Reldiff_label_avgspeed = r'$\langle v_{t,i}  \rangle/\langle v_{t, \Lambda \mathrm{CDM}} \rangle - 1$'
			ylimits = (-0.5,0.5)
			ylim_mass = (0, 0)
			ylim_mass_rld = (-1, 1)
			ylim_len = (0, 0)
			ylim_len_rld = (-1, 1)
			ylim_thick = (0, 0)
			ylim_thick_rld = (-1, 1)
		elif speedtype == 'Density':
			Average_speed_label = r'$\langle \rho \rangle - [M_\odot h^2/\mathrm{Mpc}^3] \times 1e14$'
			Average_speed_label_nounit = r'$\langle \rho \rangle$'
			#Reldiff_label_avgspeed = r'$(\langle \rho_i \rangle - \langle \rho_{\Lambda \mathrm{CDM}} \rangle)/(\langle \rho_{\Lambda \mathrm{CDM}} \rangle)$'
			Reldiff_label_avgspeed = r'$\langle \rho_i \rangle/\langle \rho_{\Lambda \mathrm{CDM}} \rangle - 1$'
			#do_fb = True
			ylimits = (-0.1,0.4)
			Special_y_scale_density = 'log'
			self.Dist_mass_filter = 0
			ylim_mass = (0, 0)
			ylim_mass_rld = (-0.1, 0.5)
			ylim_len = (0, 0)
			ylim_len_rld = (-0.1, 0.5)
			ylim_thick = (0, 0)
			ylim_thick_rld = (-0.2, 0.55)
		
		### Further filtering based on masses
		for i in range(NModels):
			Mass_filter = self.Mass_filters[i]
			All_speeds[i] = All_speeds[i][Mass_filter]
			if self.Dist_mass_filter == 0:
				Part_distances[i] = Part_distances[i][Mass_filter]
		self.Dist_mass_filter = 1
		#SymmLCDM = np.array([0,1,2,3,4])
		#FofrLCDM = np.array([0,5,6,7])
		#Symm_only = np.array([1,2,3,4])
		#Fofr_only = np.array([5,6,7])
		SymmLCDM = self.SymmLCDM_ids
		FofrLCDM = self.FofrLCDM_ids
		Symm_only = self.SymmLCDM_ids[1:]
		Fofr_only = self.FofrLCDM_ids[1:]
		All_speeds = np.asarray(All_speeds)
		Symm_filenames = ['LCDM', 'Symm_A', 'Symm_B', 'Symm_C', 'Symm_D']
		Fofr_filenames = ['LCDM', 'fofr4', 'fofr5', 'fofr6']
		All_filenames = ['LCDM', 'Symm_A', 'Symm_B', 'Symm_C', 'Symm_D', 'fofr4', 'fofr5', 'fofr6']
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
		Common_bin_length = OF.Get_common_bin_logX(self.FilLengths, binnum=binnum)

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
		AverageSpeed_AllFils = plt.figure(figsize=(12,5))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax1 = plt.subplot(2,1,1)
		for i in SymmLCDM:
			plt.plot(Common_bin_distances_normalized, Mean_speed_allFils[i], color=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		if speedtype == 'Density':
			plt.yscale('log')
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax2 = plt.subplot(2,1,2, sharex=ax1)
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Common_bin_distances_normalized, Mean_speed_allFils[i], color=self.Plot_colors_fofr[ij])
		plt.legend(self.fofr_legends)
		AverageSpeed_AllFils.text(0.5, 0, Distance_normalized_label, ha='center', fontsize=10)
		AverageSpeed_AllFils.text(0, 0.5, Average_speed_label, va='center', ha='center', rotation='vertical', fontsize=10)
		plt.xscale('log')
		if speedtype == 'Density':
			plt.yscale('log')
		plt.subplots_adjust(wspace=0.01, hspace=0.04)
		#plt.tight_layout()
		####################
		#################### SIMILAR MASSES
		####################
		### Average speed of filaments with similar masses, comparing LCDM + Symmetron
		Mass_titles = ['$M \in [10^{12}, 10^{13}]M_\odot/h$', '$M \in [10^{13}, 10^{14}]M_\odot/h$', '$M \in [10^{14}, 10^{15}]M_\odot/h$']
		Reldiff_speed_xlimits = (Common_bin_distances_normalized[1], Common_bin_distances_normalized[-1])
		Mass_ranges = [1e12, 1e13, 1e14, 1e15]   # Units of M_sun/h, maybe use min and max of mass bin?
		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Mass_ranges, self.Filament_masses, 'Mass', binnum)
		AverageSpeed_SimilarMass_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, title=Mass_titles,
														linestyles=self.Linestyles, xlim=(0,1), yscale=Special_y_scale_density)
		AverageSpeed_SimilarMass_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label,
														Average_speed_label, self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, 
														title=Mass_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)
		### Average speed of filaments with similar masses, comparing LCDM + f(R)
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Mass_ranges, self.Filament_masses, 'Mass', binnum)
		AverageSpeed_SimilarMass_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb, title=Mass_titles,
														linestyles=self.Linestyles, xlim=(0,1), yscale=Special_y_scale_density)
		AverageSpeed_SimilarMass_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, 
														Average_speed_label, self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb,
														title=Mass_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)
		### Relative difference of the average speed, for similar masses.
		RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass = store_reldiff_data(self.Filament_masses, 'Mass', Mass_ranges, 
																						Common_bin_distances_normalized, binnum)
		### Plotting relative differences with error
		### Symmetron differences to LCDM
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass)
		RelDiff_SimilarMass_plot_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
														self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb, ylim=ylimits, 
														title=Mass_titles, linestyles=self.Linestyles, reldiff=True, xlim=(0,1))
		RelDiff_SimilarMass_plot_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
														self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb, ylim=ylimits, 
														title=Mass_titles, linestyles=self.Linestyles, reldiff=True, xscale='log')
		### f(R) difference to LCDM
		Yplot, Yerror = get_data_reldiffs(Fofr_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass)
		RelDiff_SimilarMass_plot_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
														self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, ylim=ylimits, 
														title=Mass_titles, linestyles=self.Linestyles, reldiff=True, xlim=(0,1))
		RelDiff_SimilarMass_plot_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
														self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, ylim=ylimits, 
														title=Mass_titles, linestyles=self.Linestyles, reldiff=True, xscale='log')
		####################
		#################### SIMILAR LENGTHS
		####################
		### Average speed of filament of similar length. Symmetron + LCDM comparison
		Length_titles = ['$L \in [1,5]$ Mpc/h', '$L \in [5,10]$ Mpc/h', '$L \in [10,20]$ Mpc/h']
		Length_ranges = [1, 5, 10, 20]   # Units of Mpc/h, maybe use min and max of mass bin?

		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Length_ranges, self.FilLengths, 'Length', binnum)
		AverageSpeed_SimilarLength_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, title=Length_titles,
														linestyles=self.Linestyles, yscale=Special_y_scale_density)
		AverageSpeed_SimilarLength_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, 
														Average_speed_label, self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, 
														title=Length_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)
		#### Average speed of filament of similar length. f(R) + LCDM comparison
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Length_ranges, self.FilLengths, 'Length', binnum)
		AverageSpeed_SimilarLength_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb, title=Length_titles,
														linestyles=self.Linestyles, yscale=Special_y_scale_density)
		AverageSpeed_SimilarLength_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, 
														Average_speed_label, self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb, 
														title=Length_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)
		### Relative differences of filaments of similar length
		RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength = store_reldiff_data(self.FilLengths, 'Length', Length_ranges, 
																						Common_bin_distances_normalized, binnum)

		### Plotting relative differences with error
		### Symmetron differences to LCDM
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Length_ranges, RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength)
		RelDiff_SimilarLen_plot_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb, 
													ylim=ylimits, title=Length_titles, linestyles=self.Linestyles, reldiff=True, xlim=(0,1))
		RelDiff_SimilarLen_plot_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb,  ylim=ylimits, 
													title=Length_titles, linestyles=self.Linestyles, reldiff=True, xscale='log')
		### f(R) difference to LCDM
		Yplot, Yerror = get_data_reldiffs(Fofr_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength)
		RelDiff_SimilarLen_plot_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, 
													ylim=ylimits, title=Length_titles, linestyles=self.Linestyles, reldiff=True, xlim=(0,1))
		RelDiff_SimilarLen_plot_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, ylim=ylimits, 
													title=Length_titles, linestyles=self.Linestyles, reldiff=True, xscale='log')
		####################
		#################### SIMILAR THICKNESS
		####################
		### Average speed of filament of similar length. Symmetron + LCDM comparison
		Thickness_titles = ['$R_T \in [0.1,1]$ Mpc/h', '$R_T \in [1,5]$ Mpc/h', '$R_T \in [5,10]$ Mpc/h']
		Thickness_ranges = [0.1, 1, 5, 10]   # Units of Mpc/h, maybe use min and max of mass bin?
		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Thickness_ranges, self.Thresholds, 'Thickness', binnum)
		AverageSpeed_SimilarThickness_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, title=Thickness_titles,
														linestyles=self.Linestyles, yscale=Special_y_scale_density)
		AverageSpeed_SimilarThickness_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances, Mean_profiles, Distance_normalized_label, 
														Average_speed_label, self.Symm_legends, self.Plot_colors_symm, error=Mean_stds, fillbetween=do_fb, 
														title=Thickness_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)

		#### Average speed of filament of similar length. f(R) + LCDM comparison
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Thickness_ranges, self.Thresholds, 'Thickness', binnum)
		AverageSpeed_SimilarThickness_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Mean_profiles, Distance_normalized_label, Average_speed_label,
														self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb, title=Thickness_titles,
														linestyles=self.Linestyles, yscale=Special_y_scale_density)
		AverageSpeed_SimilarThickness_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances, Mean_profiles, Distance_normalized_label, 
														Average_speed_label, self.fofr_legends, self.Plot_colors_fofr, error=Mean_stds, fillbetween=do_fb, 
														title=Thickness_titles, linestyles=self.Linestyles, xscale='log', yscale=Special_y_scale_density)
		### Relative differences of filaments of similar thickness
		RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness = store_reldiff_data(self.Thresholds, 'Thickness', Thickness_ranges, 
																						Common_bin_distances_normalized, binnum)
		### Plotting relative differences with error
		### Symmetron differences to LCDM
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Thickness_ranges, RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness)
		RelDiff_SimilarThickness_plot_Symm = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb, title=Thickness_titles,
													linestyles=self.Linestyles, reldiff=True, xlim=(0,1), ylim=ylimits)
		RelDiff_SimilarThickness_plot_Symm_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, 
													Reldiff_label_avgspeed, self.Symm_legends[1:], self.Plot_colors_symm[1:], error=Yerror, fillbetween=do_fb,
													title=Thickness_titles, linestyles=self.Linestyles, reldiff=True, xscale='log', 
													xlim=(0,1), ylim=ylimits)
		### f(R) difference to LCDM
		Yplot, Yerror = get_data_reldiffs(Fofr_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness)
		RelDiff_SimilarThickness_plot_fofr = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, Reldiff_label_avgspeed,
													self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, title=Thickness_titles,
													linestyles=self.Linestyles, reldiff=True, xlim=(0,1), ylim=ylimits)

		RelDiff_SimilarThickness_plot_fofr_logx = pf.Do_subplots_sameX(Common_bin_distances_normalized, Yplot, Distance_normalized_label, 
													Reldiff_label_avgspeed, self.fofr_legends[1:], self.Plot_colors_fofr[1:], error=Yerror, fillbetween=do_fb, 
													title=Thickness_titles, linestyles=self.Linestyles, reldiff=True, xscale='log',
													xlim=(0,1), ylim=ylimits)

		####################
		#################### OVER DIFFERENT MASS BINS
		####################
		def Do_bin_speed_gridspec(xdata, ydata1, ydata2, err1, err2, colors, legend1, linestyle1, xlabel, ylabel1, ylabel2, rowcol=[2,2], ylim1=(0,0), ylim2=(0,0)):
			figure = plt.figure()
			gs = gridspec.GridSpec(rowcol[0],rowcol[1])
			for i in range(len(ydata1)):
				ax0 = plt.subplot(gs[0,i]) if i == 0 else plt.subplot(gs[0,i], sharey=ax0)
				plt.setp(ax0.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax0.get_yticklabels(), visible=True)
				plt.setp(ax0.get_xticklabels(), visible=False) 
				#plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				for j in range(len(ydata1[i])):
					plt.plot(xdata, ydata1[i][j], color=colors[i][j], label=legend1[i][j], linestyle=linestyle1[j])
					plt.fill_between(xdata, ydata1[i][j]-err1[i][j], ydata1[i][j]+err1[i][j], alpha=0.3, facecolor=colors[i][j])
				plt.ylabel(ylabel1) if i == 0 else plt.ylabel('')
				ax1 = plt.subplot(gs[1,i], sharex=ax0) if i == 0 else plt.subplot(gs[1,i], sharex=ax0, sharey=ax1)
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				plt.setp(ax1.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax1.get_xticklabels(), visible=True)
				for j in range(len(ydata2[i])):
					plt.plot(xdata, ydata2[i][j], color=colors[i][j+1], label=legend1[i][j+1], linestyle=linestyle1[j+1])
					plt.fill_between(xdata, ydata2[i][j]-err2[i][j], ydata2[i][j]+err2[i][j], alpha=0.3, facecolor=colors[i][j+1])
				plt.ylabel(ylabel2) if i == 0 else plt.ylabel('')
				h,l=ax0.get_legend_handles_labels()
				ax0.legend(h,l, fontsize=9)
				plt.xscale('log')
			if not (ylim1[0] == 0 and ylim1[1] == 0):
				ax0.set_ylim(ylim1)
			if not (ylim2[0] == 0 and ylim2[1] == 0):
				ax1.set_ylim(ylim2)
			figure.text(0.5, 0, xlabel, ha='center', fontsize=10)
			#figure.text(0, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=10)
			#ax0.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True)
			#ax1.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True)
			plt.tight_layout()
			return figure
		#### Average speed for different mass bins
		xlim_mass = (Common_bin_mass[0], Common_bin_mass[-1])
		Compare_both = [SymmLCDM, FofrLCDM]
		Compare_both_legends = [self.Symm_legends, self.fofr_legends]
		Compare_both_colors = [self.Plot_colors_symm, self.Plot_colors_fofr]
		Average_speed_massbin = []
		Average_speed_massbin_std = []
		Mid_common_bin_mass = np.array([np.linspace(Common_bin_mass[i], Common_bin_mass[i+1],3)[1] for i in range(len(Common_bin_mass)-1)])
		AverageSpeed_MassBins = plt.figure(figsize=(8,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(1,2,1)
		for j in range(len(Compare_both)):
			ij = 0
			if j > 0:
				ax = plt.subplot(1,2,j+1, sharey=ax)
				plt.setp(ax.get_yticklabels(), visible=False)
			for k in Compare_both[j]:
				Average_speed, Error_speed = self.Compute_average_speeds_Propertybinned(All_speeds[k], self.Filament_masses[k], Common_bin_mass)
				if not ((j == 1) and (k == 0)):
					Average_speed_massbin.append(Average_speed)
					Average_speed_massbin_std.append(Error_speed)
				#plt.errorbar(Common_bin_mass[1:], Average_speed, Error_speed)
				plt.plot(Mid_common_bin_mass, Average_speed, color=Compare_both_colors[j][ij], linestyle=self.Linestyles[ij])
				plt.fill_between(Mid_common_bin_mass, Average_speed - Error_speed, Average_speed + Error_speed, alpha=0.3, facecolor=Compare_both_colors[j][ij])
				plt.xscale('log')
				ij += 1
			plt.legend(Compare_both_legends[j])
		plt.xlim(xlim_mass)
		AverageSpeed_MassBins.text(0.5, 0, Mass_label, ha='center', fontsize=10)
		AverageSpeed_MassBins.text(0, 0.5, Average_speed_label, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		Average_speed_massbin = np.array(Average_speed_massbin)
		Average_speed_massbin_std = np.array(Average_speed_massbin_std)
		####
		#### Relative difference of the average speeds, using mass bins
		####
		#RelDiffs_AvgSpeed_massbins = [OF.relative_deviation(Average_speed_massbin, i) for i in range(1, NModels)]
		RelDiffs_AvgSpeed_massbins = np.array([OF.relative_deviation_singular(Average_speed_massbin[0], Average_speed_massbin[i]) for i in range(1, NModels)])
		PropErr_AvgSpeed_massbins = np.array([OF.Propagate_error_reldiff(Average_speed_massbin[0], Average_speed_massbin[i],
															Average_speed_massbin_std[0], Average_speed_massbin_std[i]) for i in range(1, NModels)])
		AverageSpeed_RelativeDifference_MassBins = plt.figure(figsize=(12,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		set_y_limit = 0
		ax = plt.subplot(1,2,1)
		plt.plot(Mid_common_bin_mass, np.zeros(len(Mid_common_bin_mass)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Symm_only:
			if np.nanmax(RelDiffs_AvgSpeed_massbins[i-1]+PropErr_AvgSpeed_massbins[i-1]) > 1:
				set_y_limit = 1
			plt.plot(Mid_common_bin_mass, RelDiffs_AvgSpeed_massbins[i-1], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i-1])
			plt.fill_between(Mid_common_bin_mass, RelDiffs_AvgSpeed_massbins[i-1]-PropErr_AvgSpeed_massbins[i-1],
							 RelDiffs_AvgSpeed_massbins[i-1]+PropErr_AvgSpeed_massbins[i-1], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.plot(Mid_common_bin_mass, np.zeros(len(Mid_common_bin_mass)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for ij in range(len(Fofr_only)):
			i = Fofr_only[ij]
			plt.plot(Mid_common_bin_mass, RelDiffs_AvgSpeed_massbins[i-1], color=self.Plot_colors_fofr[ij+1], linestyle=self.Linestyles[ij])
			plt.fill_between(Mid_common_bin_mass, RelDiffs_AvgSpeed_massbins[i-1]-PropErr_AvgSpeed_massbins[i-1],
							 RelDiffs_AvgSpeed_massbins[i-1]+PropErr_AvgSpeed_massbins[i-1], alpha=0.4, facecolor=self.Plot_colors_fofr[ij+1])
		plt.legend(self.fofr_legends)
		plt.xscale('log')
		plt.xlim(xlim_mass)
		#if set_y_limit:
		plt.ylim(-0.1,0.3)
		AverageSpeed_RelativeDifference_MassBins.text(0.5, 0, Mass_label, ha='center', fontsize=10)
		AverageSpeed_RelativeDifference_MassBins.text(0, 0.5, Reldiff_label_avgspeed, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		# ==============
		Average_speed_massbins_both = Do_bin_speed_gridspec(Mid_common_bin_mass, [Average_speed_massbin[SymmLCDM], Average_speed_massbin[FofrLCDM]],
											[RelDiffs_AvgSpeed_massbins[Symm_only-1], RelDiffs_AvgSpeed_massbins[Fofr_only-1]], 
											[Average_speed_massbin_std[SymmLCDM], Average_speed_massbin_std[FofrLCDM]],
											[PropErr_AvgSpeed_massbins[Symm_only-1], PropErr_AvgSpeed_massbins[Fofr_only-1]],
											[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
											self.Linestyles, Mass_label, Average_speed_label, Reldiff_label_avgspeed, ylim1=ylim_mass, ylim2=ylim_mass_rld)
		####
		#### Average speed for different length bins, with LCDM and Symmetron models
		####
		Average_speed_lengthbins = []
		Average_speed_lengthbins_std = []
		Mid_common_bin_length = np.array([np.linspace(Common_bin_length[i], Common_bin_length[i+1],3)[1] for i in range(len(Common_bin_length)-1)])
		AverageSpeed_LengthBins = plt.figure(figsize=(8,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(1,2,1)
		for j in range(len(Compare_both)):
			ij = 0
			if j > 0:
				ax = plt.subplot(1,2,j+1, sharey=ax)
				plt.setp(ax.get_yticklabels(), visible=False)
			for k in Compare_both[j]:
				Average_speed, Error_speed = self.Compute_average_speeds_Propertybinned(All_speeds[k], self.FilLengths[k], Common_bin_length)  # Maybe a new bin with 21 binning?
				if not ((j == 1) and (k == 0)):
					Average_speed_lengthbins.append(Average_speed)
					Average_speed_lengthbins_std.append(Error_speed)
				plt.plot(Mid_common_bin_length, Average_speed, color=Compare_both_colors[j][ij], linestyle=self.Linestyles[ij])
				plt.fill_between(Mid_common_bin_length, Average_speed - Error_speed, Average_speed + Error_speed, alpha=0.3, facecolor=Compare_both_colors[j][ij])
				plt.xscale('log')
				ij += 1
			plt.legend(Compare_both_legends[j])
		AverageSpeed_LengthBins.text(0.5, 0, Mass_label, ha='center', fontsize=10)
		AverageSpeed_LengthBins.text(0, 0.5, Average_speed_label, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		Average_speed_lengthbins = np.array(Average_speed_lengthbins)
		Average_speed_lengthbins_std = np.array(Average_speed_lengthbins_std)
		###
		### Relative difference of different length bins
		###
		RelDiffs_AvgSpeed_lengthbins = np.array([OF.relative_deviation(Average_speed_lengthbins, i) for i in range(1, NModels)])
		PropErr_AvgSpeed_lengthbins = np.array([OF.Propagate_error_reldiff(Average_speed_lengthbins[0], Average_speed_lengthbins[i],
															Average_speed_lengthbins_std[0], Average_speed_lengthbins_std[i]) for i in range(1, NModels)])
		AverageSpeed_RelativeDifference_LengthBins = plt.figure(figsize=(12,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		set_y_limit = 0
		ax = plt.subplot(1,2,1)
		plt.plot(Mid_common_bin_length, np.zeros(len(Mid_common_bin_length)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Symm_only:
			if np.nanmax(RelDiffs_AvgSpeed_lengthbins[i-1]+PropErr_AvgSpeed_lengthbins[i-1]) > 1:
				set_y_limit = 1
			plt.plot(Mid_common_bin_length, RelDiffs_AvgSpeed_lengthbins[i-1], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i-1])
			plt.fill_between(Mid_common_bin_length, RelDiffs_AvgSpeed_lengthbins[i-1]-PropErr_AvgSpeed_lengthbins[i-1],
							 RelDiffs_AvgSpeed_lengthbins[i-1]+PropErr_AvgSpeed_lengthbins[i-1], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.plot(Mid_common_bin_length, np.zeros(len(Mid_common_bin_length)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for ij in range(len(Fofr_only)):
			i = Fofr_only[ij]
			plt.plot(Mid_common_bin_length, RelDiffs_AvgSpeed_lengthbins[i-1], color=self.Plot_colors_fofr[ij+1], linestyle=self.Linestyles[ij])
			plt.fill_between(Mid_common_bin_length, RelDiffs_AvgSpeed_lengthbins[i-1]-PropErr_AvgSpeed_lengthbins[i-1],
							 RelDiffs_AvgSpeed_lengthbins[i-1]+PropErr_AvgSpeed_lengthbins[i-1], alpha=0.4, facecolor=self.Plot_colors_fofr[ij+1])
		plt.legend(self.fofr_legends)
		plt.xscale('log')
		if set_y_limit:
			plt.ylim(-1,1)
		AverageSpeed_RelativeDifference_LengthBins.text(0.5, 0, Mass_label, ha='center', fontsize=10)
		AverageSpeed_RelativeDifference_LengthBins.text(0, 0.5, Reldiff_label_avgspeed, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		# =====
		Average_speed_lengthbins_both = Do_bin_speed_gridspec(Mid_common_bin_length, [Average_speed_lengthbins[SymmLCDM], Average_speed_lengthbins[FofrLCDM]],
											[RelDiffs_AvgSpeed_lengthbins[Symm_only-1], RelDiffs_AvgSpeed_lengthbins[Fofr_only-1]], 
											[Average_speed_lengthbins_std[SymmLCDM], Average_speed_lengthbins_std[FofrLCDM]],
											[PropErr_AvgSpeed_lengthbins[Symm_only-1], PropErr_AvgSpeed_lengthbins[Fofr_only-1]],
											[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
											self.Linestyles, Length_label, Average_speed_label, Reldiff_label_avgspeed, ylim1=ylim_len, ylim2=ylim_len_rld)		
		###
		### Average speeds over different thickness bins
		###
		Average_speed_thicknessbin = []
		Average_speed_thicknessbin_std = []
		Mid_common_bin_thickness = np.array([np.linspace(Common_bin_thickness[i], Common_bin_thickness[i+1],3)[1] for i in range(len(Common_bin_thickness)-1)])
		AverageSpeed_ThicknessBins = plt.figure(figsize=(8,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(1,2,1)
		for j in range(len(Compare_both)):
			ij = 0
			if j > 0:
				ax = plt.subplot(1,2,j+1, sharey=ax)
				plt.setp(ax.get_yticklabels(), visible=False)
			for k in Compare_both[j]:
				Average_speed, Error_speed = self.Compute_average_speeds_Propertybinned(All_speeds[k], self.Thresholds[k], Common_bin_thickness)
				if not ((j == 1) and (k == 0)):
					Average_speed_thicknessbin.append(Average_speed)
					Average_speed_thicknessbin_std.append(Error_speed)
				plt.plot(Mid_common_bin_thickness, Average_speed, '-', color=Compare_both_colors[j][ij], linestyle=self.Linestyles[ij])
				plt.fill_between(Mid_common_bin_thickness, Average_speed-Error_speed, Average_speed+Error_speed, alpha=0.3, facecolor=Compare_both_colors[j][ij])
				plt.xscale('log')
				ij += 1
			plt.legend(Compare_both_legends[j])
		AverageSpeed_ThicknessBins.text(0.5, 0, Thickness_label, ha='center', fontsize=10)
		AverageSpeed_ThicknessBins.text(0, 0.5, Average_speed_label, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		Average_speed_thicknessbin = np.array(Average_speed_thicknessbin)
		Average_speed_thicknessbin_std = np.array(Average_speed_thicknessbin_std)
		#### Relative difference of the average speeds, using mass bins
		RelDiffs_AvgSpeed_thicknessbins = np.array([OF.relative_deviation(Average_speed_thicknessbin, i) for i in range(1, NModels)])
		PropErr_AvgSpeed_thicknessbins = np.array([OF.Propagate_error_reldiff(Average_speed_thicknessbin[0], Average_speed_thicknessbin[i],
															Average_speed_thicknessbin_std[0], Average_speed_thicknessbin_std[i]) for i in range(1, NModels)])
		AverageSpeed_RelativeDifference_ThicknessBins = plt.figure(figsize=(12,6))
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		set_y_limit = 0
		ax = plt.subplot(1,2,1)
		plt.plot(Mid_common_bin_thickness, np.zeros(len(Mid_common_bin_thickness)), 'k-', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Symm_only:
			if np.nanmax(RelDiffs_AvgSpeed_thicknessbins[i-1]+PropErr_AvgSpeed_thicknessbins[i-1]) > 1:
				set_y_limit = 1
			plt.plot(Mid_common_bin_thickness, RelDiffs_AvgSpeed_thicknessbins[i-1], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i-1])
			plt.fill_between(Mid_common_bin_thickness, RelDiffs_AvgSpeed_thicknessbins[i-1]-PropErr_AvgSpeed_thicknessbins[i-1],
							 RelDiffs_AvgSpeed_thicknessbins[i-1]+PropErr_AvgSpeed_thicknessbins[i-1], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.plot(Mid_common_bin_thickness, np.zeros(len(Mid_common_bin_thickness)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for ij in range(len(Fofr_only)):
			i = Fofr_only[ij]
			plt.plot(Mid_common_bin_thickness, RelDiffs_AvgSpeed_thicknessbins[i-1], color=self.Plot_colors_fofr[ij+1], linestyle=self.Linestyles[ij])
			plt.fill_between(Mid_common_bin_thickness, RelDiffs_AvgSpeed_thicknessbins[i-1]-PropErr_AvgSpeed_thicknessbins[i-1],
							 RelDiffs_AvgSpeed_thicknessbins[i-1]+PropErr_AvgSpeed_thicknessbins[i-1], alpha=0.4, facecolor=self.Plot_colors_fofr[ij+1])
		plt.legend(self.fofr_legends)
		plt.xscale('log')
		if set_y_limit:
			plt.ylim(-2,2)
		AverageSpeed_RelativeDifference_ThicknessBins.text(0.5, 0, Thickness_label, ha='center', fontsize=10)
		AverageSpeed_RelativeDifference_ThicknessBins.text(0, 0.5, Reldiff_label_avgspeed, ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()
		# =====
		Average_speed_thicknessbins_both = Do_bin_speed_gridspec(Mid_common_bin_thickness, 
											[Average_speed_thicknessbin[SymmLCDM], Average_speed_thicknessbin[FofrLCDM]],
											[RelDiffs_AvgSpeed_thicknessbins[Symm_only-1], RelDiffs_AvgSpeed_thicknessbins[Fofr_only-1]], 
											[Average_speed_thicknessbin_std[SymmLCDM], Average_speed_thicknessbin_std[FofrLCDM]],
											[PropErr_AvgSpeed_thicknessbins[Symm_only-1], PropErr_AvgSpeed_thicknessbins[Fofr_only-1]],
											[self.Plot_colors_symm, self.Plot_colors_fofr], [self.Symm_legends, self.fofr_legends],
											self.Linestyles, Thickness_label, Average_speed_label, Reldiff_label_avgspeed, ylim1=ylim_thick, ylim2=ylim_thick_rld)
		#######
		####### Gridspec plots
		#######
		def Do_density_reldiff_gridspec(xdata, ydata1, ydata2, err1, err2, color1, color2, legend1, legend2, linestyle1, ylim1, ylim2, titlestop, rowcol=[2,3]):
			figure = plt.figure(figsize=(8,6))
			gs = gridspec.GridSpec(rowcol[0],rowcol[1])
			for i in range(len(titlestop)):
				ax0 = plt.subplot(gs[0,i]) if i == 0 else plt.subplot(gs[0,i], sharey=ax0)
				plt.setp(ax0.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax0.get_yticklabels(), visible=True)
				plt.setp(ax0.get_xticklabels(), visible=False) 
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				for j in range(len(ydata1[i])):
					plt.plot(xdata, ydata1[i][j], color=color1[j], label=legend1[j], linestyle=linestyle1[j])
					plt.fill_between(xdata, ydata1[i][j]-err1[i][j], ydata1[i][j]+err1[i][j], alpha=0.3, facecolor=color1[j])
				plt.title(titlestop[i])
				ax1 = plt.subplot(gs[1,i], sharex=ax0) if i == 0 else plt.subplot(gs[1,i], sharex=ax0, sharey=ax1)
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3, facecolor='k')
				plt.setp(ax1.get_yticklabels(), visible=False) if i > 0 else plt.setp(ax1.get_xticklabels(), visible=True)
				for j in range(len(ydata2[i])):
					plt.plot(xdata, ydata2[i][j], color=color2[j], label=legend2[j], linestyle=linestyle1[j])
					plt.fill_between(xdata, ydata2[i][j]-err2[i][j], ydata2[i][j]+err2[i][j], alpha=0.3, facecolor=color2[j])
			plt.xscale('log')
			ax0.set_ylim(ylim1)
			ax1.set_ylim(ylim2)
			figure.text(0.5, 0, Distance_normalized_label, ha='center', fontsize=10)
			figure.text(0, 0.5, Reldiff_label_avgspeed, ha='center', va='center', rotation='vertical', fontsize=10)
			#ax0.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True)
			#ax1.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True)
			h0,l0=ax0.get_legend_handles_labels()
			ax0.legend(h0,l0, fontsize=9)
			h1,l1=ax1.get_legend_handles_labels()
			ax1.legend(h1,l1, fontsize=9)
			plt.tight_layout()
			return figure
		GS_xscale = 'log'
		GS_yscale = 'linear'
		if speedtype == 'Orthogonal':
			GS_yscale = 'linear'
		if speedtype == 'Parallel':
			special_limit = (-90, 50)
		else:
			special_limit = (0,0)
		## Average speed of filaments with similar mass
		RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass = store_reldiff_data(self.Filament_masses, 'Mass', Mass_ranges, 
																						Common_bin_distances_normalized, binnum)
		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Mass_ranges, self.Filament_masses, 'Mass', binnum)
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass)
		AverageSpeed_SimilarMass_Symm_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot, Distance_normalized_label, 
																Average_speed_label, Reldiff_label_avgspeed, self.Symm_legends, self.Plot_colors_symm,
																Primerror=Mean_stds, Secerror=Yerror, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Mass_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False, 
																ylim=special_limit)
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Mass_ranges, self.Filament_masses, 'Mass', binnum)
		Yplot2, Yerror2 = get_data_reldiffs(Fofr_only-1, Mass_ranges, RelDiff_AvgSpeed_SimilarMass, PropErr_AvgSpeed_SimilarMass)
		AverageSpeed_SimilarMass_fofr_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot2, Distance_normalized_label,
																Average_speed_label, Reldiff_label_avgspeed, self.fofr_legends, self.Plot_colors_fofr,
																Primerror=Mean_stds, Secerror=Yerror2, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Mass_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False,
																ylim=special_limit)
		if speedtype == 'Density':
			Density_reldiff_mass = Do_density_reldiff_gridspec(Common_bin_distances_normalized, Yplot, Yplot2, Yerror, Yerror2, self.Plot_colors_symm[1:], 
																self.Plot_colors_fofr[1:], self.Symm_legends[1:], self.fofr_legends[1:],
																self.Linestyles, (-0.05, 0.4), (-0.05, 0.4), Mass_titles)
			self.savefigure(Density_reldiff_mass, 'Reldiff_speed_similar_mass_ALL_gridspec', velocity_results_dir, dpi_mult=2)

		## Average speed of filaments with similar length
		if speedtype == 'Speed':
			specific_legend_axis = 0
		else:
			specific_legend_axis = -1
		RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength = store_reldiff_data(self.FilLengths, 'Length', Length_ranges, 
																						Common_bin_distances_normalized, binnum)
		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Length_ranges, self.FilLengths, 'Length', binnum)
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Length_ranges, RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength)
		AverageSpeed_SimilarLength_Symm_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot, Distance_normalized_label, 
																Average_speed_label, Reldiff_label_avgspeed, self.Symm_legends, self.Plot_colors_symm,
																Primerror=Mean_stds, Secerror=Yerror, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Length_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False)
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Length_ranges, self.FilLengths, 'Length', binnum)
		Yplot2, Yerror2 = get_data_reldiffs(Fofr_only-1, Length_ranges, RelDiff_AvgSpeed_SimilarLength, PropErr_AvgSpeed_SimilarLength)
		AverageSpeed_SimilarLength_fofr_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot2, Distance_normalized_label, 
																Average_speed_label,Reldiff_label_avgspeed, self.fofr_legends, self.Plot_colors_fofr,
																Primerror=Mean_stds, Secerror=Yerror2, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Length_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False,
																LegendAxis=specific_legend_axis)
		if speedtype == 'Density':
			Density_reldiff_length = Do_density_reldiff_gridspec(Common_bin_distances_normalized, Yplot, Yplot2, Yerror, Yerror2, self.Plot_colors_symm[1:], 
																self.Plot_colors_fofr[1:], self.Symm_legends[1:], self.fofr_legends[1:], 
																self.Linestyles, (-0.05, 0.4), (-0.05, 0.3), Length_titles)
			self.savefigure(Density_reldiff_length, 'Reldiff_speed_similar_length_ALL_gridspec', velocity_results_dir, dpi_mult=2)
		## Average speed of filaments with similar thickness
		RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness = store_reldiff_data(self.Thresholds, 'Thickness', Thickness_ranges, 
																						Common_bin_distances_normalized, binnum)
		Mean_profiles, Mean_stds = get_data(SymmLCDM, Symm_filenames, Common_bin_distances_normalized, Thickness_ranges, self.Thresholds, 'Thickness', binnum)
		Yplot, Yerror = get_data_reldiffs(Symm_only-1, Thickness_ranges, RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness)
		AverageSpeed_SimilarThickness_Symm_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot, Distance_normalized_label, 
																Average_speed_label, Reldiff_label_avgspeed, self.Symm_legends, self.Plot_colors_symm,
																Primerror=Mean_stds, Secerror=Yerror, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Thickness_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False)
		Mean_profiles, Mean_stds = get_data(FofrLCDM, Fofr_filenames, Common_bin_distances_normalized, Thickness_ranges, self.Thresholds, 'Thickness', binnum)
		Yplot2, Yerror2 = get_data_reldiffs(Fofr_only-1, Thickness_ranges, RelDiff_AvgSpeed_SimilarThickness, PropErr_AvgSpeed_SimilarThickness)
		AverageSpeed_SimilarThickness_fofr_GRIDSPEC = pf.Do_gridspec_sameX(Common_bin_distances_normalized, Mean_profiles, Yplot2, Distance_normalized_label, 
																Average_speed_label,Reldiff_label_avgspeed, self.fofr_legends, self.Plot_colors_fofr,
																Primerror=Mean_stds, Secerror=Yerror2, linestyles=self.Linestyles, reldiff=True,
																fillbetween=True, rowcol=[2,3], title=Thickness_titles, ylim_diff=ylimits,
																xscale=GS_xscale, yscale=GS_yscale, xscale_diff=GS_xscale, legend_anchor=False)
		if speedtype == 'Density':
			Density_reldiff_thickness = Do_density_reldiff_gridspec(Common_bin_distances_normalized, Yplot, Yplot2, Yerror, Yerror2, 
																self.Plot_colors_symm[1:], self.Plot_colors_fofr[1:], self.Symm_legends[1:], 
																self.fofr_legends[1:], self.Linestyles, (-0.1, 0.4), (-0.1, 0.4), Thickness_titles)
			self.savefigure(Density_reldiff_thickness, 'Reldiff_speed_similar_thickness_ALL_gridspec', velocity_results_dir, dpi_mult=2)

		####### Save figures #######
		print '--- SAVING IN: ', velocity_results_dir, ' ---'
		self.savefigure(AverageSpeed_AllFils, 'Average_speed_all_filaments', velocity_results_dir)
		### SIMILAR MASS
		self.savefigure(AverageSpeed_SimilarMass_Symm, 'Average_speed_similar_mass_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarMass_fofr, 'Average_speed_similar_mass_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarMass_plot_Symm, 'RelDiff_speed_similar_mass_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarMass_plot_fofr, 'RelDiff_speed_similar_mass_cFofr', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarMass_Symm_logx, 'Average_speed_similar_mass_logx_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarMass_fofr_logx, 'Average_speed_similar_mass_logx_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarMass_plot_Symm_logx, 'RelDiff_speed_similar_mass_logx_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarMass_plot_fofr_logx, 'RelDiff_speed_similar_mass_logx_cFofr', velocity_results_dir)
		### SIMILAR LENGTH
		self.savefigure(AverageSpeed_SimilarLength_Symm, 'Average_speed_similar_length_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarLength_fofr, 'Average_speed_similar_length_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarLen_plot_Symm, 'RelDiff_speed_similar_length_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarLen_plot_fofr, 'RelDiff_speed_similar_length_cFofr', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarLength_Symm_logx, 'Average_speed_similar_length_logx_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarLength_fofr_logx, 'Average_speed_similar_length_logx_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarLen_plot_Symm_logx, 'RelDiff_speed_similar_length_logx_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarLen_plot_fofr_logx, 'RelDiff_speed_similar_length_logx_cFofr', velocity_results_dir)
		### SIMILAR THICKNESS
		self.savefigure(AverageSpeed_SimilarThickness_Symm, 'Average_speed_similar_thickness_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarThickness_fofr, 'Average_speed_similar_thickness_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarThickness_plot_Symm, 'RelDiff_speed_similar_thickness_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarThickness_plot_fofr, 'RelDiff_speed_similar_thickness_cFofr', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarThickness_Symm_logx, 'Average_speed_similar_thickness_logx_cSymmetron', velocity_results_dir)
		self.savefigure(AverageSpeed_SimilarThickness_fofr_logx, 'Average_speed_similar_thickness_logx_cFofr', velocity_results_dir)
		self.savefigure(RelDiff_SimilarThickness_plot_Symm_logx, 'RelDiff_speed_similar_thickness_logx_cSymmetron', velocity_results_dir)
		self.savefigure(RelDiff_SimilarThickness_plot_fofr_logx, 'RelDiff_speed_similar_thickness_logx_cFofr', velocity_results_dir)
		### DIFFERENT MASS BINS
		self.savefigure(AverageSpeed_MassBins, 'Average_speed_massbins', velocity_results_dir)
		self.savefigure(AverageSpeed_RelativeDifference_MassBins, 'Reldiff_AverageSpeed_massbins', velocity_results_dir)
		self.savefigure(Average_speed_massbins_both, 'Average_speed_massbins_gridspec', velocity_results_dir)
		self.savefigure(AverageSpeed_LengthBins, 'Average_speed_lengthbins', velocity_results_dir)
		self.savefigure(AverageSpeed_RelativeDifference_LengthBins, 'Reldiff_AverageSpeed_lengthbins', velocity_results_dir)
		self.savefigure(Average_speed_thicknessbins_both, 'Average_speed_thicknessbins_gridspec', velocity_results_dir)
		self.savefigure(AverageSpeed_ThicknessBins, 'Average_speed_thicknessbins', velocity_results_dir)
		self.savefigure(AverageSpeed_RelativeDifference_ThicknessBins, 'Reldiff_AverageSpeed_thicknessbins', velocity_results_dir)
		self.savefigure(Average_speed_lengthbins_both, 'Average_speed_lengthbins_gridspec', velocity_results_dir)
		### Gridspec plots
		self.savefigure(AverageSpeed_SimilarMass_Symm_GRIDSPEC, 'Average_speed_similar_mass_cSymmetron_gridspec', velocity_results_dir, dpi_mult=2)
		self.savefigure(AverageSpeed_SimilarMass_fofr_GRIDSPEC, 'Average_speed_similar_mass_cFofr_gridspec', velocity_results_dir, dpi_mult=2)
		self.savefigure(AverageSpeed_SimilarLength_Symm_GRIDSPEC, 'Average_speed_similar_length_cSymmetron_gridspec', velocity_results_dir, dpi_mult=2)
		self.savefigure(AverageSpeed_SimilarLength_fofr_GRIDSPEC, 'Average_speed_similar_length_cFofr_gridspec', velocity_results_dir, dpi_mult=2)
		self.savefigure(AverageSpeed_SimilarThickness_Symm_GRIDSPEC, 'Average_speed_similar_thickness_cSymmetron_gridspec', velocity_results_dir, dpi_mult=2)
		self.savefigure(AverageSpeed_SimilarThickness_fofr_GRIDSPEC, 'Average_speed_similar_thickness_cFofr_gridspec', velocity_results_dir, dpi_mult=2)
		plt.close('all')	# Clear all current windows to free memory

	def Other_profiles(self):
		""" Other plots """
		#### Computing relevant stuff
		NModels = len(self.Filament_masses)
		Number_distribution_points = 100
		#Mass_label = 'Filament mass - [$M_\odot / h$]'
		#Number_label = '$N$ filaments'
		#Thickness_label = 'Filament thickness - [Mpc/h]'
		#Length_label = 'Filament length - [Mpc/h]'
		#Mean_Mass_label = 'Mean filament mass - [$M_\odot / h$]'
		#Mean_Thickness_label = 'Mean filament thickness - [Mpc/h]'
		
		Mass_label = '$M$ - [$M_\odot / h$]'
		Number_label = '$N$'
		Thickness_label = '$R_T$ - [Mpc/$h$]'
		Length_label = '$L$ - [Mpc/$h$]'
		Mean_Mass_label = r'$\bar{M} - [\mathrm{Mpc}/h]$'
		Mean_Thickness_label = r'$\bar{R_T} - [\mathrm{Mpc}/h]$'
		Distance_normalized_label = '$r/R_T$'

		Mass_N_label = '$N(>M)$'
		Thickness_N_label = '$N(>R_T)$'
		#SymmLCDM = np.array([0,1,2,3,4])
		#FofrLCDM = np.array([0,5,6,7])
		#Symm_only = np.array([1,2,3,4])
		#Fofr_only = np.array([5,6,7])
		SymmLCDM = self.SymmLCDM_ids
		FofrLCDM = self.FofrLCDM_ids
		Symm_only = self.SymmLCDM_ids[1:]
		Fofr_only = self.FofrLCDM_ids[1:]
		Symm_filenames = ['LCDM', 'Symm_A', 'Symm_B', 'Symm_C', 'Symm_D']
		Fofr_filenames = ['LCDM', 'fofr4', 'fofr5', 'fofr6']
		All_filenames = ['LCDM', 'Symm_A', 'Symm_B', 'Symm_C', 'Symm_D', 'fofr4', 'fofr5', 'fofr6']
		
		#### Number of filaments larger than a given mass N(>M). Also includes relative differences with LCDM as base model
		Mass_values = OF.Get_common_bin_logX(self.Filament_masses, binnum=Number_distribution_points)
		Mass_distribution = []
		Mass_distribution_error = []
		for Model_masses in self.Filament_masses:
			temp = []
			temp_err = []
			for masses in Mass_values:
				Number_count = len(np.where(Model_masses > masses)[0])
				Error = np.sqrt(Number_count)
				temp.append(float(Number_count))
				temp_err.append(Error)
			Mass_distribution.append(np.array(temp))
			Mass_distribution_error.append(np.array(temp_err))
		RelDiff_mass_distribution = [OF.relative_deviation_singular(Mass_distribution[0], Mass_distribution[i]) for i in range(1,NModels)]
		PropErr_mass_distribution = [OF.Propagate_error_reldiff(Mass_distribution[0], Mass_distribution[i], Mass_distribution_error[0],
															 Mass_distribution_error[i]) for i in range(1, NModels)]

		#### Number of filaments larger than a given thickness N(>T). Also includes relative differences with LCDM as base model
		Thickness_values = OF.Get_common_bin_logX(self.Thresholds, binnum=Number_distribution_points)
		Thickness_distribution = []
		Thickness_distribution_error = []
		for Model_thickness in self.Thresholds:
			temp = []
			temp_err = []
			for thickness in Thickness_values:
				Number_count = len(np.where(Model_thickness > thickness)[0])
				Error = np.sqrt(Number_count)
				temp.append(float(Number_count))
				temp_err.append(Error)
			Thickness_distribution.append(np.array(temp))
			Thickness_distribution_error.append(np.array(temp_err))
		RelDiff_thickness_distribution = [OF.relative_deviation_singular(Mass_distribution[0], Mass_distribution[i]) for i in range(1,NModels)]
		PropErr_thickness_distribution = [OF.Propagate_error_reldiff(Thickness_distribution[0], Thickness_distribution[i], Thickness_distribution_error[0],
															 Thickness_distribution_error[i]) for i in range(1, NModels)]
		######## Plotting stuff ########
		########
		######## MASS DISTRIBUTION
		########
		### Number of filaments larger than a given mass N(>M)
		Number_filaments_larger_mass = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(1,2,1)
		for i in SymmLCDM:
			plt.plot(Mass_values, Mass_distribution[i], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i])
			plt.fill_between(Mass_values, Mass_distribution[i]-Mass_distribution_error[i], Mass_distribution[i]+Mass_distribution_error[i],
							 alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		plt.yscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Mass_values, Mass_distribution[i], color=self.Plot_colors_fofr[ij], linestyle=self.Linestyles[ij])
			plt.fill_between(Mass_values, Mass_distribution[i]-Mass_distribution_error[i], Mass_distribution[i]+Mass_distribution_error[i],
							 alpha=0.4, facecolor=self.Plot_colors_fofr[ij])
		plt.legend(self.fofr_legends)
		plt.xscale('log')
		plt.yscale('log')
		Number_filaments_larger_mass.text(0.5, 0.01, Mass_label, ha='center', fontsize=10)
		Number_filaments_larger_mass.text(0.02, 0.55, Mass_N_label, ha='center', rotation='vertical', fontsize=10)
		### Relative differences of the data above
		Number_filaments_larger_mass_reldiff = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(2,1,1)
		plt.setp(ax.get_xticklabels(), visible=False)
		plt.plot(Mass_values, np.zeros(len(Mass_values)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Symm_only:
			plt.plot(Mass_values, RelDiff_mass_distribution[i-1], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i-1], label=self.Symm_legends[i])
			plt.fill_between(Mass_values, RelDiff_mass_distribution[i-1]-PropErr_mass_distribution[i-1],
							 RelDiff_mass_distribution[i-1]+PropErr_mass_distribution[i-1], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.xscale('log')
		ax2 = plt.subplot(2,1,2, sharey=ax)
		plt.plot(Mass_values, np.zeros(len(Mass_values)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Fofr_only:
			plt.plot(Mass_values, RelDiff_mass_distribution[i-1], color=self.Plot_colors_fofr[i-4], linestyle=self.Linestyles[i-5], label=self.fofr_legends[i-4])
			plt.fill_between(Mass_values, RelDiff_mass_distribution[i-1]-PropErr_mass_distribution[i-1],
							 RelDiff_mass_distribution[i-1]+PropErr_mass_distribution[i-1], alpha=0.4, facecolor=self.Plot_colors_fofr[i-4])
		plt.xscale('log')
		#ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True, fontsize=9)
		#ax2.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True, fontsize=9)
		h0,l0=ax.get_legend_handles_labels()
		ax.legend(h0,l0, fontsize=9)
		h2,l2=ax2.get_legend_handles_labels()
		ax2.legend(h2,l2, fontsize=9)
		Number_filaments_larger_mass_reldiff.text(0.5, 0, Mass_label, ha='center', fontsize=10)
		Number_filaments_larger_mass_reldiff.text(0, 0.5, r'$N_i(>M)/N_{\Lambda \mathrm{CDM}}(>M) - 1$', va='center', 
												ha='center', rotation='vertical', fontsize=10)
		ax.set_ylim(-0.1,0.5)
		ax2.set_ylim(-0.1,0.5)
		plt.tight_layout()

		########
		######## THICKNESS DISTRIBUTION
		########
		### Number of filaments larger than a given mass N(>M)
		Number_filaments_larger_thickness = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(1,2,1)
		for i in SymmLCDM:
			plt.plot(Thickness_values, Thickness_distribution[i], color=self.Plot_colors_symm[i], linestyle=self.Linestyles[i])
			plt.fill_between(Thickness_values, Thickness_distribution[i]-Thickness_distribution_error[i],
							 Thickness_distribution[i]+Thickness_distribution_error[i], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.legend(self.Symm_legends)
		plt.xscale('log')
		plt.yscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for ij in range(len(FofrLCDM)):
			i = FofrLCDM[ij]
			plt.plot(Thickness_values, Thickness_distribution[i], color=self.Plot_colors_fofr[ij], linestyle=self.Linestyles[ij])
			plt.fill_between(Thickness_values, Thickness_distribution[i]-Thickness_distribution_error[i],
							 Thickness_distribution[i]+Thickness_distribution_error[i], alpha=0.4, facecolor=self.Plot_colors_fofr[ij])
		plt.legend(self.fofr_legends)
		plt.xscale('log')
		plt.yscale('log')
		Number_filaments_larger_thickness.text(0.5, 0.01, Thickness_label, ha='center', fontsize=10)
		Number_filaments_larger_thickness.text(0.02, 0.55, Thickness_N_label, ha='center', rotation='vertical', fontsize=10)
		### Relative differences of the data above
		Number_filaments_larger_thickness_reldiff = plt.figure()
		plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
		ax = plt.subplot(2,1,1)
		plt.setp(ax.get_xticklabels(), visible=False)
		plt.plot(Thickness_values, np.zeros(len(Thickness_values)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Symm_only:
			plt.plot(Thickness_values, RelDiff_thickness_distribution[i-1], color=self.Plot_colors_symm[i], 
					 linestyle=self.Linestyles[i-1], label=self.Symm_legends[i])
			plt.fill_between(Thickness_values, RelDiff_thickness_distribution[i-1]-PropErr_thickness_distribution[i-1],
							 RelDiff_thickness_distribution[i-1]+PropErr_thickness_distribution[i-1], alpha=0.4, facecolor=self.Plot_colors_symm[i])
		plt.xscale('log')
		ax2 = plt.subplot(2,1,2, sharey=ax)
		plt.plot(Thickness_values, np.zeros(len(Thickness_values)), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), alpha=0.6)
		for i in Fofr_only:
			plt.plot(Thickness_values, RelDiff_thickness_distribution[i-1], color=self.Plot_colors_fofr[i-4],
					 linestyle=self.Linestyles[i-5], label=self.fofr_legends[i-4])
			plt.fill_between(Thickness_values, RelDiff_thickness_distribution[i-1]-PropErr_thickness_distribution[i-1],
							 RelDiff_thickness_distribution[i-1]+PropErr_thickness_distribution[i-1], alpha=0.4, facecolor=self.Plot_colors_fofr[i-4])
		plt.xscale('log')
		#ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True, fontsize=9)
		#ax2.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.3), ncol=1, fancybox=True, fontsize=9)
		h0,l0=ax.get_legend_handles_labels()
		ax.legend(h0,l0, fontsize=9)
		h2,l2=ax2.get_legend_handles_labels()
		ax2.legend(h2,l2, fontsize=9)
		Number_filaments_larger_thickness_reldiff.text(0.5, 0, Thickness_label, ha='center', fontsize=10)
		Number_filaments_larger_thickness_reldiff.text(0, 0.5, r'$N_i(>R_T)/N_{\Lambda \mathrm{CDM}}(>R_T) - 1$', va='center', 
												ha='center', rotation='vertical', fontsize=10)
		ax.set_ylim(-0.1,0.5)
		ax2.set_ylim(-0.1,0.8)
		plt.tight_layout()


		print '--- SAVING IN: ', self.results_dir, ' ---'
		###### Mass distribution plots
		self.savefigure(Number_filaments_larger_mass, 'Mass_distribution_number')
		self.savefigure(Number_filaments_larger_mass_reldiff, 'Mass_distribution_number_RelDiff')
		self.savefigure(Number_filaments_larger_thickness, 'Thickness_distribution_number')
		self.savefigure(Number_filaments_larger_thickness_reldiff, 'Thickness_distribution_number_RelDiff')
		
		

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
	parser.add_argument("-DoFilter", "--DoFilter", help="If True, further filters filament based on mass, thickness and length. Default = False",
						 type=bool, default=False)
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
	DoFilter = parsed_arguments.DoFilter
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
	Dist_accepted_sorted = []

	Filament_lengths = []
	All_speed_list = []
	Par_speed_list = []
	Orth_speed_list = []
	Density_prof = []

	N_filament_connections = []

	def Append_data(Fid, D_thres, OK_part, OK_dist, modelname):
		""" Appends common data to a common list """
		Filament_ids.append(Fid)
		Dist_thresholds.append(D_thres)
		Part_accepted.append(OK_part)
		Dist_accepted.append(OK_dist)
		Models_included.append(modelname)
		Dist_accepted_sorted_temp = []
		for i in range(len(OK_dist)):
			Dist_accepted_sorted_temp.append(np.sort(OK_dist[i]))
		Dist_accepted_sorted.append(np.array(Dist_accepted_sorted_temp))

	def Append_data_speeds(AllS, OrthS, ParaS):
		All_speed_list.append(AllS)
		Par_speed_list.append(ParaS)
		Orth_speed_list.append(OrthS)

	if N_parts == 64:
		for p_model in Models_to_be_run:
			if not p_model == 'symmC':
				Instance = FilterParticlesAndFilaments(p_model, N_parts, N_sigma)
				OK_fils, thresholds, OK_particles, OK_distances, SegIDs = Instance.Get_threshold_and_noise()
				Append_data(OK_fils, thresholds, OK_particles, OK_distances, p_model)
				Speeds, Ospeed, Pspeed = Instance.Get_speed_components(OK_particles, SegIDs, OK_fils)
				Append_data_speeds(Speeds, Ospeed, Pspeed)
				Fillens = Instance.Get_filament_length()[OK_fils]
				Filament_lengths.append(Fillens)
				Densities = Instance.Compute_density_profile(OK_distances, Fillens)
				Density_prof.append(Densities*Mpc**3/Solmass*1e-14)	# Convert from kg h^2/m^3 to M_sun h^2/Mpc^3
				N_filament_connections.append(Instance.Number_filament_connections()[OK_fils])		
	else:
		for p_model in Models_to_be_run:
			#if p_model == 'lcdm' or p_model == 'symmA' or p_model == 'fofr4':
			Instance = FilterParticlesAndFilaments(p_model, N_parts, N_sigma)
			OK_fils, thresholds, OK_particles, OK_distances, SegIDs = Instance.Get_threshold_and_noise()
			Speeds, Ospeed, Pspeed = Instance.Get_speed_components(OK_particles, SegIDs, OK_fils)
			Fillens = Instance.Get_filament_length()[OK_fils]
			Densities = Instance.Compute_density_profile(OK_distances, Fillens)
			print "Number of filaments after density and length filtering 1:", len(Fillens)
			if DoFilter:
				# Do a secondary filter	of the length
				#Large_filament_filt = Instance.Secondary_filter(OK_fils)
				Large_filament_filt = Fillens <= 20
				# Filter out thickness larger than 10
				Large_thickness_filt = thresholds < 10
				Total_filt = Large_filament_filt*Large_thickness_filt
				OK_fils = OK_fils[Total_filt]
				thresholds = thresholds[Total_filt]
				OK_particles = OK_particles[Total_filt]
				OK_distances = OK_distances[Total_filt]
				SegIDs = SegIDs[Total_filt]
				Fillens = Fillens[Total_filt]

				Speeds = Speeds[Total_filt]
				Ospeed = Ospeed[Total_filt]
				Pspeed = Pspeed[Total_filt]
				Densities = Densities[Total_filt]
			Append_data(OK_fils, thresholds, OK_particles, OK_distances, p_model)
			Append_data_speeds(Speeds, Ospeed, Pspeed)
			Filament_lengths.append(Fillens)
			Density_prof.append(Densities*Mpc**3/Solmass*1e-14)	# Convert from kg h^2/m^3 to M_sun h^2/Mpc^3
			N_filament_connections.append(Instance.Number_filament_connections()[OK_fils])
			print "Number of filaments after filter for " + p_model + ":", len(Fillens)

	#Plot_instance = Plot_results(Models_included, N_sigma, 'ModelComparisons/ParticleAnalysis/', filetype=Filetype)
	#Plot_instance.Particle_profiles(Dist_thresholds, Part_accepted, Filament_lengths)
	#Plot_instance.Velocity_profiles(All_speed_list, Dist_accepted, speedtype='Speed')
	#Plot_instance.Velocity_profiles(Orth_speed_list, Dist_accepted, speedtype='Orthogonal')
	#Plot_instance.Velocity_profiles(Par_speed_list, Dist_accepted, speedtype='Parallel')
	#Plot_instance.Velocity_profiles(Density_prof, Dist_accepted_sorted, speedtype='Density')
	#Plot_instance.Other_profiles()

	#savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
	#CompI = HComp.CompareModels(savefile=1, foldername='ModelComparisons/FilteredGlobalProperties/',
	#							 savefile_directory=savefile_directory, filetype='.pdf', nPart=N_parts, Nsigma=N_sigma)
	#CompI.Compare_disperse_data_clean(N_filament_connections, Filament_lengths, [])