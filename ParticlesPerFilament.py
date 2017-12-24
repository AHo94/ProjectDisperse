import numpy as np
import ReadGadgetFile
import os
import cPickle as pickle

# Global variables
# Box boundaries assumed to be from 0 to 256.0 Mpc/h
lower_boundary = 0.0
upper_boundary = 256.0
class particles_per_filament():
	"""	A class that computes the distance of the dark matter particles to the filaments. """
	def __init__(self, model, maskdirs, masklimits, box_expand):
		self.model = model
		RGF_instance = ReadGadgetFile.Read_Gadget_file(maskdirs, masklimits)
		self.particlepos = RGF_instance.Get_3D_particles(model)
		self.box_expand = box_expand

	def filament_box(self, filament):
		""" 
		Determines the max and min x,y and z coordinates based on the filament.
		The bounding box is defined by these coordinates.
		"""
		xmin = np.min(filament[:,0])
		xmax = np.max(filament[:,0])
		ymin = np.min(filament[:,1])
		ymax = np.max(filament[:,1])
		zmin = np.min(filament[:,2])
		zmax = np.max(filament[:,2])
		return xmin, xmax, ymin, ymax, zmin, zmax

	def particle_mask(self, box, ParticlePos):
		""" Mask the dark matter particles based on the box """
		xmin = box[0]
		xmax = box[1]
		ymin = box[2]
		ymax = box[3]
		zmin = box[4]
		zmax = box[5]
		maskx = (ParticlePos[:,0] > xmin) & (ParticlePos[:,0] < xmax)
		masky = (ParticlePos[:,1] > ymin) & (ParticlePos[:,1] < ymax)
		maskz = (ParticlePos[:,2] > zmin) & (ParticlePos[:,2] < zmax)
		return maskx*masky*maskz

	def particle_box(self, filament_box, distance_threshold):
		"""
		Creates a particle box based on the filament box. 
		The size is increased a little to include more particles. Size increase can be changed.
		If the box increase surpasses the simulation box, the box will 'continue' on the other side.
		This will create two seperate boxes. Particles in the second box will be moved to the other boundary.

		Example: Filament is near the left edge of the x-axis, i.e close to xmin.
		The box_expand causes xmin - box_expand to be less than the boundary, i.e xmin - box_expand < 0.
		We create a second box at the other side of the boundary, i.e xmax, with the size corresponding of whats left of the difference xmin - lower_boundary.
		Particles in the other box, the ones close to xmax, will be moved so they are close to xmin.
		See notes for more details.
		"""
		xmin = filament_box[0] - self.box_expand
		xmax = filament_box[1] + self.box_expand
		ymin = filament_box[2] - self.box_expand
		ymax = filament_box[3] + self.box_expand
		zmin = filament_box[4] - self.box_expand
		zmax = filament_box[5] + self.box_expand
		box = [xmin, xmax, ymin, ymax, zmin, zmax]
		box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
		MovePartx = 0
		MoveParty = 0
		MovePartz = 0
		Atboundary = 0
		if np.abs((xmin + self.box_expand) - lower_boundary) < 1e-3:
			box[0] = xmin + self.box_expand
			Atboundary = 1
		elif np.abs((xmax - self.box_expand) - upper_boundary) < 1e-3:
			box[1] = xmax - self.box_expand
			Atboundary = 1
		if np.abs((ymin + self.box_expand) - lower_boundary) < 1e-3:
			box[2] = ymin + self.box_expand
			Atboundary = 1
		elif np.abs((ymax - self.box_expand) - upper_boundary) < 1e-3:
			box[3] = ymax - self.box_expand
			Atboundary = 1
		if np.abs((zmin + self.box_expand) - lower_boundary) < 1e-3:
			box[4] = zmin + self.box_expand
			Atboundary = 1
		elif np.abs((zmax - self.box_expand) - upper_boundary) < 1e-3:
			box[5] = zmax - self.box_expand
			Atboundary = 1
			
		if Atboundary:
			box2 = box
		else:
			if xmin < lower_boundary:
				box[0] = lower_boundary
				box2[1] = upper_boundary
				box2[0] = upper_boundary + (xmin - lower_boundary)
				MovePartx = -256.0
			elif xmax > upper_boundary:
				box[1] = upper_boundary
				box2[0] = lower_boundary
				box2[1] = lower_boundary + (xmax - upper_boundary)
				MovePartx = 256.0
				
			if ymin < lower_boundary:
				box[2] = lower_boundary
				box2[3] = upper_boundary
				box2[2] = upper_boundary + (ymin - lower_boundary)
				MoveParty = -256.0
			elif ymax > upper_boundary:
				box[3] = upper_boundary
				box2[2] = lower_boundary
				box2[3] = lower_boundary + (ymax - upper_boundary)
				MoveParty = 256.0
				
			if zmin < lower_boundary:
				box[4] = lower_boundary
				box2[5] = upper_boundary
				box2[4] = upper_boundary + (zmin - lower_boundary)
				MovePartz = -256.0
			elif zmax > upper_boundary:
				box[5] = upper_boundary
				box2[4] = lower_boundary
				box2[5] = lower_boundary + (zmax - upper_boundary)
				MovePartz = 256.0

		if box == box2:
			mask = self.particle_mask(box, self.particlepos)
			return self.particlepos[mask],1# self.particleIDs[mask]
		else:
			mask1 = self.particle_mask(box, self.particlepos)
			mask2 = self.particle_mask(box2, self.particlepos)
			Particles1 = self.particlepos[mask1]
			Particles2 = self.particlepos[mask2]
			#PartIDs1 = self.particleIDs[mask1]
			#PartIDs2 = self.particleIDs[mask2]
			Particles2[:,0] += MovePartx
			Particles2[:,1] += MoveParty
			Particles2[:,2] += MovePartz
			Masked_particles = np.concatenate([Particles1, Particles2])
			Masked_particleIDs = np.concatenate([PartIDs1, PartIDs2])
			return Masked_particles,1#  Masked_particleIDs

	def particle_box2(self, filament, masked_ids):
		"""
		Creates a particle box based on the filament box. 
		The size is increased a little to include more particles. Size increase can be changed.
		If the box increase surpasses the simulation box, the box will 'continue' on the other side.
		This will create two seperate boxes. Particles in the second box will be moved to the other boundary.

		Example: Filament is near the left edge of the x-axis, i.e close to xmin.
		The box_expand causes xmin - box_expand to be less than the boundary, i.e xmin - box_expand < 0.
		We create a second box at the other side of the boundary, i.e xmax, with the size corresponding of whats left of the difference xmin - lower_boundary.
		Particles in the other box, the ones close to xmax, will be moved so they are close to xmin.
		See notes for more details.
		"""
		xmin = np.min(filament[:,0]) - self.box_expand
		xmax = np.max(filament[:,0]) + self.box_expand
		ymin = np.min(filament[:,1]) - self.box_expand
		ymax = np.max(filament[:,1]) + self.box_expand
		zmin = np.min(filament[:,2]) - self.box_expand
		zmax = np.max(filament[:,2]) + self.box_expand
		MovePartx = 0
		MoveParty = 0
		MovePartz = 0
		if xmin < lower_boundary:
			MovePartx = -256.0
		elif xmax > upper_boundary:
			MovePartx = 256.0
				
		if ymin < lower_boundary:
			MoveParty = -256.0
		elif ymax > upper_boundary:
			MoveParty = 256.0
				
		if zmin < lower_boundary:
			MovePartz = -256.0
		elif zmax > upper_boundary:
			MovePartz = 256.0

		if not len(masked_ids) == 2:
			return self.particlepos[masked_ids]
		else:
			Particles1 = self.particlepos[masked_ids[0]]
			Particles2 = self.particlepos[masked_ids[1]]
			Particles2[:,0] += MovePartx
			Particles2[:,1] += MoveParty
			Particles2[:,2] += MovePartz
			Masked_particle_box = np.concatenate([Particles1, Particles2])
			return Masked_particles_box

	def masked_particle_indices(self, filament_box):
		"""
		Creates a mask of the particles based on the filament box.
		The function will return particle IDs based on the corresponding mask.
		The particle IDs will correspond to the particle box.
		This function does NOT create a particle box for computation.
		"""
		xmin = filament_box[0] - self.box_expand
		xmax = filament_box[1] + self.box_expand
		ymin = filament_box[2] - self.box_expand
		ymax = filament_box[3] + self.box_expand
		zmin = filament_box[4] - self.box_expand
		zmax = filament_box[5] + self.box_expand
		box = [xmin, xmax, ymin, ymax, zmin, zmax]
		box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
		MovePartx = 0
		MoveParty = 0
		MovePartz = 0
		Atboundary = 0
		if np.abs((xmin + self.box_expand) - lower_boundary) < 1e-3:
			box[0] = xmin + self.box_expand
			Atboundary = 1
		elif np.abs((xmax - self.box_expand) - upper_boundary) < 1e-3:
			box[1] = xmax - self.box_expand
			Atboundary = 1
		if np.abs((ymin + self.box_expand) - lower_boundary) < 1e-3:
			box[2] = ymin + self.box_expand
			Atboundary = 1
		elif np.abs((ymax - self.box_expand) - upper_boundary) < 1e-3:
			box[3] = ymax - self.box_expand
			Atboundary = 1
		if np.abs((zmin + self.box_expand) - lower_boundary) < 1e-3:
			box[4] = zmin + self.box_expand
			Atboundary = 1
		elif np.abs((zmax - self.box_expand) - upper_boundary) < 1e-3:
			box[5] = zmax - self.box_expand
			Atboundary = 1
			
		if Atboundary:
			box2 = box
		else:
			if xmin < lower_boundary:
				box[0] = lower_boundary
				box2[1] = upper_boundary
				box2[0] = upper_boundary + (xmin - lower_boundary)
			elif xmax > upper_boundary:
				box[1] = upper_boundary
				box2[0] = lower_boundary
				box2[1] = lower_boundary + (xmax - upper_boundary)
				
			if ymin < lower_boundary:
				box[2] = lower_boundary
				box2[3] = upper_boundary
				box2[2] = upper_boundary + (ymin - lower_boundary)
			elif ymax > upper_boundary:
				box[3] = upper_boundary
				box2[2] = lower_boundary
				box2[3] = lower_boundary + (ymax - upper_boundary)
				
			if zmin < lower_boundary:
				box[4] = lower_boundary
				box2[5] = upper_boundary
				box2[4] = upper_boundary + (zmin - lower_boundary)
			elif zmax > upper_boundary:
				box[5] = upper_boundary
				box2[4] = lower_boundary
				box2[5] = lower_boundary + (zmax - upper_boundary)

		if box == box2:
			mask = self.particle_mask(box, self.particlepos)
			return np.where(mask)[0]
		else:
			mask1 = self.particle_mask(box, self.particlepos)
			mask2 = self.particle_mask(box2, self.particlepos)
			indices1 = np.where(mask1)[0]
			indices2 = np.where(mask2)[0]
			#Masked_indices = np.concatenate([indices1, indices2])
			return np.array([indices1, indices2])

	def get_distance_old(self, filament, part_box):
		""" 
		Computes the distance from a particle to every segment in the filament.
		Only the shortest distance is included.
		Old version, use the ones below.
		"""
		distances = []
		checker = 0
		Filament_point_diff = np.linalg.norm(filament[1:] - filament[:-1], axis=1)
		for part_pos in part_box:
			diff1 = part_pos - filament[1:]
			diff2 = part_pos - filament[:-1]
			cross_prod = np.cross(diff1, diff2)
			cross_prod_norm = np.linalg.norm(cross_prod, axis=1)
			d = cross_prod_norm/Filament_point_diff
			distances.append(np.min(d))
		return np.array(distances)

	def get_distance(self, filament, masked_ids):
		""" 
		Computes the distance from a particle to every segment in the filament.
		Only the shortest distance is included.
		Input must be filament 3D position and the corresponding particle mask.
		This version loops through all particles in the masked particle box.
		Particles distances are thus individually calculated to each segment.
		Does not scale for large number of particles.
		"""
		part_box = self.particle_box2(filament, masked_ids)
		Filament_point_diff = np.linalg.norm(filament[1:] - filament[:-1], axis=1)
		distances = []
		for part_pos in part_box:
			diff1 = part_pos - filament[1:]
			diff2 = part_pos - filament[:-1]
			cross_prod = np.cross(diff1, diff2)
			cross_prod_norm = np.linalg.norm(cross_prod, axis=1)
			d = cross_prod_norm/Filament_point_diff
			distances.append(np.min(d))
		return np.array(distances)

	def get_distance_scaled(self, filament, part_box):
		"""
		Computes the distance from a particle to every segment in the filament.
		Only the shortest distance is included.
		Input must be filament 3D position and the corresponding particle masked box.
		This version loops through every segment and finds the distance from said segment to every particle in the masked box.
		Once all segemnts are looped over, we find the smallest distance from one particle to each segment.
		Scales a lot better for larger number of particles.
		"""
		dist_temp = []
		for i in range(len(filament)-1):
			d1 = part_box - filament[i]
			d2 = part_box - filament[i+1]
			cross_prod = np.cross(d1, d2)
			cross_prod_norm = np.linalg.norm(cross_prod, axis=1)
			d = cross_prod_norm/np.linalg.norm(filament[i] - filament[i+1])
			dist_temp.append(d)
		dist_temp = np.swapaxes(np.asarray(dist_temp), 0, 1)
		distances = np.min(dist_temp, axis=1)
		return distances


	def get_distance_seginterp(self, filament, part_box):
		"""
		For each segment, 'interpolate' a set of points between the two connection points in the segment.
		Create N new 3D coordinate points in the segment.
		For each point in the segment, compute distance from the point to each particle.
		The current output should be a (N x Pn)-matrix, where Pn = number of particles.
		Reshape the output such that it becomes a (Pn x N)-matrix.
		Find the minimum along the N axis of the matrix.
		This gives an Pn-array shortest distance for each particle to any of the interpolated points.
		Repeat this for other segments, will have multiple arrays of 'shortest distance'.
		Repeat until we have an (Sn x Pn)-matrix, where Sn = number of segments
		Reshape so we have an (Pn x Sn)-matrix, that is a matrix where each array contains distances of one particle to each seg.
		Find minimum of each of those to get shortest distance from particles to filament

		This version is used as the previous two versions assumed the segment to be infinitely long.
		"""
		true_dist = []
		for i in range(len(filament)-1):
			segpoints = []
			distances = []
			xlims = np.linspace(filament[i][0], filament[i+1][0], 100)
			ylims = np.linspace(filament[i][1], filament[i+1][1], 100)
			zlims = np.linspace(filament[i][2], filament[i+1][2], 100)
			# Create 3D position array of the segments, based on the interpolated values
			for i in range(len(xlims)):
				segpoints.append(np.column_stack((xlims[i], ylims[i], zlims[i])))
			# Find distance from each point in the segment to every particle
			for pts in segpoints:
				distances.append(np.linalg.norm(pts - part_box, axis=1))
			# Selects the shortest distance between a particle and every segment point.
			distances = np.swapaxes(np.asarray(distances), 0, 1)
			shortest = np.min(distances, axis=1)
			true_dist.append(shortest)
		# Selects the shortest distance from a particle to every segment
		true_dist = np.swapaxes(np.asarray(true_dist), 0, 1)
		dist = np.min(true_dist, axis=1)
		return dist


	def get_masked_ids(self, filament):
		filbox = self.filament_box(filament)
		masked_ids = self.masked_particle_indices(filbox)
		return masked_ids

	def solve(self, filament, distance_threshold):
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticleData/'
		cache_model = cachedir + self.model + '_particleIDs.p'

		####
		# Note to self, maybe return the mask instead of all Ids and lengths.
		# This will prevent sending a lot of data arouund multi processing.
		# One caveat is that we have to read the gadget file twice, but should not be an issue.
		# Returning True/False array mask would probably be excessive in terms of memory usage.
		# Try returning indices where the mask is True.
		# Masking is the time consuming part anyway.
		# See solve function above
		####
		#if os.path.isfile(cache_model):
		#	self.particleIDs = pickle.load(open(cache_model, 'rb'))
		#else:
		#	raise ValueError("Pickle file containing particle IDs does not exist!")

		filbox = self.filament_box(filament)
		partbox, masked_part_ids = self.particle_box(filbox, distance_threshold)
		distances = self.get_distance(filament, partbox)
		accepted_particle_ids = np.where(distances <= distance_threshold)[0]
		return len(accepted_particle_ids),1# masked_part_ids[accepted_particle_ids]

# Everything below is old code for -periodic argument in delaunaya
def get_filament_box(filament, filament_p2=np.array([])):
	Box_expand = 1
	if not filament_p2.any():
		xmin = np.min(filament[:,0]) - Box_expand
		xmax = np.max(filament[:,0]) + Box_expand
		ymin = np.min(filament[:,1]) - Box_expand
		ymax = np.max(filament[:,1]) + Box_expand
		zmin = np.min(filament[:,2]) - Box_expand
		zmax = np.max(filament[:,2]) + Box_expand
		box = [xmin, xmax, ymin, ymax, zmin, zmax]
		return [box]
	else:
		xmin = np.min(filament[:,0]) - Box_expand
		xmax = np.max(filament[:,0]) + Box_expand
		ymin = np.min(filament[:,1]) - Box_expand
		ymax = np.max(filament[:,1]) + Box_expand
		zmin = np.min(filament[:,2]) - Box_expand
		zmax = np.max(filament[:,2]) + Box_expand
		
		xmin_p2 = np.min(filament_p2[:,0]) - Box_expand
		xmax_p2 = np.max(filament_p2[:,0]) + Box_expand
		ymin_p2 = np.min(filament_p2[:,1]) - Box_expand
		ymax_p2 = np.max(filament_p2[:,1]) + Box_expand
		zmin_p2 = np.min(filament_p2[:,2]) - Box_expand
		zmax_p2 = np.max(filament_p2[:,2]) + Box_expand
		if (np.abs(np.min(filament[:,0]) - min_boundary) < 1e-3) or (np.abs(np.max(filament[:,0]) - max_boundary) < 1e-3):
			xmin = np.min(filament[:,0])
			xmax = np.max(filament[:,0])
			xmin_p2 = np.min(filament_p2[:,0])
			xmax_p2 = np.max(filament_p2[:,0])
		elif (np.abs(np.min(filament[:,1]) - min_boundary) < 1e-3) or (np.abs(np.max(filament[:,1]) - max_boundary) < 1e-3):
			ymin = np.min(filament[:,1])
			ymax = np.max(filament[:,1])
			ymin_p2 = np.min(filament_p2[:,1])
			ymax_p2 = np.max(filament_p2[:,1])
		elif (np.abs(np.min(filament[:,2]) - min_boundary) < 1e-3) or (np.abs(np.max(filament[:,2]) - max_boundary) < 1e-3):
			zmin = np.min(filament[:,2])
			zmax = np.max(filament[:,2])
			zmin_p2 = np.min(filament_p2[:,2])
			zmax_p2 = np.max(filament_p2[:,2])
			
		box1 = [xmin, xmax, ymin, ymax, zmin, zmax]
		box2 = [xmin_p2, xmax_p2, ymin_p2, ymax_p2, zmin_p2, zmax_p2]
		return [box1, box2]
