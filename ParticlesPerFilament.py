# Basic modules
import numpy as np
import os
import sys
import argparse
import cPickle as pickle

# ZMQ modules
import zmq
import ZMQArraySending as ZMQAS

# Own modules
import ReadGadgetFile


# Global variables
# Box boundaries assumed to be from 0 to 256.0 Mpc/h
lower_boundary = 0.0
upper_boundary = 256.0
halfbox = 256.0/2.0
class particles_per_filament():
	"""	
	A class that computes the distance of the dark matter particles to the filaments.
	Used with Python's Multiprocessing module. Check non-class function for ZMQ computing.
	"""
	def __init__(self, model, maskdirs, masklimits, box_expand):
		self.model = model
		RGF_instance = ReadGadgetFile.Read_Gadget_file(maskdirs, masklimits)
		self.particlepos = RGF_instance.Get_3D_particles(model)
		self.box_expand = np.float32(box_expand)

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
			return Masked_particle_box

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
			return (np.where(mask)[0]).astype(np.int32)
		else:
			mask1 = self.particle_mask(box, self.particlepos)
			mask2 = self.particle_mask(box2, self.particlepos)
			indices1 = (np.where(mask1)[0]).astype(np.int32)
			indices2 = (np.where(mask2)[0]).astype(np.int32)
			Masked_indices = np.concatenate([indices1, indices2])
			return Masked_indices
			#return np.array([indices1, indices2])

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

	def get_distance_seginterp(self, filament, masked_ids):
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
		#print 'Actually computing now'
		true_dist = []
		part_box = self.particle_box2(filament, masked_ids)
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
		#print 'done'
		return dist


	def get_masked_ids(self, filament):
		filbox = self.filament_box(filament)
		masked_ids = self.masked_particle_indices(filbox)
		return masked_ids

	def get_partbox_and_ids(self, filament):
		filbox = self.filament_box(filament)
		masked_ids = self.masked_particle_indices(filbox)
		partbox = self.particle_box2(filament, masked_ids)
		return masked_ids, partbox

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

def Argument_parser():
	""" Parses optional argument when program is run from the command line """
	parser = argparse.ArgumentParser()
	# Optional arguments
	parser.add_argument("-BoxExp", "--BOXEXPAND", help="Determines how far the filament masking box increases from the filament edges. Default = 3.0", type=float, default=3.0)
	parser.add_argument("-euclid21", "--Euclid21Run", help="Connects worker to localhost instead of infiniband if True. Use for workers in euclid21. False by default.", type=int, default=0)
	# Parse arguments
	args = parser.parse_args()
	return args

def get_interpolation_points(filament, BoxSize):
	"""
	Gives a set amount of interpolation points between the segments based on their length.
	The longer the segment is, the more interpolation points it gets.
	Total number of interpolation points is 100 per segment.
	Each segment will have at least 5 points.

	If the number of points in the end does not add up, after assigning each segment a set number of points, the last few points goes to the largest segment.
	If the numbeo of left over points are larger than the number of interpolated points for the largest segment, then it will be distributed evenly.
	"""
	nsegs = len(filament) - 1
	Segment_lengths = []
	diffx = filament[1:,0] - filament[:-1,0]
	diffy = filament[1:,1] - filament[:-1,1]
	diffz = filament[1:,2] - filament[:-1,2]
	diffx[diffx <= -BoxSize/2.0] += BoxSize
	diffx[diffx >= BoxSize/2.0] -= BoxSize
	diffy[diffy <= -BoxSize/2.0] += BoxSize
	diffy[diffy >= BoxSize/2.0] -= BoxSize
	diffz[diffz <= -BoxSize/2.0] += BoxSize
	diffz[diffz >= BoxSize/2.0] -= BoxSize
	for i in range(len(diffx)):
		Segment_lengths.append(np.sqrt(diffx[i]**2 + diffy[i]**2 + diffz[i]**2))

	numpts = nsegs*100
	number_points = np.ones(nsegs, dtype=np.int32)*5
	numpts_new = numpts - 5*nsegs
	
	Total_length = np.sum(Segment_lengths)
	Percentage = np.array(Segment_lengths)/Total_length
	numpts_percentage = np.round(Percentage*numpts_new).astype(np.int32)
	number_points += numpts_percentage
	total_number_points = np.sum(number_points)

	if not total_number_points == numpts:
		diff = numpts - total_number_points
		if np.abs(diff) > np.max(number_points):
			diff_share = diff/nsegs
			diff_left = diff - diff_share*nsegs
			number_points += diff_share
			left = numpts - np.sum(number_points)
			if left:
				if left < numpts:
					number_points[:left] += 1
				elif left > numpts:
					largest = np.where(numpts_percentage == np.max(numpts_percentage))[0][0]
					number_points[largest] += left
		else:
			largest = np.where(numpts_percentage == np.max(numpts_percentage))[0][0]
			number_points[largest] += diff
	return number_points

def particle_difference_periodic(segpoint, part_box):
	"""
	Computes the difference between the filament segment point and the particle point(s).
	Does not compute the distance between them. This is done outside the function with numpy.linalg.norm
	Assuming box size is the same as the particle box size, i.e 256.0 Mpc/h every direction
	"""
	Difference = segpoint - part_box
	diffx = Difference[:,0]
	diffy = Difference[:,1]
	diffz = Difference[:,2]
	Difference[diffx <= -halfbox,0] += 256.0
	Difference[diffx >= halfbox,0] -= 256.0
	Difference[diffy <= -halfbox,1] += 256.0
	Difference[diffy >= halfbox,1] -= 256.0
	Difference[diffz <= -halfbox,2] += 256.0
	Difference[diffz >= halfbox,2] -= 256.0
	return Difference   

def Compute_distance(filament, part_box, BoxSize):
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
	
	Also finds the index of the interpolated segment points, which corresponds to the axis along the filament.
	In the case where particles has the exact same distance to two segment points, choose the first one.
	
	This version is used as the previous two versions assumed the segment to be infinitely long.
	This version is also outside the class, to not load all the particle positions.
	"""
	true_dist = []
	filaxis_temp = []
	N_seg_pts = get_interpolation_points(filament, BoxSize)
	for i in range(len(filament)-1):
		segpoints = []
		distances = []
		xlims = np.linspace(filament[i][0], filament[i+1][0], N_seg_pts[i])
		ylims = np.linspace(filament[i][1], filament[i+1][1], N_seg_pts[i])
		zlims = np.linspace(filament[i][2], filament[i+1][2], N_seg_pts[i])
		# Create 3D position array of the segments, based on the interpolated values
		for j in range(len(xlims)):
			segpoints.append(np.column_stack((xlims[j], ylims[j], zlims[j])))
		# Find distance from each point in the segment to every particle
		for pts in segpoints:
			differences = particle_difference_periodic(pts, part_box)
			distances.append(np.linalg.norm(differences, axis=1))
		# Selects the shortest distance between a particle and every segment point.
		distances2 = np.swapaxes(np.asarray(distances), 0, 1)
		shortest = np.min(distances2, axis=1)
		true_dist.append(shortest)
		# Selects the corresponding coordinate index from the segpoints
		index_segpoint = []
		for k in range(len(shortest)):
			identical = np.where(np.abs(shortest[k] - distances2[k]) < 1e-16)[0]
			index_segpoint.append(identical[0])
		filaxis_temp.append(np.array(index_segpoint) + i*100)
	# Selects the shortest distance from a particle to every segment
	true_dist = np.swapaxes(np.asarray(true_dist), 0, 1)
	dist = np.min(true_dist, axis=1)
	# Finds the corresponding segpoint index
	indices = []
	for k in range(len(dist)):
		indices.append(np.where(np.abs(dist[k] - true_dist[k]) < 1e-16)[0])
	filaxis_swapped = np.swapaxes(np.asarray(filaxis_temp), 0, 1)
	Segpoint_index = []
	# Some particles may have similar distance to multiple segment points. Choose the first one
	for l in range(len(indices)):
		Segpoint_index.append(filaxis_swapped[l][indices[l][0]])
	return dist.astype(np.float32), np.array(Segpoint_index, np.int32)

def get_distance_analytic(filpoint1, filpoint2, part_box):
	""" 
	Computes the distance from a set of particles to a given segment 
	Let the segment be parametrized by s(t) = s_1 + t*(s_2 - s_1), where s_1 and s_2 are the filament point of a segment
	We only allow solutions of t to be t = [0, 1]
	Distance from particle position p to segment is D = |s(t) - p|
	To find a solution of t, solve the derivative of D^2 = 0, i.e. dD^2/dt = 0, which gives
	t_solution = [(p-s_1)*(s_2-s_1)]/|s_2 - s_1|^2, dot product in the numerator
	Ensure that t_solution is between 0 and 1 and use this to compute the distance from particle to segment
	Takes into account periodic boundary when the differences in 3D are computed
	"""
	segment = filpoint2 - filpoint1
	segment[segment >= halfbox] -= 256.0
	segment[segment <= -halfbox] += 256.0
	particle_dist = part_box - filpoint1
	particle_dist[particle_dist[:,0] >= halfbox, 0] -= 256.0
	particle_dist[particle_dist[:,0] <= -halfbox, 0] += 256.0
	particle_dist[particle_dist[:,1] >= halfbox, 1] -= 256.0
	particle_dist[particle_dist[:,1] <= -halfbox, 1] += 256.0
	particle_dist[particle_dist[:,2] >= halfbox, 2] -= 256.0
	particle_dist[particle_dist[:,2] <= -halfbox, 2] += 256.0
	#t_solution = np.dot(part_box - filpoint1, segment)/(np.linalg.norm(segment)**2.0)
	t_solution = np.dot(particle_dist, segment)/(np.linalg.norm(segment)**2.0)
	# Ensure the solution of t is between 0 and 1
	Smaller0 = t_solution < 0.0
	Greater1 = t_solution > 1.0
	t_solution[Smaller0] = 0.0
	t_solution[Greater1] = 1.0
	#distances = filpoint1 + t_solution.reshape(len(t_solution), 1)*segment - part_box
	distances = t_solution.reshape(len(t_solution), 1)*segment - particle_dist
	# Apply periodic boundary condition to distances
	distances[distances <= -halfbox] += 256.0
	distances[distances >= halfbox] -= 256.0
	return np.linalg.norm(distances, axis=1), t_solution
    
def Compute_distance_analytic(filament, part_box):
	""" 
	Analytical method of finding particle distance to a segment line
	Computes particle distance to each segment. Returns the shortest distance of all computed distances.
	Also returns the segment number the particle is closest to.
	"""
	true_dist = []
	filaxis_temp = []
	t_solutions_all = []
	Nsegs = len(filament)
	for i in range(Nsegs-1):
		distances, t_sol = get_distance_analytic(filament[i], filament[i+1], part_box)
		true_dist.append(distances)
		t_solutions_all.append(t_sol)
	true_dist = np.array(true_dist).swapaxes(0,1)
	t_solutions_all = np.array(t_solutions_all).swapaxes(0,1)
	Shortest_dist = np.min(true_dist, axis=1)
	Segment_ids = true_dist.argmin(axis=1)
	# Sort the soltuions of t with respect to the corresponding segment that is shortest to a given particle
	Shortest_t_sol = []
	for i in range(len(Segment_ids)):
		Shortest_t_sol.append(t_solutions_all[i][Segment_ids[i]])
	return Shortest_dist.astype(np.float32), Segment_ids.astype(np.int32), np.array(Shortest_t_sol, np.float32)
    
def ZMQ_get_distances(euclid21check):
	""" Multiprocessing part for the use of ZMQ """
	context = zmq.Context()
	context.linger = 0

	# Socket to receive data from, i.e. the ventilator
	receiver = context.socket(zmq.PULL)
	receiver.connect("tcp://euclid21.ib.intern:5070")
	#receiver.RCVTIMEO = 1000000
    
	# Socket to send computed data to, i.e. back to ventilator
	sender = context.socket(zmq.PUSH)
	sender.connect("tcp://euclid21.ib.intern:5072")

	# Socket controller, ensures the worker is killed
	controller = context.socket(zmq.PULL)
	controller.connect("tcp://euclid21.ib.intern:5074")

	# Only poller for receiver as receiver 2 has similar size
	poller = zmq.Poller()
	poller.register(receiver, zmq.POLLIN)
	poller.register(controller, zmq.POLLIN)
	while True:
		socks = dict(poller.poll(500000))
		""" # Numerical brute froce method
		# Computes context when data is recieved
		if socks.get(receiver) == zmq.POLLIN:
			data = ZMQAS.recv_zipped_pickle(receiver)
			FilamentPos = data[0]
			ParticleBox = data[1]
			ID = data[2]
			BoxSize = data[3]
			Distances, Segpoint_index = Compute_distance(FilamentPos, ParticleBox, BoxSize)
			ZMQAS.send_zipped_pickle(sender, [Distances, Segpoint_index, ID])
		# Closes the context when data computing is done
		if socks.get(controller) == zmq.POLLIN:
			control_message = controller.recv()
			if control_message == "FINISHED":
				break
		if not socks:
			break
		"""
		# Analytical method
		# Computes context when data is recieved
		if socks.get(receiver) == zmq.POLLIN:
			data = ZMQAS.recv_zipped_pickle(receiver)
			FilamentPos = data[0]
			ParticleBox = data[1]
			ID = data[2]
			Distances, Segment_id, t_solution = Compute_distance_analytic(FilamentPos, ParticleBox)
			ZMQAS.send_zipped_pickle(sender, [Distances, Segment_id, ID, t_solution])
		# Closes the context when data computing is done
		if socks.get(controller) == zmq.POLLIN:
			control_message = controller.recv()
			if control_message == "FINISHED":
				break
		if not socks:
			break

	# Finished, closing context
	receiver.close()
	sender.close()
	controller.close()
	context.term()
	
if __name__ == '__main__':
	# This is only called when the script is called from a command line directly.
	# Will then run ZeroMQ paralellziation function.
	# Importing the module will not call this.
	parsed_arguments = Argument_parser()
	Run_euclid21_workers = parsed_arguments.Euclid21Run
	#box_expand = parsed_arguments.BOXEXPAND
	ZMQ_get_distances(Run_euclid21_workers)