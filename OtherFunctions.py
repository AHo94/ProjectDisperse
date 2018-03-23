# Ohter functions, used for JupyterNotebook testing
# Will be applied mostly for the particle property analysis
import numpy as np
import matplotlib.pyplot as plt


def compute_avg_bins(bins_all, bin_values, binnum=30):
	""" 
	Computes the average bins of a given set of bin values.
	The bin range will be the minimum to maximum between all the bin values given.
	The average will be based on the bin values and number of samples.
	Nan values will be ignored, and the sampling size will be reduced accordingly.
	"""
	bins_average = np.linspace(np.min(bins_all), np.max(bins_all), binnum)
	idx_bin_avg = np.digitize(bins_all, bins_average)
	avg_len_bins = []
	for j in range(len(idx_bin_avg)):
		avg_len_bins.append(np.array([np.average(bin_values[j, i == idx_bin_avg[j]]) for i in range(len(bins_average))]))
	avg_len_bins = np.asarray(avg_len_bins)
	average_bin_value = []
	std = [np.std(bin_values[:,i]) for i in range(len(bins_average))]
	
	for i in range(len(bins_average)):
		Numbers = avg_len_bins[:,i]
		samples = len(Numbers)
		Total = 0.0
		for data in Numbers:
			if np.isnan(data):
				samples -= 1
			else:
				Total += data
		if Total:
			average_bin_value.append(Total/float(samples))
		else:
			average_bin_value.append(np.nan)
	
	return np.array(average_bin_value), bins_average, np.array(std)

def Compute_speed(vx, vy, vz):
	""" Computes the speed of a particle, given velocity components """
	velocity_3d = Get_3D_vel(vx, vy, vz)
	Speeds = np.linalg.norm(velocity_3d, axis=1)
	return np.array(Speeds)

def Get_3D(px, py, pz, vx, vy, vz):
	""" Returns 3D velocity AND particle position coordinates"""
	position_3d = np.column_stack((px, py, pz))
	velocity_3d = np.column_stack((vx, vy, vz))
	return position_3d, velocity_3d

def Get_3D_vel(vx, vy, vz):
	""" Returns 3D velocity coordinates """
	velocity_3d = np.column_stack((vx, vy, vz))
	return velocity_3d

def Get_3D_partpos(px, py, pz):
	""" Returns 3D particle position coordinates """
	position_3d = np.column_stack((px, py, pz))
	return position_3d

def Get_filament_length(FilamentPos):
	""" Returns length of a filament. Periodic boundary are taken into account """
	boxsize = 256.0
	bbb = boxsize/2.0
	FilPos_diff = FilamentPos[1:] - FilamentPos[:-1]
	FilPos_diff[FilPos_diff[:,1] >= bbb,1] -= boxsize
	FilPos_diff[FilPos_diff[:,1] <= -bbb,1] += boxsize
	FilPos_diff[FilPos_diff[:,0] >= bbb,0] -= boxsize
	FilPos_diff[FilPos_diff[:,0] <= -bbb,0] += boxsize
	FilPos_diff[FilPos_diff[:,2] >= bbb,2] -= boxsize
	FilPos_diff[FilPos_diff[:,2] <= -bbb,2] += boxsize

	seglens = np.linalg.norm(FilPos_diff, axis=1)
	fillength = np.sum(seglens)
	return fillength

def Get_segaxis(FilamentPos, SegPoints):
	""" 
	Receives segment axis, of a filament, based on the pre-determined segment points.
	Segment point based on the percentage of the segment length vs the total filament length.
	"""
	Nsegs = len(FilamentPos)-1
	xlims = np.array([])
	ylims = np.array([])
	zlims = np.array([])
	for i in range(Nsegs):
		xlims = np.concatenate((xlims, np.linspace(FilamentPos[i][0], FilamentPos[i+1][0], SegPoints[i])))
		ylims = np.concatenate((ylims, np.linspace(FilamentPos[i][1], FilamentPos[i+1][1], SegPoints[i])))
		zlims = np.concatenate((zlims, np.linspace(FilamentPos[i][2], FilamentPos[i+1][2], SegPoints[i])))
	Segment_axis = []
	for j in range(len(xlims)):
		Segment_axis.append(np.array([xlims[j], ylims[j], zlims[j]]))
	Segment_axis = np.asarray(Segment_axis)
	return Segment_axis

def Histogram_average(Bins, Bin_values):
	""" 
	Computes the average of a set number of histogram, for the same bin range.
	This computes the average for a given set of data at a given bin
	Use this to get average of multiple filament velocities etc.
	"""
	Average_binned_values = []
	for i in range(len(Bins)):
		Numbers = Bin_values[:,i]
		samples = len(Numbers)
		Total = 0.0
		for data in Numbers:
			if np.isnan(data):
				samples -= 1
			else:
				Total += data
		if Total:
			Average_binned_values.append(Total/float(samples))
		else:
			Average_binned_values.append(np.nan)
	return np.array(Average_binned_values)

def Bin_mean(X_value, Y_value, binnum=30):
	""" Computes the mean bin value of a given Y-value vs the X-value, e.g distance from filament vs particle speed """
	#hist, bins = np.histogram(X_value, bins=binnum)
	bins = np.linspace(np.min(X_value), np.max(X_value), binnum)
	index_bin = np.digitize(X_value, bins)
	binval = np.array([ np.mean(Y_value[i == index_bin]) for i in range(len(bins))])
	binstd = [ np.std(Y_value[i == index_bin]) for i in range(len(bins)) ]
	return bins, binval, binstd

def Bin_numbers(X_value, Y_value, binnum=30):
	""" Computes the number of data sets of a given Y-value vs the X-value. E.g number of filaments for a given filament length. """
	#hist, bins = np.histogram(X_value, bins=binnum)
	bins = np.linspace(np.min(X_value), np.max(X_value), binnum)
	index_bin = np.digitize(X_value, bins)
	binval = np.array([ len(Y_value[i == index_bin]) for i in range(len(bins))])
	binstd = [ np.std(Y_value[i == index_bin]) for i in range(len(bins)) ]
	return bins, binval, binstd

def Bin_numbers_logX(X_value, Y_value, binnum=30):
	""" Computes the number of data sets of a given Y-value vs the X-value. E.g number of filaments for a given filament length. """
	#hist, bins = np.histogram(X_value, bins=binnum)
	bins = np.exp(np.linspace(np.log(np.min(X_value)), np.log(np.max(X_value)), binnum))
	index_bin = np.digitize(X_value, bins)
	binval = np.array([ len(Y_value[i == index_bin]) for i in range(len(bins))])
	binstd = [ np.std(Y_value[i == index_bin]) for i in range(len(bins)) ]
	return bins, binval, binstd

def Get_common_bin(Data, binnum=30):
	""" Gets a common binning based on the dataset given. Assumes multiple dataset in one variable, e.g [[data1], [data2], ...] """
	Minima = np.array([np.min(Data[i]) for i in range(len(Data))])
	Maxima = np.array([np.max(Data[i]) for i in range(len(Data))])
	Max_value = np.max(Maxima)
	Min_value = np.min(Minima)
	Bin = np.linspace(Min_value, Max_value, binnum)
	return Bin

def Get_common_bin_logX(Data, binnum=30):
	""" Gets a common binning based on the dataset given. Assumes multiple dataset in one variable, e.g [[data1], [data2], ...] """
	Minima = np.array([np.min(Data[i]) for i in range(len(Data))])
	Maxima = np.array([np.max(Data[i]) for i in range(len(Data))])
	Max_value = np.max(Maxima)
	Min_value = np.min(Minima)
	Bin = np.exp(np.linspace(np.log(Min_value), np.log(Max_value), binnum))
	return Bin

def Bin_numbers_common(X_value, Y_value, bins, std='std'):
	""" 
	Computes number of data set of a given Y-value vs the X-value. Uses common binning. 
	Two different errors to be chosen:
	1) Standard deviation
	2) Poisson distribution, where error = sqrt(N), N = number of data in a given bin.
	"""
	index_bin = np.digitize(X_value, bins)
	binval = np.array([float(len(Y_value[i == index_bin])) for i in range(len(bins))])
	if std == 'std':
		binstd = np.array([ np.std(Y_value[i == index_bin]) for i in range(len(bins)) ])
	elif std == 'poisson':
		binstd = np.array([np.sqrt(numbers) for numbers in binval])
	else:
		raise ValueError("Argument std not set properly! Try std='std' or std='poisson'")
	return binval, binstd

def Bin_mean_common(X_value, Y_value, bins):
	index_bin = np.digitize(X_value, bins)
	binval = np.array([ np.mean(Y_value[i == index_bin]) for i in range(len(bins))])
	binstd = np.array([ np.std(Y_value[i == index_bin]) for i in range(len(bins)) ])
	Number_points = np.array([len(Y_value[i == index_bin]) for i in range(len(bins))]).astype(np.float32)
	return binval, binstd/np.sqrt(Number_points) 		# FIXME? Divide by N or sqrt(N)?

def Mean_per_bin(bins, bin_values):
	Npts = len(bins)
	std = np.array([np.std(bin_values[:,i]) for i in range(Npts)])
	Mean = np.array([np.mean(bin_values[:,i]) for i in range(Npts)])
	Number_points = np.array([len(bin_values[:,i]) for i in range(Npts)])
	return Mean, std/Number_points


def filament_box(filament):
	""" Gives filament box, plotting purposes """
	xmin = np.min(filament[:,0])
	xmax = np.max(filament[:,0])
	ymin = np.min(filament[:,1])
	ymax = np.max(filament[:,1])
	zmin = np.min(filament[:,2])
	zmax = np.max(filament[:,2])
	return xmin, xmax, ymin, ymax, zmin, zmax

def particle_mask(box):
	""" Masks particles, plotting purposes """
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

def particle_box_plot(filament_box, dist):
	""" Masks the particle box around filament, plotting purposes """
	xmin = filament_box[0] - dist
	xmax = filament_box[1] + dist
	ymin = filament_box[2] - dist
	ymax = filament_box[3] + dist
	zmin = filament_box[4] - dist
	zmax = filament_box[5] + dist
	box = [xmin, xmax, ymin, ymax, zmin, zmax]
	box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
	MovePartx = 0
	MoveParty = 0
	MovePartz = 0
	Atboundary = 0
	if np.abs((xmin + dist) - lower_boundary) < 1e-3:
		box[0] = xmin + dist
		Atboundary = 1
	elif np.abs((xmax - dist) - lower_boundary) < 1e-3:
		box[1] = xmax - dist
		Atboundary = 1
	if np.abs((ymin + dist) - lower_boundary) < 1e-3:
		box[2] = ymin + dist
		Atboundary = 1
	elif np.abs((ymax - dist) - lower_boundary) < 1e-3:
		box[3] = ymax - dist
		Atboundary = 1
	if np.abs((zmin + dist) - lower_boundary) < 1e-3:
		box[4] = zmin + dist
		Atboundary = 1
	elif np.abs((zmax - dist) - lower_boundary) < 1e-3:
		box[5] = zmax - dist
		Atboundary = 1
	#Atboundary = 0
	if Atboundary:
		box2 = box
	else:
		if xmin < lower_boundary and not Atboundary:
			box[0] = lower_boundary
			box2[1] = upper_boundary
			box2[0] = upper_boundary + (xmin - lower_boundary)
			MovePartx = -256.0
		elif xmax > upper_boundary and not Atboundary:
			box[1] = upper_boundary
			box2[0] = lower_boundary
			box2[1] = lower_boundary + (xmax - upper_boundary)
			MovePartx = 256.0

		if ymin < lower_boundary and not Atboundary:
			box[2] = lower_boundary
			box2[3] = upper_boundary
			box2[2] = upper_boundary + (ymin - lower_boundary)
			MoveParty = -256.0
		elif ymax > upper_boundary and not Atboundary:
			box[3] = upper_boundary
			box2[2] = lower_boundary
			box2[3] = lower_boundary + (ymax - upper_boundary)
			MoveParty = 256.0

		if zmin < lower_boundary and not Atboundary:
			box[4] = lower_boundary
			box2[5] = upper_boundary
			box2[4] = upper_boundary + (zmin - lower_boundary)
			MovePartz = -256.0
		elif zmax > upper_boundary and not Atboundary:
			box[5] = upper_boundary
			box2[4] = lower_boundary
			box2[5] = lower_boundary + (zmax - upper_boundary)
			MovePartz = 256.0
	
	# Still need to move particles to the other side of the boundary
	if box == box2:
		mask = particle_mask(box)
		return ParticlePos[mask]
	else:
		mask1 = particle_mask(box)
		mask2 = particle_mask(box2)
		Particles1 = ParticlePos[mask1]
		Particles2 = ParticlePos[mask2]
		Particles2[:,0] += MovePartx
		Particles2[:,1] += MoveParty
		Particles2[:,2] += MovePartz
		Masked_particles = np.concatenate([Particles1, Particles2])
		return Masked_particles

def find_middle_point(filament, fillen):
	""" 
	Finds the middle point of a filament 
	First find which segment is the middle segment. Done by considering segment length and comparing that to total filament length.
	If added segment lengths are greater than 0.5*filament_length, then we have found the middle segment
	Parametrize s(t) = s_i + t*(s_{i+1} + s_i), find a t_value such that |s(t)| + segment_lengths = 0.5*filament_length
	"""
	seglens = np.linalg.norm(filament[1:] - filament[:-1], axis=1)
	currlen = 0
	segnum = 0
	for i in range(len(seglens)):
		currlen += seglens[i]
		if currlen >= fillen/2.0:
			segnum = i
			break
		else:
			continue
	Difference = filament[segnum+1] - filament[segnum]
	diffx = Difference[0]
	diffy = Difference[1]
	diffz = Difference[2]
	diffx += 256.0*(diffx <= -256.0/2.0)
	diffx -= 256.0*(diffx >= 256.0/2.0)
	diffy += 256.0*(diffy <= -256.0/2.0)
	diffy -= 256.0*(diffy >= 256.0/2.0)
	diffz += 256.0*(diffz <= -256.0/2.0)
	diffz -= 256.0*(diffz >= 256.0/2.0)
	
	a = diffx*diffx + diffy*diffy + diffz*diffz
	c =  -(fillen*fillen)/4.0
	t_sol = np.roots([a,0,c])
	real_t1 = t_sol >= 0
	
	t = t_sol[real_t1]	# Only select positive valued t
	if t > 1.0:
		t = 1
	midpoint = filament[segnum] + t*(filament[segnum+1] - filament[segnum])
	return midpoint

def get_filament_distances(midpts):
	""" Computes filament distances from the midpoints """
	Filament_distances = []
	for i in range(len(midpts)):
		distances = np.zeros(len(midpts)-1)
		if i > 0:
			for k in range(i):
				distances[k] = Filament_distances[k][-1] if (i == (len(midpts)-1)) else Filament_distances[k][i-1]
		midpts_diff = midpts[i+1:] - midpts[i]
		midpts_diff[midpts_diff >= 256.0/2.0] -= 256.0
		midpts_diff[midpts_diff <= -256.0/2.0] += 256.0
		distances[i:] = np.linalg.norm(midpts_diff, axis=1)
		Filament_distances.append(distances)
	return np.array(Filament_distances)

def Periodic_slicing(limit_crossing, boundary, Yslices):
	""" Used to determine which particle slice box the filament requires. """
	min_crossing = limit_crossing[0]
	max_crossing = limit_crossing[1]
	if boundary == 'upper':
		Max_cross = max_crossing - 256.0
		Min_cross = 0.0
		max_crossing = 256.0
		Mask_t = (Yslices[:,0] <= max_crossing) & (Yslices[:,1] >= min_crossing)
		Mask_b = (Yslices[:,0] <= Max_cross) & (Yslices[:,1] >= Min_cross)
		indices1 = np.where(Mask_t)[0]
		indices2 = np.where(Mask_b)[0]
		idx_cross = np.concatenate((indices1, indices2))
	elif boundary == 'lower':
		Min_cross = min_crossing + 256.0
		Max_cross = 256.0
		min_crossing = 0.0
		Mask_t = (Yslices[:,0] <= max_crossing) & (Yslices[:,1] >= min_crossing)
		Mask_b = (Yslices[:,0] <= Max_cross) & (Yslices[:,1] >= Min_cross)
		indices1 = np.where(Mask_t)[0]
		indices2 = np.where(Mask_b)[0]
		idx_cross = np.concatenate((indices1, indices2))
	elif boundary == 'none':
		Mask_t = (Yslices[:,0] <= max_crossing) & (Yslices[:,1] >= min_crossing)
		idx_cross = np.where(Mask_t)[0]
	else:	   
		raise ValueError("Boundary not set properly!")
	return idx_cross


def get_indices_slicing(Filament_slice, Yslices, box_expand):
	""" Gets indices of the slice array. Used to determine which slice box is sent for masking. """
	Ymax = np.max(Filament_slice[:,1]) + box_expand
	Ymin = np.min(Filament_slice[:,1]) - box_expand
	Zmax = np.max(Filament_slice[:,2]) + box_expand
	Zmin = np.min(Filament_slice[:,2]) - box_expand

	if Ymax > 256.0:
		idY = Periodic_slicing([Ymin, Ymax], 'upper', Yslices)
	elif Ymin < 0.0:
		idY = Periodic_slicing([Ymin, Ymax], 'lower', Yslices)
	else:
		idY = Periodic_slicing([Ymin, Ymax], 'none', Yslices)
	if Zmax > 256.0:
		idZ = Periodic_slicing([Zmin, Zmax], 'upper', Yslices)
	elif Zmin < 0.0:
		idZ = Periodic_slicing([Zmin, Zmax], 'lower', Yslices)
	else:
		idZ = Periodic_slicing([Zmin, Zmax], 'none', Yslices)

	return idY, idZ

def Get_particle_box_slices(ParticlePos):
	""" Creates slices of the particle box, for a speed up in masking. """
	Nslices = 20
	SliceSize = 256.0/Nslices
	Yparts = []
	IDparts = []
	PartIDs = np.linspace(0, len(ParticlePos)-1, len(ParticlePos), dtype=np.int32)
	for i in range(Nslices):
		MM = (ParticlePos[:,1] > SliceSize*i) & (ParticlePos[:,1] < SliceSize*(i+1))
		Yparts.append(ParticlePos[MM])
		IDparts.append(PartIDs[MM])
	SlicedParts = []
	SlicedIDs = []
	for i in range(Nslices):
		TempSlices = []
		TempIDSlices = []
		for j in range(Nslices):
			Ybox = Yparts[i]
			MM = (Ybox[:,2] > SliceSize*j) & (Ybox[:,2] < SliceSize*(j+1))
			TempSlices.append(Yparts[i][MM])
			TempIDSlices.append(IDparts[i][MM])
		SlicedParts.append(np.array(TempSlices))
		SlicedIDs.append(np.array(TempIDSlices))
	SlicedParts = np.asarray(SlicedParts)
	SlicedIDs = np.asarray(SlicedIDs)
	# Ranges of each slice
	SliceRanges = []
	for i in range(Nslices):
		SliceRanges.append([i*SliceSize, (i+1)*SliceSize])
	SliceRanges = np.asarray(SliceRanges)
	return SlicedParts, SlicedIDs, SliceRanges

def relative_deviation(data, index):
	""" 
	Computes relative deviation of a given physical quantity.
	This function assumes that the base model is in the first element.
	"""
	delta = (data[index] - data[0])/data[0]
	return delta

def Propagate_error_reldiff(data_0, data, error_0, error):
	""" 
	Propagated error of the relative difference F = (B - A)/A, A = base model
	This function assumes the derivatives is with the data itself
	Derivatives
	dFdA = -B/A**2
	dfdB = 1/A
	"""
	dfdDi = 1/data_0
	dfdD0 = -data/data_0**2
	Prop_error = np.sqrt((error_0**2)*(dfdD0**2) + (error**2)*(dfdDi)**2)
	return Prop_error