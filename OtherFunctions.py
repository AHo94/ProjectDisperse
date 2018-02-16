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
	velocity_3d = []
	for i in range(len(vx)):
		velocity_3d.append(np.array([vx[i], vy[i], vz[i]]))
	Speeds = []
	for vel in velocity_3d:
		Speeds.append(np.linalg.norm(vel))
	return np.array(Speeds)

def Get_3D(px, py, pz, vx, vy, vz):
	""" Returns 3D velocity AND particle position coordinates"""
	velocity_3d = []
	position_3d = []
	for i in range(len(vx)):
		velocity_3d.append(np.array([vx[i], vy[i], vz[i]]))
		position_3d.append(np.array([px[i], py[i], pz[i]]))
	return np.asarray(position_3d), np.asarray(velocity_3d)

def Get_3D_vel(vx, vy, vz):
	""" Returns 3D velocity coordinates """
	velocity_3d = []
	for i in range(len(vx)):
		velocity_3d.append(np.array([vx[i], vy[i], vz[i]]))
	return np.asarray(velocity_3d)

def Get_3D_partpos(px, py, pz):
	""" Returns 3D particle position coordinates """
	position_3d = []
	for i in range(len(px)):
		position_3d.append(np.array([px[i], py[i], pz[i]]))
	return np.asarray(position_3d)

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

def Average_hist(Bins, Bin_values):
	""" Alternative average histogram computer """
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