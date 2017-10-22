import numpy as np
from scipy import stats
import time

class InterpolateDensity():
	"""
	This class is used to interpolate, using scipy's guassian_kde function, the positions of the dark matter particles.
	This is used to find the density profile of the dark matter particles, and compares the filamentary structure with disperse.
	One may also compute a zoomed in section of the box, given a zoom in limit.
	"""
	def __init__(self, bandwidth, particlepositions, grid, npoints=100j):
		self.PartPosX = particlepositions[0]
		self.PartPosY = particlepositions[1]
		self.PartPosZ = particlepositions[2]
		self.xmin = grid[0]
		self.xmax = grid[1]
		self.ymin = grid[2]
		self.ymax = grid[3]
		self.npoints = npoints
		if bandwidth != 'Scott':
			self.method = float(bandwidth)
		else:
			self.method = bandwidth

	def Compute_density(self, particle_pos, grid_pos, x_array):
		kernel = stats.gaussian_kde(particle_pos, bw_method=self.method)
		Density = np.reshape(kernel(grid_pos).T, x_array.shape)
		Log_density = np.log(Density/np.average(Density))
		return Density, Log_density

	def Interpolate_DM_particles(self):
		""" Interpolates dark matter particles within the whole grid """
		print 'Interpolating dark matter particle density'
		time_start = time.clock()
		X, Y = np.mgrid[self.xmin:self.xmax:self.npoints, self.ymin:self.ymax:self.npoints]
		positions = np.vstack([X.ravel(), Y.ravel()])
		particle_positions = np.vstack([self.PartPosX, self.PartPosY])
		Interpolated_Z, Logarithmic_density = self.Compute_density(particle_positions, positions, X)
		print 'Interplation time:', time.clock() - time_start, 's'
		return Interpolated_Z, Logarithmic_density

	def compute_zoomed_density(self, Zoom_areas):
		""" 
		Interpolates the dark matter particles within a given zoomed in section
		Input must be an array or list, which contains the boundaries of x and y direction
		For example: Zoom_areas = [xmin, xmax, ymin, ymax]
		Multiple zoom ins: Zoom_areas = [[xmin1, xmax1, ymin1, ymax1], [xmin2, xmax2, ymin2, ymax2]]
		"""
		print 'Interpolating dark matter density in a zoomed in section'
		def Density_zoomed(ZoomArea):
			""" Calculates density profile of the zoomed in section. Assumes UnitConverter is on by default """
			Xzoom, Yzoom = np.mgrid[ZoomArea[0]:ZoomArea[1]:self.npoints, ZoomArea[2]:ZoomArea[3]:self.npoints]
			position_zoom = np.vstack([Xzoom.ravel(), Yzoom.ravel()])
			Particle_positions = np.column_stack([self.PartPosX, self.PartPosY, self.PartPosZ])
			mask = np.logical_and(
				np.logical_and(Particle_positions[:,0] > ZoomArea[0], Particle_positions[:,0] < ZoomArea[1]),
				np.logical_and(Particle_positions[:,1] > ZoomArea[2], Particle_positions[:,1] < ZoomArea[3])
				)
			Masked_positions = Particle_positions[mask]
			MaskedX = Masked_positions[:,0]
			MaskedY = Masked_positions[:,1]
			partpositions = np.vstack([MaskedX, MaskedY])
			Density, Log_density = self.Compute_density(partpositions, position_zoom, Xzoom)
			return Density, Log_density

		Zoomed_density_list = []
		Log_zoomed_density_list = []	
		for zoom_grid in Zoom_areas:
			Zoomed_density, Log_zoomed_density = Density_zoomed(zoom_grid)
			Zoomed_density_list.append(Zoomed_density)
			Log_zoomed_density_list.append(Log_zoomed_density)
		return Zoomed_density_list, Log_zoomed_density_list