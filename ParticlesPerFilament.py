import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Global variables
# Box boundaries assumed to be from 0 to 256.0 Mpc/h
lower_boundary = 0.0
upper_boundary = 256.0

def filament_box(filament):
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

def particle_mask(box, ParticlePos):
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

def particle_box(filament_box, ParticlePos, box_expand):
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
	xmin = filament_box[0] - box_expand
	xmax = filament_box[1] + box_expand
	ymin = filament_box[2] - box_expand
	ymax = filament_box[3] + box_expand
	zmin = filament_box[4] - box_expand
	zmax = filament_box[5] + box_expand
	box = [xmin, xmax, ymin, ymax, zmin, zmax]
	box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
	MovePartx = 0
	MoveParty = 0
	MovePartz = 0
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
	# Still need to move particles to the other side of the boundary
	if box == box2:
		mask = particle_mask(box)
		return ParticlePos[mask]
	else:
		mask1 = particle_mask(box, ParticlePos)
		mask2 = particle_mask(box2, ParticlePos)
		Particles1 = ParticlePos[mask1]
		Particles2 = ParticlePos[mask2]
		Particles2[:,0] += MovePartx
		Particles2[:,1] += MoveParty
		Particles2[:,2] += MovePartz
		Masked_particles = np.concatenate([Particles1, Particles2])
		return Masked_particles

def get_distance(filament, part_box):
	""" 
	Computes the distance from a particle to every segment in the filament.
	Only the shortest distance is included.
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

def solve(filament, particlepos, distance_threshold, box_expand):
	filbox = filament_box(filament)
	partbox = particle_box(filbox, particlepos, box_expand)
	distances = get_distance(filament, partbox)
	return len(np.where(distances <= distance_threshold)[0])

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


"""
def PartPerFilament(filaments, filamentID, particlepos):
	fil = filaments[0]
	xmin = np.min(fil[:,0])
	xmax = np.max(fil[:,0])
	ymin = np.min(fil[:,1])
	ymax = np.max(fil[:,1])
	zmin = np.min(fil[:,2])
	zmax = np.max(fil[:,2])

	FilamentPos = np.column_stack((fil[:,0], fil[:,1]))
	#x = fil[:,0]
	#y = fil[:,1]
	#FilamentPos = []
	#for i in range(len(x)):
	#	FilamentPos.append(np.array([x[i], y[i]]))
	#print FilamentPos
	fig, ax = plt.subplots()
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	line = LineCollection([FilamentPos], linestyle='solid')
	ax.add_collection(line)
	plt.show()
"""