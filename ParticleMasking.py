# Basic modules
import numpy as np

# ZMQ modules
import zmq
import ZMQArraySending as ZMQAS


def filament_box(filament):
	""" Creates a box around the filament """
	xmin = np.min(filament[:,0])
	xmax = np.max(filament[:,0])
	ymin = np.min(filament[:,1])
	ymax = np.max(filament[:,1])
	zmin = np.min(filament[:,2])
	zmax = np.max(filament[:,2])
	return xmin, xmax, ymin, ymax, zmin, zmax

def particle_mask(box, ParticlePos_):
	""" Masks particles based on the particle box received """
	xmin = box[0]
	xmax = box[1]
	ymin = box[2]
	ymax = box[3]
	zmin = box[4]
	zmax = box[5]
	maskx = (ParticlePos_[:,0] > xmin) & (ParticlePos_[:,0] < xmax)
	masky = (ParticlePos_[:,1] > ymin) & (ParticlePos_[:,1] < ymax)
	maskz = (ParticlePos_[:,2] > zmin) & (ParticlePos_[:,2] < zmax)
	return maskx*masky*maskz

def masked_particle_indices(filament_box, ParticlePos_):
	""" 
	Masks particle indices based on the filament box and particle positions.
	The indices may be used to identify velocities etc.
	"""
	xmin = filament_box[0] - 3.0
	xmax = filament_box[1] + 3.0
	ymin = filament_box[2] - 3.0
	ymax = filament_box[3] + 3.0
	zmin = filament_box[4] - 3.0
	zmax = filament_box[5] + 3.0
	box = [xmin, xmax, ymin, ymax, zmin, zmax]
	box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
	MovePartx = 0
	MoveParty = 0
	MovePartz = 0
	Atboundary = 0
	# Check which boundary it crosses
	if np.abs((xmin + 3.0) - 0.0) < 1e-3:
		box[0] = xmin + 3.0
		Atboundary = 1
	elif np.abs((xmax - 3.0) - 256.0) < 1e-3:
		box[1] = xmax - 3.0
		Atboundary = 1
	if np.abs((ymin + 3.0) - 0.0) < 1e-3:
		box[2] = ymin + 3.0
		Atboundary = 1
	elif np.abs((ymax - 3.0) - 256.0) < 1e-3:
		box[3] = ymax - 3.0
		Atboundary = 1
	if np.abs((zmin + 3.0) - 0.0) < 1e-3:
		box[4] = zmin + 3.0
		Atboundary = 1
	elif np.abs((zmax - 3.0) - 256.0) < 1e-3:
		box[5] = zmax - 3.0
		Atboundary = 1
		
	# If boundary is crossed, create a new box
	if Atboundary:
		box2 = box
	else:
		if xmin < 0.0:
			box[0] = 0.0
			box2[1] = 256.0
			box2[0] = 256.0 + (xmin - 0.0)
		elif xmax > 256.0:
			box[1] = 256.0
			box2[0] = 0.0
			box2[1] = 0.0 + (xmax - 256.0)
			
		if ymin < 0.0:
			box[2] = 0.0
			box2[3] = 256.0
			box2[2] = 256.0 + (ymin - 0.0)
		elif ymax > 256.0:
			box[3] = 256.0
			box2[2] = 0.0
			box2[3] = 0.0 + (ymax - 256.0)
			
		if zmin < 0.0:
			box[4] = 0.0
			box2[5] = 256.0
			box2[4] = 256.0 + (zmin - 0.0)
		elif zmax > 256.0:
			box[5] = 256.0
			box2[4] = 0.0
			box2[5] = 0.0 + (zmax - 256.0)

	if box == box2:
		mask = particle_mask(box, ParticlePos_)
		return (np.where(mask)[0]).astype(np.int32)
	else:
		mask1 = particle_mask(box, ParticlePos_)
		mask2 = particle_mask(box2, ParticlePos_)
		indices1 = (np.where(mask1)[0]).astype(np.int32)
		indices2 = (np.where(mask2)[0]).astype(np.int32)
		Masked_indices = np.concatenate([indices1, indices2])
		return Masked_indices

def particle_box(filament_box, masked_ids, ParticlePos_):
	""" 
	Creates a particle box and masks off the particles within the filament box.
	Returns particle positions and not IDs.
	"""
	xmin = filament_box[0] - 3.0
	xmax = filament_box[1] + 3.0
	ymin = filament_box[2] - 3.0
	ymax = filament_box[3] + 3.0
	zmin = filament_box[4] - 3.0
	zmax = filament_box[5] + 3.0
	MovePartx = 0
	MoveParty = 0
	MovePartz = 0
	# Check which boundary is crossed. If crossed, moves the other side
	if xmin < 0.0:
		MovePartx = -256.0
	elif xmax > 256.0:
		MovePartx = 256.0
			
	if ymin < 0.0:
		MoveParty = -256.0
	elif ymax > 256.0:
		MoveParty = 256.0
			
	if zmin < 0.0:
		MovePartz = -256.0
	elif zmax > 256.0:
		MovePartz = 256.0

	if not len(masked_ids) == 2:
		return ParticlePos_[masked_ids]
	else:
		Particles1 = ParticlePos_[masked_ids[0]]
		Particles2 = ParticlePos_[masked_ids[1]]
		Particles2[:,0] += MovePartx
		Particles2[:,1] += MoveParty
		Particles2[:,2] += MovePartz
		Masked_particle_box = np.concatenate([Particles1, Particles2])
		return Masked_particle_box


def ZMQ_mask_particles():
	""" Uses ZMQ to mask the particles. An alternative to the multiprocessing method. """
	context = zmq.Context()
	context.linger = 0

	# Socket to receive data from
	receiver = context.socket(zmq.PULL)
	receiver.connect("tcp://euclid21.ib.intern:6070")
	
	# Socket to send computed data to
	sender = context.socket(zmq.PUSH)
	sender.connect("tcp://euclid21.ib.intern:6072")

	# Socket controller, ensures the worker is killed
	controller = context.socket(zmq.PULL)
	controller.connect("tcp://euclid21.ib.intern:6074")

	# Only poller for receiver as receiver 2 has similar size
	poller = zmq.Poller()
	poller.register(receiver, zmq.POLLIN)
	poller.register(controller, zmq.POLLIN)
	while True:
		if socks.get(receiver) == zmq.POLLIN:
			data = ZMQAS.recv_zipped_pickle(receiver)
			FilamentPos = data[0]
			ParticleBox = data[1]
			ID = data[2]
			filbox = filament_box(FilamentPos)
			Masked_ids = masked_particle_indices(filbox, ParticleBox)
			Masked_part_box = particle_box(filbox, Masked_ids, ParticleBox)
			ZMQAS.send_zipped_pickle(sender, [Masked_ids, Masked_part_box, ID])
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
	ZMQ_mask_particles()