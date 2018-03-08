# Basic modules
import numpy as np
import os
import sys
import argparse

# ZMQ modules
import zmq
import ZMQArraySending as ZMQAS

def Argument_parser():
	""" Parses optional argument when program is run from the command line """
	parser = argparse.ArgumentParser()
	# Optional arguments
	parser.add_argument("-ID", "--WorkID", help="Worker ID", type=int, default=0)
	# Parse arguments
	args = parser.parse_args()
	return args

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

def masked_particle_indices(filament_box, ParticlePos_, box_expand):
	""" 
	Masks particle indices based on the filament box and particle positions.
	The indices may be used to identify velocities etc.
	"""
	""" 
	Masks particle indices based on the filament box and particle positions.
	The indices may be used to identify velocities etc.
	"""
	xmin = filament_box[0] - box_expand
	xmax = filament_box[1] + box_expand
	ymin = filament_box[2] - box_expand
	ymax = filament_box[3] + box_expand
	zmin = filament_box[4] - box_expand
	zmax = filament_box[5] + box_expand
	box = [xmin, xmax, ymin, ymax, zmin, zmax]
	box2 = [xmin, xmax, ymin, ymax, zmin, zmax]
	SingleBox = 0
	leftovers = np.array([0.0, 0.0, 0.0])
	# Check if filament crosses beyond the particle box. Limits the first box
	if xmin < 0.0:
		box[0] = 0.0
		leftovers[0] = xmin
	elif xmax > 256.0:
		box[1] = 256.0
		leftovers[0] = xmax - 256.0
	if ymin < 0.0:
		box[2] = 0.0
		leftovers[1] = ymin
	elif ymax > 256.0:
		box[3] = 256.0
		leftovers[1] = ymax - 256.0
	if zmin < 0.0:
		box[4] = 0.0
		leftovers[2] = zmin
	elif zmax > 256.0:
		box[5] = 256.0
		leftovers[2] = zmax - 256.0

	# If boundary is crossed, create a new box that covers the other side of the boundary
	if not leftovers.any():
		SingleBox = 1
	else:
		if xmin < 0.0:
			box2[1] = 256.0
			box2[0] = 256.0 + leftovers[0]
		elif xmax > 256.0:
			box2[0] = 0.0
			box2[1] = leftovers[0]
		if ymin < 0.0:
			box2[3] = 256.0
			box2[2] = 256.0 + leftovers[1]
		elif ymax > 256.0:
			box2[2] = 0.0
			box2[3] = leftovers[1]
		if zmin < 0.0:
			box2[5] = 256.0
			box2[4] = 256.0 + leftovers[2]
		elif zmax > 256.0:
			box2[4] = 0.0
			box2[5] = leftovers[2]
            
	if SingleBox:
		mask = particle_mask(box, ParticlePos_)
		return (np.where(mask)[0]).astype(np.int32)
	else:
		mask1 = particle_mask(box, ParticlePos_)
		mask2 = particle_mask(box2, ParticlePos_)
		indices1 = (np.where(mask1)[0]).astype(np.int32)
		indices2 = (np.where(mask2)[0]).astype(np.int32)
		Masked_indices = np.concatenate([indices1, indices2])
		return Masked_indices


def particle_box(filament_box, masked_ids, ParticlePos_, box_expand):
	""" 
	Creates a particle box and masks off the particles within the filament box.
	Returns particle positions and not IDs.
	"""
	xmin = filament_box[0] - box_expand
	xmax = filament_box[1] + box_expand
	ymin = filament_box[2] - box_expand
	ymax = filament_box[3] + box_expand
	zmin = filament_box[4] - box_expand
	zmax = filament_box[5] + box_expand
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
	receiver.connect("tcp://euclid21.ib.intern:5080")
	
	# Socket to send computed data to
	sender = context.socket(zmq.PUSH)
	sender.connect("tcp://euclid21.ib.intern:5082")

	# Socket controller, ensures the worker is killed
	controller = context.socket(zmq.PULL)
	controller.connect("tcp://euclid21.ib.intern:5084")

	# Only poller for receiver as receiver 2 has similar size
	poller = zmq.Poller()
	poller.register(receiver, zmq.POLLIN)
	poller.register(controller, zmq.POLLIN)
	#print 'WORKER READY', WORKER_ID
	while True:
		socks = dict(poller.poll(100000))
		if socks.get(receiver) == zmq.POLLIN:
			data = ZMQAS.recv_zipped_pickle(receiver)
			FilamentPos = data[0]
			ParticleBox = data[1]
			ID = data[2]
			box_expand = data[3]
			filbox = filament_box(FilamentPos)
			Masked_ids = masked_particle_indices(filbox, ParticleBox, box_expand)
			Masked_part_box = particle_box(filbox, Masked_ids, ParticleBox, box_expand)
			#print 'ID = ', ID, WORKER_ID
			ZMQAS.send_zipped_pickle(sender, [Masked_ids, Masked_part_box, ID])
			#ZMQAS.send_zipped_pickle(sender, [FilamentPos*ParticleBox, ID*box_expand, ID])
		# Closes the context when data computing is done
		if socks.get(controller) == zmq.POLLIN:
			control_message = controller.recv()
			if control_message == "FINISHED":
				break
		if not socks:
			print 'Did not receive shit, closing...'#, WORKER_ID
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
	WORKER_ID = parsed_arguments.WorkID
	ZMQ_mask_particles()