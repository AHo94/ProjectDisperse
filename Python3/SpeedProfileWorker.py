# Common modules
#import numpy as np

# ZeroMQ
import zmq
import ZMQArraySending as ZMQAS

# Own functions
import OtherFunctions as OF


def Velocity_bins():
	context = zmq.Context()
	context.linger = 0

	# Socket to receive data from, i.e. the ventilator
	receiver = context.socket(zmq.PULL)
	receiver.connect("tcp://euclid21.ib.intern:4070")
	#receiver.RCVTIMEO = 1000000

	# Socket to send computed data to, i.e. back to ventilator
	sender = context.socket(zmq.PUSH)
	sender.connect("tcp://euclid21.ib.intern:4072")

	# Socket controller, ensures the worker is killed
	controller = context.socket(zmq.PULL)
	controller.connect("tcp://euclid21.ib.intern:4074")
	poller = zmq.Poller()
	poller.register(receiver, zmq.POLLIN)
	poller.register(controller, zmq.POLLIN)
	while True:
		socks = dict(poller.poll(500000))
		# Analytical method
		# Computes context when data is recieved
		if socks.get(receiver) == zmq.POLLIN:
			data = ZMQAS.recv_zipped_pickle(receiver)
			Parts_included = data[0]
			Speeds_included = data[1]
			bins = data[2]
			bin_value, bin_std = OF.Bin_mean_common(Parts_included, Speeds_included, bins)
			ZMQAS.send_zipped_pickle(sender, [bin_value])
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
	Velocity_bins()
