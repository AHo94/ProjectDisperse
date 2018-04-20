# Common modules
import numpy as np
import time

def Read_halo_data(model):
	""" Reads halo data from  ahf catalogue"""
	toggle = False
	Model_check = ['lcdm', 'symm_A', 'symm_B', 'symm_C', 'symm_D', 'fofr4', 'fofr5', 'fofr6']
	for models in Model_check:
		if modelfile == models:
			toggle = True
	if not toggle:
		raise ValueError('Model input name %s not correctly set into the halo file reader.' %modelfile)

	print "Reading halo catalogue file..."
	filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/' + model + '/z0.00/ahf_halocatalog/'
	HaloID = []
	Position = []
	Radius_vir = []
	timer = time.time()
	datafiles = open(filedir, 'r')
	skipfirst = 0
	for line in datafiles:
		data_sets = line.split()
		if skipfirst:
			data_sets = line.split()
			HaloID.append(data_sets[0])
			Position.append(data_sets[5:8])
			Radius_vir.append(data_sets[11])
		else:
			skipfirst = 1
	print "Time took reading halo file:", time.time() - timer, "s"
	HaloID = np.asarray(HaloID).astype(np.int32)
	Position = np.asarray(Position).astype(np.float32)/100.0 	# Convert Kpc/h to Mpc/h
	Radius_vir = np.asarray(Radius_vir).astype(np.float32)/100.0 	#Conver to Kpc/h to Mpc/h
	return HaloID, Position, Radius_vir
