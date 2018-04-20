# Common modules
import numpy as np
import time

def Read_halo_data(model):
	""" Reads halo data from  ahf catalogue"""
	toggle = False
	Model_check = ['lcdm', 'symm_A', 'symm_B', 'symm_C', 'symm_D', 'fofr4', 'fofr5', 'fofr6']
	for models in Model_check:
		if model == models:
			toggle = True
	if not toggle:
		raise ValueError('Model input name %s not correctly set into the halo file reader.' %modelfile)

	print "Reading halo catalogue file..."
	if model == 'lcdm':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/lcdm/z0.000/ahf_halocatalog/output.z0.005.AHF_halos'
	elif model == 'symm_A':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/symm_A/z0.000/ahf_halocatalog/output.z0.004.AHF_halos'	
	elif model == 'symm_B':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/symm_B/z0.000/ahf_halocatalog/output.z0.005.AHF_halos'	
	elif model == 'symm_C':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/symm_C/z0.000/ahf_halocatalog/output.0009.z0.003.AHF_halos'	
	elif model == 'symm_D':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/symm_D/z0.000/ahf_halocatalog/output.z0.003.AHF_halos'	
	elif model == 'fofr4':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/fofr4/z0.000/ahf_halocatalog/output.z0.002.AHF_halos'
	elif model == 'fofr5':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/fofr5/z0.000/ahf_halocatalog/output.z0.001.AHF_halos'
	elif model == 'fofr6':
		filedir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/fofr6/z0.000/ahf_halocatalog/output.0007.z0.004.AHF_halos'
	HaloID = []
	Position = []
	Radius_vir = []
	HaloMass = []
	timer = time.time()
	datafiles = open(filedir, 'r')
	skipfirst = 0
	for line in datafiles:
		data_sets = line.split()
		if skipfirst:
			data_sets = line.split()
			HaloID.append(data_sets[0])
			HaloMass.append(data_sets[3])
			Position.append(data_sets[5:8])
			Radius_vir.append(data_sets[11])
		else:
			skipfirst = 1
	print "Time took reading halo file:", time.time() - timer, "s"
	HaloID = np.asarray(HaloID).astype(np.int32)
	HaloMass = np.asarray(HaloMass).astype(np.float32)
	Position = np.asarray(Position).astype(np.float32)/100.0 	# Convert Kpc/h to Mpc/h
	Radius_vir = np.asarray(Radius_vir).astype(np.float32)/100.0 	#Conver to Kpc/h to Mpc/h
	return HaloID, Position, HaloMass, Radius_vir
