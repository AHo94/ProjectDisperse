# Common modules
import numpy as np
import cPickle as pickle

# Own modules
import ReadHaloData as RHD

class AnalyseCritPts():
	def __init__(self, model, npart, sigma):
		self.npart = npart
		self.sigma = sigma
		if model == 'lcdm':
			Dispersemodel = 'LCDM'
		elif model == 'symmA':
			model = 'symm_A'
			Dispersemodel = 'SymmA'
		elif model == 'symmB':
			model = 'symm_B'
			Dispersemodel = 'SymmB'
		elif model == 'symmC':
			model = 'symm_C'
			Dispersemodel = 'SymmC'
		elif model == 'symmD':
			model = 'symm_D'
			Dispersemodel = 'SymmD'
		else:
			Dispersemodel = model

		self.read_filament_data(Dispersemodel)
		self.HaloID, self.HaloPos, self.VirRadius = RHD.Read_halo_data(model)

	def read_filament_data(self, model, npart, sigma):
		""" Reads filament data from DisPerSE """
		cachedir_foldername_extra = model + 'npart'+str(npart) + 'nsig'+str(sigma)
		cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/FilamentData/' + cachedir_foldername_extra + '/'
		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'Conditions_check.p'
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'
		# Read critical points from pickle and creates 3D coordinates of critical point positions
		Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		self.Unpack_filament_data(CP_coordinate_data)
		CP_3DPos = []
		for j in range(self.CritPointXpos):
			CP_3DPos.append(np.column_stack(self.CritPointXpos[j], self.CritPointYpos[j], self.CritPointZpos[j]))
		self.CP_3DPos = np.array(CP_3DPos)

		Filter_coordinates = self.CP_type >= 2 	# Filters out CP of type 1 or less. Filaments only connect to type 2 and 3
		self.CP_3DPos = self.CP_3DPos[Filter_coordinates]
		self.Number_filaments_connecting_to_CP = self.Number_filaments_connecting_to_CP[Filter_coordinates]
		
	def Unpack_filament_data(self, CP_coord):
		""" Unpacks critial point data from the read filament data module. """
		self.CritPointXpos = CP_coord[0]
		self.CritPointYpos = CP_coord[1]
		self.CritPointZpos = CP_coord[2]
		self.CP_type = CP_coord[3]
		self.CP_persistent_pair = CP_coord[4] 
		self.Critpts_filamentID = CP_coord[5] 
		self.CP_id_of_connecting_filament = CP_coord[6] 
		self.Number_filaments_connecting_to_CP = CP_coord[7]

	def CP_in_halo(self, HaloCord, CPCord, radius):
		norm_axis = 1
		if len(CPCord) == 1:
			norm_axis = 0
		Diff_len = np.linalg.norm(HaloCord - CPCord, axis=norm_axis)**2
		InHaloes = Diff_len <= radius**2
		return InHaloes

	def Check_cps(self):
		NumCP_in_halo = []
		CP_ids_in_halo = []
		for i in range(len(self.HaloPos)):
			InHalo = self.CP_in_halo(self.HaloPos[i], self.CP_3DPos, self.VirRadius[i])
			OK_cps = np.where(InHalo)[0]
			NumCP_in_halo.append(len(OK_cps))
			CP_ids_in_halo.append(OK_cps)
		return np.array(NumCP_in_halo), np.array(CP_ids_in_halo)



def Argument_parser():
	""" Parses optional argument when program is run from the command line """
	#print 'Run python code with -h argument to see extra optional arguments'
	parser = argparse.ArgumentParser()
	# Optional arguments
	parser.add_argument("-Model", "--Model", help="Determines which model to run."\
	 				+ "Models may be: lcdm, symmX (X=A,B,C,D) or fofrY (Y=4,5,6). Use argument 'all' to run all models. Runs none by default."\
	 				+ "Do not run seperate models and 'all' at once! Defaults at lcdm", type=str, default='lcdm')
	parser.add_argument("-Nparts", "--NumberParticles", help="Run with a set number of particles. Default 64.", type=int, default=64)
	parser.add_argument("-Nsigma", "--Nsigma_disperse", help="Sigma value used for DisPerSE. Default at 3.", type=int, default=3)
	parser.add_argument("-filetype", "--Filetype", help="Filetype the figures are saved in. Can be pdf or png (default)", type=str, default='png')
	# Parse arguments
	args = parser.parse_args()
	# Some checks
	Model_ok = False
	Model_disperse = False
	Model_check = ['lcdm', 'symmA', 'symmB', 'symmC', 'symmD', 'fofr4', 'fofr5', 'fofr6', 'all']
	Models_to_be_run = []
	if args.Model:
		for models in Model_check:
			if args.Model == models:
				Model_ok = True
				if args.Model == 'all':
					Models_to_be_run = [Model_check[i] for i in range(len(Model_check)-1)]
				else:
					Models_to_be_run.append(args.Model)
		if not Model_ok:
			raise ValueError('Model input name %s not correctly set -Model argument.' %args.Model)
	return args, Models_to_be_run

if __name__ == '__main__':
	savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
	# Parse the arguments
	parsed_arguments, Models_to_be_run = Argument_parser()
	p_model = parsed_arguments.Model
	N_parts = parsed_arguments.NumberParticles
	N_sigma = parsed_arguments.Nsigma_disperse
	Filetype = parsed_arguments.Filetype
	print '====== Information ======'
	print 'Running for model:', p_model
	print 'Number of particles: ', N_parts
	print 'Sigma value: ', N_sigma
	print '========================='

	for p_model in Models_to_be_run:
		Instance = AnalyseCritPts(p_model, N_parts, N_sigma)
		NumCP_in_halo, CP_ids_in_halo = Instance.Check_cps()
		