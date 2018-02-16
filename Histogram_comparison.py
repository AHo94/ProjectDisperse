import numpy as np
import matplotlib.pyplot as plt
import os
import OtherFunctions as OF

class Histogram_Comparison():
	def __init__(self, savefile, savefigDirectory, savefile_directory, filetype, redshift, LCDM=False, SymmA=False, SymmB=False, nsigComparison=False):
		self.savefile = savefile
		self.LCDM_check = LCDM
		self.SymmA_check = SymmA
		self.SymmB_check = SymmB
		self.filetype = filetype

		self.ParticleComparison = False
		self.ModelComparison = False
		self.SigmaComparison = False

		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir) and savefile == 1:
			os.makedirs(self.results_dir)

		if LCDM and not SymmA and not SymmB:
			self.ModelFilename = 'LCDM' + 'z' + str(redshift)
			self.ParticleComparison = True
		elif not LCDM and SymmA and not SymmB:
			self.ModelFilename = 'SymmA' + 'z' + str(redshift)
			self.ParticleComparison = True
		elif not LCDM and not SymmA and SymmB:
			self.ModelFilename = 'SymmB' + 'z' + str(redshift)
			self.ParticleComparison = True
		elif LCDM and SymmA and not SymmB == 0:
			self.ModelFilename = 'LCDM_SymmB' + 'z' + str(redshift)
			self.ModelComparison = True
		elif LCDM and not SymmA and SymmB:
			self.ModelFilename = 'LCDM_SymmB' + 'z' + str(redshift)
			self.ModelComparison = True
		elif LCDM and SymmA and not SymmB:
			self.ModelFilename = 'LCDM_SymmB' + 'z' + str(redshift)
			self.ModelComparison = True
		elif LCDM and SymmA and SymmB:
			self.ModelFilename = 'LCDM_SymmA_SymmB' + 'z' + str(redshift)
		else:
			raise ValueError('At least one model must be set to compare!')

		if nsigComparison:
			self.ModelComparison = 0
			self.ParticleComparison = 0
			self.SigmaComparison = 1

			self.ModelFilename += 'SigmaCompare'

	def Run(self, NumberConnections, FilamentLengths, NPointsPerFilament, nPart=False):
		if not nPart:
			self.ModelFilename += self.filetype
		elif nPart == 64:
			self.ModelFilename += '64Part' + self.filetype
		elif nPart == 128:
			self.ModelFilename += '128Part' + self.filetype
		elif nPart == 256:
			self.ModelFilename += '256Part' + self.filetype
		elif nPart == 512:
			self.ModelFilename += '512Part' + self.filetype
		self.nParticles = nPart

		if type(NumberConnections) != list:
			raise ValueError('Argument NumberConnections must be a list!')
		elif type(FilamentLengths) != list:
			raise ValueError('Argument FilamentLengths must be a list!')
		
		if len(NumberConnections) < 2:
			raise ValueError('Nothing to compare because NumberConnections has less than two arrays!')
		elif len(FilamentLengths) < 2:
			raise ValueError('Nothing to compare because FilamentLengths has less than two arrays!')

		if len(NumberConnections) != len(FilamentLengths):
			raise ValueError('List size of NumberConnections and FilamentLengths must be equal!')

		if not self.SigmaComparison:
			for i in range(len(NumberConnections)-1):
				if len(NumberConnections[i+1]) < len(NumberConnections[i]):
					raise ValueError('List of NumberConnections must be in order to increasing number of particles')
				elif len(FilamentLengths[i+1]) < len(FilamentLengths[i]):
					raise ValueError('List of FilamentLengths must be in order to increasing number of particles')
				elif len(NPointsPerFilament[i+1]) < len(NPointsPerFilament[i]):
					raise ValueError('List of Number of points per filament must be in order to increasing number of particles')
				

		self.NumberConnections = NumberConnections
		self.FilamentLengths = FilamentLengths
		self.NPointsPerFilament = NPointsPerFilament

		self.N = len(self.NumberConnections)
		self.Check_Number_Comparisons()
		self.Plot_Histograms_particleComparison()


	def Check_Number_Comparisons(self):
		if self.ParticleComparison:
			if self.N == 2:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles']
			elif self.N == 3:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', '$\mathregular{256^3}$ particles']
			elif self.N == 4:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', \
									'$\mathregular{256^3}$ particles', '$\mathregular{512^3}$ particles']
		elif self.ModelComparison:
			self.LegendText = []
			if self.LCDM_check:
				self.LegendText.append('LCDM')
			if self.SymmA_check:
				self.LegendText.append('Symm_A')
			if self.SymmB_check:
				self.LegendText.append('Symm_B')
		elif self.SigmaComparison:
			self.LegendText = ['$\mathregular{\sigma=3}$', '$\mathregular{\sigma=4}$', '$\mathregular{\sigma=5}$']


	def Plot_Histograms_particleComparison(self):
		alphas = [0.4, 0.5, 0.6, 0.7]
		#alphas = [0.6, 0.5, 0.4, 0.3]
		if self.nParticles == 512:
			subsample_text = 's'
		else:
			subsample_text = ' subsample'

		ConnectedHistComparison = plt.figure()
		plt.hold(True)
		for i in range(self.N):
			DataMin = min(self.NumberConnections[i])
			DataMax = max(self.NumberConnections[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(self.NumberConnections[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of connected filaments')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number filament connections \n with '+str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)		
		plt.legend(self.LegendText)
		plt.hold(False)
	
		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(self.N):
			plt.hist(self.FilamentLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of filament length \n with ' +str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)
		plt.legend(self.LegendText)
		plt.hold(False)

		NPointsHistComparison = plt.figure()
		plt.hold(True)
		for i in range(self.N):
			DataMin = min(self.NPointsPerFilament[i])
			DataMax = max(self.NPointsPerFilament[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(self.NPointsPerFilament[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of points per filament')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number of datapoints per filament \n with ' +str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)
		plt.legend(self.LegendText)
		plt.hold(False)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison' + self.ModelFilename)
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison' + self.ModelFilename)
			NPointsHistComparison.savefig(self.results_dir + 'HistNPointsComparison' + self.ModelFilename)
		elif self.savefile == 2:
			print 'Done! No histograms to plot.'
		else:
			plt.show()
		plt.close('all')

	def Convergence_tests(self, Nconnections, FilLengths, NptsPerFilament):
		if len(Nconnections) != len(FilLengths) or len(Nconnections) != len(NptsPerFilament) or len(FilLengths) != len(NptsPerFilament):
			raise ValueError('Lists containing the histograms are not of equal length!')

		alphas = [0.65, 0.6, 0.55, 0.5, 0.45, 0.4]
		Legends = ['$\mathregular{\sigma=5}, 64^3$ part subsample', '$\mathregular{\sigma=4},\
					 64^3$ part subsample', '$\mathregular{\sigma=3}, 64^3$ part subsample', \
					'$\mathregular{\sigma=5}, 512^3$ part', '$\mathregular{\sigma=4}, 512^3$ part', '$\mathregular{\sigma=3}, 512^3$ part']
		N = len(Nconnections)
		ConnectedHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			DataMin = min(Nconnections[i])
			DataMax = max(Nconnections[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(Nconnections[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of connected filaments')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number connections per filament')
		plt.legend(Legends)
		plt.hold(False)

		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			plt.hist(FilLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of filament lengths')
		plt.legend(Legends)
		plt.hold(False)

		NPointsHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			DataMin = min(NptsPerFilament[i])
			DataMax = max(NptsPerFilament[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(NptsPerFilament[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of points per filament')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number data points per filament')
		plt.legend(Legends)
		plt.hold(False)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison_64and512Part' + self.ModelFilename)
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison_64and512Part' + self.ModelFilename)
			NPointsHistComparison.savefig(self.results_dir + 'HistNPointsComparison_64and512Part' + self.ModelFilename)
		elif self.savefile == 2:
			print 'Done! No histograms to plot.'
		else:
			plt.show()
		plt.close('all')

	def Sigma_plot_comparison(self, Nconnections, FilLengths, NptsPerFilament, nsigma):
		if len(Nconnections) != len(FilLengths) or len(Nconnections) != len(NptsPerFilament) or len(FilLengths) != len(NptsPerFilament):
			raise ValueError('Lists containing the histograms are not of equal length!')

		Legends = ['$64^3$ part subsample', '$128^3$ part subsample', '$256^3$ part subsample', '$512^3$ particles']
		alphas = [0.7, 0.6, 0.5, 0.4]
		N = len(Nconnections)
		ConnectedHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			DataMin = min(Nconnections[i])
			DataMax = max(Nconnections[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(Nconnections[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of connected filaments')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number connections per filament for $\mathregular{\sigma=}$' + str(nsigma))
		plt.legend(Legends)
		plt.hold(False)

		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			plt.hist(FilLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of filament lengths for $\mathregular{\sigma=}$' + str(nsigma))
		plt.legend(Legends)
		plt.hold(False)

		NPointsHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			DataMin = min(NptsPerFilament[i])
			DataMax = max(NptsPerFilament[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(NptsPerFilament[i], align='mid', rwidth=1, bins=BinList, normed=False, alpha=alphas[i], histtype='step')
		plt.xlabel('Number of points per filament')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number data points per filament for $\mathregular{\sigma=}$' + str(nsigma))
		plt.legend(Legends)
		plt.hold(False)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			ConnectedHistComparison.savefig(self.results_dir 
				+ 'HistNumConnectedFilamentsComparison_AllParticles_nsig' + str(nsigma) + self.filetype)
			LengthHistComparison.savefig(self.results_dir 
				+ 'HistLengthComparison_AllParticles_nsig' + str(nsigma) + self.filetype)
			NPointsHistComparison.savefig(self.results_dir 
				+ 'HistNPointsComparison_AllParticles_nsig' + str(nsigma) + self.filetype)
		elif self.savefile == 2:
			print 'Done! No histograms to plot.'
		else:
			plt.show()
		plt.close('all')

class CompareModels():
	def __init__(self, savefile, savefigDirectory, savefile_directory, filetype, redshift, nPart):
		self.savefile = savefile
		self.filetype = filetype

		self.ParticleComparison = False
		self.ModelComparison = False
		self.SigmaComparison = False

		self.Legends = ['$\mathregular{\Lambda}$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6']

		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir) and savefile == 1:
			os.makedirs(self.results_dir)

		self.nParticles = nPart

		self.s = 0.6 	# For figure rescaling etc. Change as you wish.

	def relative_deviation(self, data, index):
		""" 
		Computes relative deviation of a given physical quantity.
		This function assumes that the base model is in the first element.
		"""
		delta = (data[index] - data[0])/data[0]
		return delta

	def savefigure(self, figure, name):
		""" Function that calls savefig based on figure instance and filename. """
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(self.results_dir + name + self.filetype)

	def Compare_disperse_data(self, Nconnections, FilLengths, FilPts):
		""" 
		Compares basic properties of the data given by disperse.
		Relative deviation compares models with respect to the LCDM model.
		The function assumes that LCDM data is the first in the data list.
		"""
		N = len(Nconnections)
		# Histogram for number of filament connections per filament
		ConnectedHistComparison = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(N):
			DataMin = min(Nconnections[i])
			DataMax = max(Nconnections[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(Nconnections[i], align='mid', rwidth=1, bins=BinList, normed=False, histtype='step')
		plt.xlabel('Number of connected filaments per filament')
		plt.ylabel('Number of filaments')
		#plt.xscale('log')
		plt.title('Histogram comparison of number of filament connections for each filmanet')
		plt.legend(self.Legends)
		# Relative deviation of number of filament connections
		"""
		NconComparison = plt.figure()
		for i in range(1,N):
			Delta = relative_deviation(Nconnections, i)
			DataMin = min(Delta)
			DataMax = max(Delta)
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(Delta, align='mid', rwidth=1, bins=BinList, normed=False, histtype='step')
		plt.xlabel('Number of connected filaments per filament')
		plt.ylabel('Number of filaments')
		plt.title('Histogram comparison of number of filament connections for each filmanet')
		plt.legend(self.Legends)
		"""

		# Histogram of filament length comparison
		LengthHistComparison = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(N):
			DataMin = min(FilLengths[i])
			DataMax = max(FilLengths[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(FilLengths[i], align='mid', rwidth=1, bins=BinList, normed=False, histtype='step')
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('Number of filaments')
		#plt.xscale('log')
		plt.title('Histogram comparison of filament lengths')
		plt.legend(self.Legends)

		# Number of filaments larger than a given length: N(>L)
		#lengths = np.linspace(np.min(np.min(FilLengths)), np.max(np.max(FilLengths)), 1000)
		distribution = []
		for fils_lens in FilLengths:
			temp_dist = []
			lengths = np.linspace(np.min(fils_lens), np.max(fils_lens), 1000)
			for lens in lengths:
				Number_count = len(np.where(fils_lens >= lens)[0])
				temp_dist.append(float(Number_count))
			distribution.append(temp_dist)
		FilLen_massfunc = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(distribution)):
			plt.semilogx(lengths, distribution[i])
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('$\mathregular{N(>L)}$')
		plt.legend(self.Legends)

		# Relative difference of the lengths. Base model is LCDM.
		RelDiff_length = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		plt.semilogx(lengths, np.zeros(len(lengths)))
		for i in range(1,len(distribution)):
			delta = self.relative_deviation(np.array(distribution), i)
			plt.semilogx(lengths, delta)
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('Relative difference of N(>L)')
		plt.legend(self.Legends)

		# Histogram lengths in bin-line form
		length_bins = []
		length_bin_values = []
		for fillens in FilLengths:
			bins_, binval_, binstd_ = OF.Bin_numbers(fillens, fillens, binnum=40)
			length_bins.append(bins_)
			length_bin_values.append(binval_)

		LengthHistComparison_digitized = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins)):
			plt.plot(length_bins[i], length_bin_values[i], 'o-')
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)

		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins)):
			plt.plot(length_bins[i], length_bin_values[i], '-')
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)

		# Cumulative distribution of the lengths
		cumulative_list = []
		for i in range(len(length_bin_values)):
			Cumulative = np.cumsum(length_bin_values[i], dtype=np.float32)
			cumulative_list.append(Cumulative)

		Cumulative_plot = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(cumulative_list)):
			plt.plot(length_bins[i], cumulative_list[i])
		plt.xlabel('Filament Length - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.title('Cumulative distribution')
		plt.legend(self.Legends)

		# Relative difference of the cumulative distribution, x-value is the same, using LCDM
		RelDiff_cumulative_length = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		plt.plot(length_bins[0], np.zeros(len(length_bins[0])))
		for i in range(1, len(cumulative_list)):
			delta = self.relative_deviation(cumulative_list,i)
			plt.plot(length_bins[0], delta)
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('Relative difference of cumulative distribution')
		plt.legend(self.Legends)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(ConnectedHistComparison, 'Number_Connected_Filaments')
			self.savefigure(LengthHistComparison, 'Filament_lengths')
			self.savefigure(FilLen_massfunc, 'Filament_lengths_massfunction')
			self.savefigure(RelDiff_length, 'Filament_lengths_relative_difference')
			self.savefigure(LengthHistComparison_digitized, 'Filament_lengths_digitized')
			self.savefigure(LengthHistComparison_digitized_nodot, 'Filament_lengths_digitized_nodot')
			self.savefigure(Cumulative_plot, 'Cumulative_plot_length')
			self.savefigure(RelDiff_cumulative_length, 'Relative_difference_cumulative_length')
		else:
			print 'Done! No figures saved.'

	def Compare_particle_properties(self, NptsPerFilament, NumParts, FilMass):
		""" Compares properties where the number of particles per filaments are computed """
		# Number of particles per filament for a given distance threshold
		NumPart_histogram = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(NumParts)):
			plt.hist(NumParts[i], align='mid', rwidth=1, bins=60, normed=False, histtype='step')
		plt.xlabel('Number of particles per filament')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)

		# Filament mass larger than a given mass: N(>M)
		Mass_array = np.linspace(np.min(np.min(FilMass)), np.max(np.max(FilMass)), 1000)
		Mass_distribution = []
		for Filament_mass in FilMass:
			temp_dist_mass = []
			for mass in Mass_array:
				Number_count = len(np.where(Filament_mass >= mass)[0])
				temp_dist_mass.append(float(Number_count))
			Mass_distribution.append(np.array(temp_dist_mass))
		FilMass_massfunc = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(Mass_distribution)):
			plt.semilogx(Mass_array, Mass_distribution[i])
		plt.xlabel('Filament mass [kg]')
		plt.ylabel('\mathregular{N(>M)}')
		plt.legend(self.Legends)

		# Relative difference of the masses. Base model is LCDM
		RelDiff_mass = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		plt.semilogx(Mass_array, np.zeros(len(Mass_array)))
		for i in range(1,len(Mass_distribution)):
			delta = self.relative_deviation(np.array(Mass_distribution), i)
			plt.semilogx(Mass_array, delta)
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('Relative difference of N(>M)')
		plt.legend(self.Legends)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(NumPart_histogram, 'NumberParticles_per_filament')
			self.savefigure(FilMass_massfunc, 'Filament_mass_massfunction')
			self.savefigure(RelDiff_mass, 'Filament_mass_relative_difference')
		else:
			print 'Done! No figures saved.'