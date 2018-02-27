import numpy as np
import matplotlib.pyplot as plt
import os
import OtherFunctions as OF
plt.rcParams.update({'font.size': 11})	# Change fontsize of figures. Default = 10

class Histogram_Comparison():
	def __init__(self, savefile, foldername, savefile_directory, filetype, redshift, LCDM=False, SymmA=False, SymmB=False, nsigComparison=False):
		self.savefile = savefile
		self.LCDM_check = LCDM
		self.SymmA_check = SymmA
		self.SymmB_check = SymmB
		self.filetype = filetype

		self.ParticleComparison = False
		self.ModelComparison = False
		self.SigmaComparison = False

		self.results_dir = os.path.join(savefile_directory, foldername)
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
	def __init__(self, savefile, foldername, savefile_directory, filetype, redshift, nPart):
		self.savefile = savefile
		self.filetype = filetype

		self.ParticleComparison = False
		self.ModelComparison = False
		self.SigmaComparison = False

		self.Legends = ['$\mathregular{\Lambda}$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD', 'fofr4', 'fofr5', 'fofr6']

		if filetype == '.png':
			foldername += 'PNG/'
		elif filetype == '.pdf':
			foldername += 'PDF/'

		self.results_dir = os.path.join(savefile_directory, foldername)
		if not os.path.isdir(self.results_dir) and savefile == 1:
			os.makedirs(self.results_dir)

		self.nParticles = nPart

		self.s = 0.8 	# For figure rescaling etc. Change as you wish.

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

	def Compute_errors(self, databins):
		""" Computes the errorbars of the length data. Assumes Poisson distribution. """
		error = [np.sqrt(data) for data in databins]
		return np.array(error)

	def Check_ok_parameters(self, parameter):
		ok_scales = np.array(['normal', 'logx', 'logy', 'loglog'])
		is_there_ok = ok_scales == parameter
		if not is_there_ok.any():
			raise ValueError("logscale parameter not entered properly!. Use either normal, logx, logy or loglog.")
		

	def Call_plot(self, xdata, ydata, xlabel, ylabel, legend, logscale='normal', style='-', title='None'):
		""" 
		Calls plotting, xdata and ydata must be same size 
		Extra parameter logscal determines the plotting method
		normal = plot(), logx = semilogx(), etc.
		"""
		if len(xdata) != len(ydata):
			raise ValueError("xdata and ydata not of same length!")
		self.Check_ok_parameters(logscale)

		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		if logscale == 'normal':
			for i in range(len(xdata)):
				plt.plot(xdata[i], ydata[i], style)
		elif logscale == 'logx':
			for i in range(len(xdata)):
				plt.semilogx(xdata[i], ydata[i], style)
		elif logscale == 'logy':
			for i in range(len(xdata)):
				plt.semilogy(xdata[i], ydata[i], style)
		elif logscale == 'loglog':
			for i in range(len(xdata)):
				plt.loglog(xdata[i], ydata[i], style)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.legend(legend)
		plt.title('') if title == 'None' else plt.title(title)
		return figure

	def Call_plot_sameX(self, xdata, ydata, xlabel, ylabel, legend, logscale='normal', style='-', title='None'):
		""" 
		Calls plotting, xdata to be common lengths/bins for all ydata
		Extra parameter logscal determines the plotting method
		normal = plot(), logx = semilogx(), etc.
		"""
		self.Check_ok_parameters(logscale)

		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		if logscale == 'normal':
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], style)
		elif logscale == 'logx':
			for i in range(len(ydata)):
				plt.semilogx(xdata, ydata[i], style)
		elif logscale == 'logy':
			for i in range(len(ydata)):
				plt.semilogy(xdata, ydata[i], style)
		elif logscale == 'loglog':
			for i in range(len(ydata)):
				plt.loglog(xdata, ydata[i], style)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.legend(legend)
		plt.title('') if title == 'None' else plt.title(title)
		return figure

	def Plot_differences(self, xdata, ydata, xlabel, ylabel, legend, logscale='normal', style='-', title='None', diff='rel'):
		""" Plots relative or absolute differences. Base data assumed to be the first element of xdata and ydata """
		if len(xdata) != len(ydata):
			raise ValueError("xdata and ydata not of same length!")
		self.Check_ok_parameters(logscale)
		
		deltas = []
		if diff == 'rel':
			for i in range(1, len(xdata)):
				deltas.append(self.relative_deviation(ydata, i))
		elif diff == 'abs':
			for i in range(1, len(xdata)):
				deltas.append(ydata[i] - ydata[0])
		else:
			raise ValueError("Argument diff not properly set! Use either diff='rel' or diff='abs'")

		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		if logscale == 'normal':
			plt.plot(xdata[0], np.zeros(len(xdata[0])), style)
			for i in range(1, len(xdata)):
				plt.plot(xdata[i], deltas[i-1], style)
		elif logscale == 'logx':
			plt.semilogx(xdata[0], np.zeros(len(xdata[0])), style)
			for i in range(1, len(xdata)):
				plt.semilogx(xdata[i], deltas[i-1], style)
		elif logscale == 'logy':
			plt.semilogy(xdata[0], np.zeros(len(xdata[0])), style)
			for i in range(1, len(xdata)):
				plt.semilogy(xdata[i], deltas[i-1], style)
		elif logscale == 'loglog':
			plt.loglog(xdata[0], np.zeros(len(xdata[0])), style)
			for i in range(1, len(xdata)):
				plt.loglog(xdata[i], deltas[i-1], style)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.legend(legend)
		plt.title('') if title == 'None' else plt.title(title)
		return figure

	def Plot_differences_sameX(self, xdata, ydata, xlabel, ylabel, legend, logscale='normal', style='-', title='None', diff='rel'):
		""" 
		Plots relative or absolute differences. Base data assumed to be the first element of xdata and ydata 
		Xdata assumed to be common for all ydata.
		"""
		self.Check_ok_parameters(logscale)
		
		deltas = []
		if diff == 'rel':
			for i in range(1, len(ydata)):
				deltas.append(self.relative_deviation(ydata, i))
		elif diff == 'abs':
			for i in range(1, len(ydata)):
				deltas.append(ydata[i] - ydata[0])
		else:
			raise ValueError("Argument diff not properly set! Use either diff='rel' or diff='abs'")

		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		if logscale == 'normal':
			plt.plot(xdata, np.zeros(len(xdata)), style)
			for i in range(1, len(ydata)):
				plt.plot(xdata, deltas[i-1], style)
		elif logscale == 'logx':
			plt.semilogx(xdata, np.zeros(len(xdata)), style)
			for i in range(1, len(ydata)):
				plt.semilogx(xdata, deltas[i-1], style)
		elif logscale == 'logy':
			plt.semilogy(xdata, np.zeros(len(xdata)), style)
			for i in range(1, len(ydata)):
				plt.semilogy(xdata, deltas[i-1], style)
		elif logscale == 'loglog':
			plt.loglog(xdata, np.zeros(len(xdata)), style)
			for i in range(1, len(ydata)):
				plt.loglog(xdata, deltas[i-1], style)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.legend(legend)
		plt.title('') if title == 'None' else plt.title(title)
		return figure

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
		lengths = []
		for fils_lens in FilLengths:
			temp_dist = []
			lengths_arr = np.linspace(np.min(fils_lens), np.max(fils_lens), 1000)
			for lens in lengths_arr:
				Number_count = len(np.where(fils_lens >= lens)[0])
				temp_dist.append(float(Number_count))
			distribution.append(temp_dist)
			lengths.append(lengths_arr)
		FilLen_massfunc = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(distribution)):
			plt.semilogx(lengths[i], distribution[i])
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('$\mathregular{N(>L)}$')
		plt.legend(self.Legends)

		# Loglog scale of the number filaments larger than a given length
		FilLen_massfunc_loglog = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(distribution)):
			plt.loglog(lengths[i], distribution[i])
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('$\mathregular{N(>L)}$')
		plt.legend(self.Legends)

		# Density of the number of filaments larger than a given length. 
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(distribution)):
			plt.semilogx(lengths[i], np.array(distribution[i])/(256.0**3))
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('$\mathregular{N(>L)}/V$')
		plt.legend(self.Legends)

		# Density of the number of filaments larger than a given length. Loglogscale
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density_loglog = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(distribution)):
			plt.loglog(lengths[i], np.array(distribution[i])/(256.0**3))
		plt.xlabel('Filament length - [Mpc/h]')
		plt.ylabel('$\mathregular{N(>L)}/V$')
		plt.legend(self.Legends)

		# Relative difference of the lengths. Base model is LCDM.
		RelDiff_length = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		plt.semilogx(lengths[0], np.zeros(len(lengths[0])))
		for i in range(1,len(distribution)):
			delta = self.relative_deviation(np.array(distribution), i)
			plt.semilogx(lengths[i], delta)
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

		# Histogram lengths in bin-line form - Logarithmic x-scale
		length_bins_logX = []
		length_bin_values_logX = []
		for fillens in FilLengths:
			bins_, binval_, binstd_ = OF.Bin_numbers_logX(fillens, fillens, binnum=80)
			length_bins_logX.append(bins_)
			length_bin_values_logX.append(binval_)
		# Computes errors
		Poisson_errors = []
		for val in length_bin_values_logX:
			Poisson_errors.append(self.Compute_errors(val))

		# Plotting, with dots
		LengthHistComparison_digitized_logX = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins_logX)):
			plt.plot(length_bins_logX[i], length_bin_values_logX[i], 'o-')
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)

		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot_logX = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins_logX)):
			plt.plot(length_bins_logX[i], length_bin_values_logX[i], '-')
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)

		# Errorbar of the lengths
		LengthHistComparison_erorbar_logX = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins_logX)):
			plt.errorbar(length_bins_logX[i], length_bin_values_logX[i], Poisson_errors[i])
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)
		
		# Without dots indicating bin location, LOGLOG scale
		LengthHistComparison_digitized_nodot_logX_loglog = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bins_logX)):
			plt.loglog(length_bins_logX[i], length_bin_values_logX[i], '-')
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

		#### SEPARATING SYMMMETRON AND f(R) MODELS
		Sep_FilLengths_S = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		
		Sep_FilLengths_S = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		
		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(ConnectedHistComparison, 'Number_Connected_Filaments')
			self.savefigure(LengthHistComparison, 'Filament_lengths')
			self.savefigure(FilLen_massfunc, 'Filament_lengths_massfunction')
			self.savefigure(FilLen_massfunc_loglog, 'Filament_lengths_massfunction_loglog')
			self.savefigure(FilLen_massfunc_density, 'Filament_lengths_massfunction_density')
			self.savefigure(FilLen_massfunc_density_loglog, 'Filament_lengths_massfunction_density_loglog')
			self.savefigure(RelDiff_length, 'Filament_lengths_relative_difference')
			self.savefigure(LengthHistComparison_digitized, 'Filament_lengths_digitized')
			self.savefigure(LengthHistComparison_digitized_nodot, 'Filament_lengths_digitized_nodot')
			self.savefigure(LengthHistComparison_digitized_logX, 'Filament_lengths_digitized_logX')
			self.savefigure(LengthHistComparison_digitized_nodot_logX, 'Filament_lengths_digitized_nodot_logX')
			self.savefigure(LengthHistComparison_erorbar_logX, 'Filament_lengths_errorbar')
			self.savefigure(LengthHistComparison_digitized_nodot_logX_loglog, 'Filament_lengths_digitized_nodot_logX_loglog')
			self.savefigure(Cumulative_plot, 'Cumulative_plot_length')
			self.savefigure(RelDiff_cumulative_length, 'Relative_difference_cumulative_length')
		else:
			print 'Done! No figures saved.'

	def Compare_disperse_data_clean(self, Nconnections, FilLengths, FilPts):
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
		Minima = np.array([np.min(FilLengths[i]) for i in range(len(FilLengths))])
		Maxima = np.array([np.max(FilLengths[i]) for i in range(len(FilLengths))])
		lengths = np.linspace(np.min(Minima), np.max(Maxima), 1000)
		distribution = []
		#lengths = []
		for fils_lens in FilLengths:
			temp_dist = []
			#lengths_arr = np.linspace(np.min(fils_lens), np.max(fils_lens), 1000)
			for lens in lengths:
				Number_count = len(np.where(fils_lens >= lens)[0])
				temp_dist.append(float(Number_count))
			distribution.append(np.array(temp_dist))
			#lengths.append(lengths_arr)

		FilLen_massfunc = self.Call_plot_sameX(lengths, distribution, 'Filament length - [Mpc/h]', '$\mathregular{N(>L)}$', self.Legends, logscale='logx', style='-')

		# Loglog scale of the number filaments larger than a given length
		FilLen_massfunc_loglog = self.Call_plot_sameX(lengths, distribution, 'Filament length - [Mpc/h]', '$\mathregular{N(>L)}$', self.Legends, logscale='loglog')
		
		# Density of the number of filaments larger than a given length. 
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density = self.Call_plot_sameX(lengths, np.array(distribution)/(256.0**3), 'Filament length - [Mpc/h]', '$\mathregular{N(>L)}/V$', self.Legends, logscale='logx')
	
		# Density of the number of filaments larger than a given length. Loglogscale
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density_loglog = self.Call_plot_sameX(lengths, np.array(distribution)/(256.0**3), 'Filament length - [Mpc/h]', '$\mathregular{N(>L)}/V$', self.Legends, logscale='loglog')

		# Relative difference of the lengths. Base model is LCDM.
		Deltas = []
		for i in range(1, len(distribution)):
			Deltas.append(self.relative_deviation(np.array(distribution), i))
		Deltas = np.array(Deltas)

		RelDiff_length = self.Plot_differences_sameX(lengths, distribution, 'Filament length - [Mpc/h]', 'Relative difference of N(>L)', self.Legends, logscale='logx')

		# Histogram lengths in bin-line form
		length_bins = OF.Get_common_bin(FilLengths, binnum=40)
		length_bin_values = []
		length_bin_error = []
		for fillens in FilLengths:
			#bins_, binval_, binstd_ = OF.Bin_numbers(fillens, fillens, binnum=40)
			#length_bins.append(bins_)
			binval_, binstd_ = OF.Bin_numbers_common(fillens, fillens, length_bins, std='poisson')
			length_bin_values.append(binval_)
			length_bin_error.append(binstd_)

		# With dots indicating bin location
		LengthHistComparison_digitized = self.Call_plot_sameX(length_bins, length_bin_values, 'Filament Length  - [Mpc/h]', 'Number of filaments', self.Legends, style='o-')

		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot = self.Call_plot_sameX(length_bins, length_bin_values, 'Filament Length  - [Mpc/h]', 'Number of filaments', self.Legends, style='-')


		# Histogram lengths in bin-line form - Logarithmic x-scale
		length_bins_logX = OF.Get_common_bin_logX(FilLengths, binnum=80)
		length_bin_values_logX = []
		length_bin_std_logX = []
		for fillens in FilLengths:
			binval_, binstd_ = OF.Bin_numbers_common(fillens, fillens, length_bins_logX)
			#length_bins_logX.append(bins_)
			length_bin_values_logX.append(binval_)
			length_bin_std_logX.append(binstd_)
		
		# Computes errors, poisson error
		Poisson_errors = []
		for val in length_bin_values_logX:
			Poisson_errors.append(self.Compute_errors(val))

		# Plotting, with dots
		LengthHistComparison_digitized_logX = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, 'Filament Length  - [Mpc/h]', 'Number of filaments', self.Legends, style='o-')

		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot_logX = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, 'Filament Length  - [Mpc/h]', 'Number of filaments', self.Legends, style='-')

		# Errorbar of the lengths
		LengthHistComparison_erorbar_logX = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(length_bin_values_logX)):
			plt.errorbar(length_bins_logX, length_bin_values_logX[i], Poisson_errors[i])
		plt.xlabel('Filament Length  - [Mpc/h]')
		plt.ylabel('Number of filaments')
		plt.legend(self.Legends)
		
		# Without dots indicating bin location, LOGLOG scale
		LengthHistComparison_digitized_nodot_logX_loglog = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, 'Filament Length  - [Mpc/h]', 'Number of filaments',
																		  self.Legends, logscale='loglog', style='-')

		# Cumulative distribution of the lengths
		cumulative_list = []
		for i in range(len(length_bin_values)):
			Cumulative = np.cumsum(length_bin_values[i], dtype=np.float32)
			cumulative_list.append(Cumulative)

		Cumulative_plot = self.Call_plot_sameX(length_bins, cumulative_list, 'Filament Length - [Mpc/h]', 'Number of filaments', self.Legends, title='Cumulative distribution')

		# Relative difference of the cumulative distribution, x-value is the same, using LCDM
		RelDiff_cumulative_length = self.Plot_differences_sameX(length_bins, cumulative_list, 'Filament length - [Mpc/h]', 'Relative difference of cumulative distribution', self.Legends)

		#### SEPARATING SYMMMETRON AND f(R) MODELS. Data includes LCDM as default
		# Length binned histogram, nodot
		Symm_legends = ['$\mathregular{\Lambda}$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD']
		fofr_legends = ['$\mathregular{\Lambda}$CDM', 'fofr4', 'fofr5', 'fofr6']
		#Symm_length_bins = length_bins_logX[0:5]
		Symm_length_values = length_bin_values_logX[0:5]
		#fofr_length_bins = np.array([length_bins_logX[i] for i in [0,5,6,7]])
		fofr_length_values = np.array([length_bin_values_logX[i] for i in [0,5,6,7]])
		Sep_FilLengths_S = self.Call_plot_sameX(length_bins_logX, Symm_length_values, 'Filament Length  - [Mpc/h]', 'Number of filaments', Symm_legends, style='-')
		Sep_FilLengths_F = self.Call_plot_sameX(length_bins_logX, fofr_length_values, 'Filament Length  - [Mpc/h]', 'Number of filaments', fofr_legends, style='-')

		# 'massfunction' of lengths
		Distribution_symm = distribution[0:5]
		Distribution_fofr = np.array([distribution[i] for i in [0,5,6,7]])
		Lengths_symm = lengths[0:5]
		Lengths_fofr = np.array([lengths[i] for i in [0,5,6,7]])

		Distribution_lengths_S = self.Call_plot_sameX(lengths, Distribution_symm, 'Filament Length  - [Mpc/h]', '$\mathregular{N(>L)}$', Symm_legends, logscale='logx')
		Distribution_lengths_F = self.Call_plot_sameX(lengths, Distribution_fofr, 'Filament Length  - [Mpc/h]', '$\mathregular{N(>L)}$', fofr_legends, logscale='logx')
		Distribution_lengths_S_loglog = self.Call_plot_sameX(lengths, Distribution_symm, 'Filament Length  - [Mpc/h]', '$\mathregular{N(>L)}$', Symm_legends, logscale='loglog')
		Distribution_lengths_F_loglog = self.Call_plot_sameX(lengths, Distribution_fofr, 'Filament Length  - [Mpc/h]', '$\mathregular{N(>L)}$', fofr_legends, logscale='loglog')


		# Relative differences
		Sep_RelDiff_length_S = self.Plot_differences_sameX(lengths, Distribution_symm, 'Filament length - [Mpc/h]', 'Relative difference of N(>L)', Symm_legends, logscale='logx')
		Sep_RelDiff_length_F = self.Plot_differences_sameX(lengths, Distribution_fofr, 'Filament length - [Mpc/h]', 'Relative difference of N(>L)', fofr_legends, logscale='logx')

		# Absolute Differences in Number of filaments and not N(>L)
		Sep_AbsDiff_Number_S = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, 'Filament length - [Mpc/h]', 'Absolute difference of N (filaments)', Symm_legends, diff='abs')
		Sep_AbsDiff_Number_F = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, 'Filament length - [Mpc/h]', 'Absolute difference of N (filaments)', fofr_legends, diff='abs')
		
		# Relative differences in Number of filaments and not N(>L) 
		Sep_RelDiff_Number_S = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, 'Filament length - [Mpc/h]', 'Relative difference of N (filaments)', Symm_legends)
		Sep_RelDiff_Number_F = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, 'Filament length - [Mpc/h]', 'Relative difference of N (filaments)', fofr_legends)

		# Absolute differences in Number of filaments and not N(>L). Logscale x
		Sep_AbsDiff_Number_S_logx = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, 'Filament length - [Mpc/h]', 'Absolute difference of N (filaments)', Symm_legends, logscale='logx', diff='abs')
		Sep_AbsDiff_Number_F_logx = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, 'Filament length - [Mpc/h]', 'Absolute difference of N (filaments)', fofr_legends, logscale='logx', diff='abs')
		
		# Relative differences in Number of filaments and not N(>L). Logscale x
		Sep_RelDiff_Number_S_logx = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, 'Filament length - [Mpc/h]', 'Relative difference of N (filaments)', Symm_legends, logscale='logx')			
		Sep_RelDiff_Number_F_logx = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, 'Filament length - [Mpc/h]', 'Relative difference of N (filaments)', fofr_legends, logscale='logx')

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(ConnectedHistComparison, 'Number_Connected_Filaments')
			self.savefigure(LengthHistComparison, 'Filament_lengths')
			self.savefigure(FilLen_massfunc, 'Length_distribution')
			self.savefigure(FilLen_massfunc_loglog, 'Length_distribution_loglog')
			self.savefigure(FilLen_massfunc_density, 'Density_distribution')
			self.savefigure(FilLen_massfunc_density_loglog, 'Density_distribution_loglog')
			self.savefigure(RelDiff_length, 'Filament_lengths_relative_difference')
			self.savefigure(LengthHistComparison_digitized, 'Filament_lengths_digitized')
			self.savefigure(LengthHistComparison_digitized_nodot, 'Filament_lengths_digitized_nodot')
			self.savefigure(LengthHistComparison_digitized_logX, 'Filament_lengths_digitized_logX')
			self.savefigure(LengthHistComparison_digitized_nodot_logX, 'Filament_lengths_digitized_nodot_logX')
			self.savefigure(LengthHistComparison_erorbar_logX, 'Filament_lengths_errorbar')
			self.savefigure(LengthHistComparison_digitized_nodot_logX_loglog, 'Filament_lengths_digitized_nodot_logX_loglog')
			self.savefigure(Cumulative_plot, 'Cumulative_plot_length')
			self.savefigure(RelDiff_cumulative_length, 'Relative_difference_cumulative_length')
			self.savefigure(Sep_FilLengths_S, 'Filament_lengths_cSymmetron')
			self.savefigure(Sep_FilLengths_F, 'Filament_lengths_cFofr')
			self.savefigure(Sep_RelDiff_length_S, 'Filament_lengths_relative_difference_cSymmetron')
			self.savefigure(Sep_RelDiff_length_F, 'Filament_lengths_relative_difference_cFofr')
			self.savefigure(Distribution_lengths_S, 'Length_distribution_cSymmetron')
			self.savefigure(Distribution_lengths_F, 'Length_distribution_cFofr')
			self.savefigure(Distribution_lengths_S_loglog, 'Length_distribution_cSymmetron_loglog')
			self.savefigure(Distribution_lengths_F_loglog, 'Length_distribution_cFofr_loglog')
			self.savefigure(Sep_AbsDiff_Number_S, 'Number_filaments_absDiff_cSymmetron')
			self.savefigure(Sep_AbsDiff_Number_F, 'Number_filaments_absDiff_cFofr')
			self.savefigure(Sep_RelDiff_Number_S, 'Number_filaments_relDiff_cSymmetron')
			self.savefigure(Sep_RelDiff_Number_F, 'Number_filaments_relDiff_cFofr')
			self.savefigure(Sep_AbsDiff_Number_S_logx, 'Number_filaments_absDiff_LOGX_cSymmetron')
			self.savefigure(Sep_AbsDiff_Number_F_logx, 'Number_filaments_absDiff_LOGX_cFofr')
			self.savefigure(Sep_RelDiff_Number_S_logx, 'Number_filaments_relDiff_LOGX_cSymmetron')
			self.savefigure(Sep_RelDiff_Number_F_logx, 'Number_filaments_relDiff_LOGX_cFofr')
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

	def Filament_distances(self, filament_pos, fil_len):
		""" 
		Computes the distances between the filaments. 
		Distances are between the center of the filaments.
		"""
		def get_binning(fildist, bins):
			index_bin = np.digitize(fildist, bins)
			bin_value = np.array([len(fildist[i == index_bin]) for i in range(len(bins))])
			return bin_value

		if len(filament_pos) != len(fil_len):
			raise ValueError("Dataset of filament position and filament lengths are not equal!")
		# Get filament centers
		Midpoints = []
		for j in range(len(filament_pos)):
			midpt_temp = []
			for i in range(len(filament_pos[j])):
				midpoint = OF.find_middle_point(filament_pos[j][i], fil_len[j][i])
				midpt_temp.append(midpoint)
			Midpoints.append(np.array(midpt_temp))
		# Compute distances from filament centers
		Distances = []
		for mid in Midpoints:
			dist = OF.get_filament_distances(mid)
			Distances.append(dist)

		# Bins the distances, assumes same bins
		bins = np.linspace(0, 256.0, 70)
		bin_values = []
		for i in range(len(Distances)):
			binval_temp = []
			for j in range(len(Distances[i])):
				values = get_binning(Distances[i][j], bins)
				binval_temp.append(values)
			bin_values.append(np.array(binval_temp))

		# Find averages of the bins
		bin_averages = []
		for i in range(len(bin_values)):
			bin_avg = OF.Average_hist(bins, bin_values[i])
			bin_averages.append(bin_avg)

		FilDistancePlot = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(bin_averages)):
			plt.plot(bins, bin_averages[i], 'o-')
		plt.xlabel('Filament distances - [Mpc/h]')
		plt.ylabel('$N$ filaments')
		plt.legend(self.Legends)

		FilDistancePlot_nodot = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		for i in range(len(bin_averages)):
			plt.plot(bins, bin_averages[i], '-')
		plt.xlabel('Filament distances - [Mpc/h]')
		plt.ylabel('$N$ filaments')
		plt.legend(self.Legends)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(FilDistancePlot, 'FilamentDistances')
			self.savefigure(FilDistancePlot_nodot, 'FilamentDistances_nodot')
		else:
			print 'Done! No figures saved.'