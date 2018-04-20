import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import OtherFunctions as OF
import PlotFunctions as pf
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
	def __init__(self, savefile, foldername, savefile_directory, filetype, nPart, Nsigma=3, redshift=0):
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
		sigma_name_folder = 'Sigma'+str(Nsigma) + '/'
		foldername += sigma_name_folder

		self.results_dir = os.path.join(savefile_directory, foldername)
		if not os.path.isdir(self.results_dir) and savefile == 1:
			os.makedirs(self.results_dir)

		self.nParticles = nPart

		self.s = 0.7 	# For figure rescaling etc. Change as you wish.

		self.Plot_colors_symm = ['k', 'orange', 'g', 'r', 'olive']
		self.Plot_colors_fofr = ['k', 'purple', 'y', 'b']
		self.Plot_colors_all = ['k', 'orange', 'g', 'r', 'olive', 'purple', 'y', 'b']
		self.Linestyles = ['-', '--', '-.', ':', (0, (5, 10))]

	def relative_deviation(self, data, index):
		""" 
		Computes relative deviation of a given physical quantity.
		This function assumes that the base model is in the first element.
		"""
		#delta = (data[index].astype(np.float32) - data[0].astype(np.float32))/data[0].astype(np.float32)
		delta = (data[index] - data[0])/data[0]
		return delta

	def savefigure(self, figure, name):
		""" Function that calls savefig based on figure instance and filename. """
		if type(name) != str:
			raise ValueError('filename not a string!')
		figure.savefig(self.results_dir + name + self.filetype, bbox_inches='tight', dpi=100)

	def Compute_errors(self, databins):
		""" Computes the errorbars of the length data. Assumes Poisson distribution. """
		error = [np.sqrt(data) for data in databins]
		return np.array(error)

	def Propagate_error_relNum(self, data_0, data, error_0, error):
		""" 
		Relative difference of N(>L) 
		Function to be F = (D[i] - D[0]/D[0]), with D = data, index 0 is LCDM
		Derivatives
		dFd(D[i]) = 1/D[0]
		dfD(D[0]) = -D[i]/D[0]**2
		"""
		dfdDi = 1/data_0
		dfdD0 = -data/data_0**2
		Prop_error = np.sqrt((error_0**2)*(dfdD0**2) + (error**2)*(dfdDi)**2)
		return Prop_error

	def Propagate_error_absNum(self, error_0, error):
		""" 
		Relative difference of the absolute difference Delta L = L2 - L1
		Let F = L2 - L1, L is the filament lengths of data 1 and 2, then
		dFdL1 = -1, dFdL2 = 1 
		"""
		Prop_error = np.sqrt(error_0**2 + error**2)
		return Prop_error

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

	def Call_plot_sameX(self, xdata, ydata, xlabel, ylabel, legend, colors, logscale='normal', style='-', title='None', Lengths=True):
		""" 
		Calls plotting, xdata to be common lengths/bins for all ydata
		Extra parameter logscal determines the plotting method
		normal = plot(), logx = semilogx(), etc.
		"""
		self.Check_ok_parameters(logscale)

		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		ax = plt.axes()
		if logscale == 'normal':
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], style, label=legend[i], color=colors[i])
		elif logscale == 'logx':
			for i in range(len(ydata)):
				plt.semilogx(xdata, ydata[i], style, label=legend[i], color=colors[i])
		elif logscale == 'logy':
			for i in range(len(ydata)):
				plt.semilogy(xdata, ydata[i], style, label=legend[i], color=colors[i])
		elif logscale == 'loglog':
			for i in range(len(ydata)):
				plt.loglog(xdata, ydata[i], style, label=legend[i], color=colors[i])
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		#plt.legend(legend)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		plt.title('') if title == 'None' else plt.title(title)
		max_yvalues = np.max(ydata)
		min_yvalues = np.min(ydata)
		if np.max(max_yvalues) > 1000 or np.min(min_yvalues) < 1:
			#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=True)
			#plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0,0))
			xfmt = plt.ScalarFormatter()
			xfmt.set_powerlimits((0,0))
			plt.gca().yaxis.set_major_formatter(xfmt)
		if Lengths:
			plt.xlim((1,np.max(xdata)))
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

	def Plot_differences_sameX(self, xdata, ydata, xlabel, ylabel, legend, colors, logscale='normal', style='-', title='None', diff='rel', Lengths=True):
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
		ax = plt.axes()
		if logscale == 'normal':
			plt.plot(xdata, np.zeros(len(xdata)), style, label='$\Lambda$CDM')
			for i in range(1, len(ydata)):
				plt.plot(xdata, deltas[i-1], style, label=legend[i-1], color=colors[i-1])
		elif logscale == 'logx':
			plt.semilogx(xdata, np.zeros(len(xdata)), style, label='$\Lambda$CDM')
			for i in range(1, len(ydata)):
				plt.semilogx(xdata, deltas[i-1], style, label=legend[i-1], color=colors[i-1])
		elif logscale == 'logy':
			plt.semilogy(xdata, np.zeros(len(xdata)), style, label='$\Lambda$CDM')
			for i in range(1, len(ydata)):
				plt.semilogy(xdata, deltas[i-1], style, label=legend[i-1], color=colors[i-1])
		elif logscale == 'loglog':
			plt.loglog(xdata, np.zeros(len(xdata)), style, label='$\Lambda$CDM')
			for i in range(1, len(ydata)):
				plt.loglog(xdata, deltas[i-1], style, label=legend[i-1], color=colors[i-1])
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		#plt.legend(legend)
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		plt.title('') if title == 'None' else plt.title(title)
		if Lengths:
			plt.xlim((1, np.max(xdata)))

		return figure

	def Plot_errobar_sameX(self, xdata, ydata, error, xlabel, ylabel, legend, colors, logscale='normal', fill_between=False, diff=False, limit_yaxis=True, Lengths=True):
		if len(ydata) != len(error):
			raise ValueError("ydata and error data not of same length!")
		self.Check_ok_parameters(logscale)

		Max_yval = np.array([np.max(ydata[i]) for i in range(len(ydata))])
		Min_yval = np.array([np.min(ydata[i]) for i in range(len(ydata))])
		figure = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		ax = plt.axes()
		if diff:
			plt.plot(xdata, np.zeros(len(xdata)), label='$\Lambda$CDM')
			if fill_between:
				plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0.3)
		if fill_between:
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], label=legend[i], color=colors[i])
				plt.fill_between(xdata, ydata[i]-error[i], ydata[i]+error[i], alpha=0.3, facecolor=colors[i])
			if limit_yaxis:
				plt.ylim((-1.0, 1))
		else:
			for i in range(len(ydata)):
				plt.errorbar(xdata, ydata[i], error[i], label=legend[i], color=colors[i])
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)

		if Lengths:
			plt.xlim(1, np.max(xdata))
		
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
		if logscale == 'logx':
			plt.xscale('log')
		elif logscale == 'logy':
			plt.yscale('log')
		elif logscale == 'loglog':
			plt.xscale('log')
			plt.yscale('log')
		#plt.legend(legend)
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
		# COMMONLY USED LABELS
		xlabel_len = 'Filament length - [Mpc/h]'

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
		lengths = np.linspace(np.min(Minima), np.max(Maxima), 500)
		distribution = []
		for fils_lens in FilLengths:
			temp_dist = []
			for lens in lengths:
				Number_count = len(np.where(fils_lens >= lens)[0])
				temp_dist.append(float(Number_count))
			distribution.append(np.array(temp_dist))
		# Error of these distribution values. Assumes Poisson distribution
		distribution_error = []
		for distribution_list in distribution:
			err_distribution = np.array([np.sqrt(number) for number in distribution_list])
			distribution_error.append(err_distribution)

		FilLen_massfunc = self.Call_plot_sameX(lengths, distribution, xlabel_len, '$\mathregular{N(>L)}$', self.Legends, self.Plot_colors_all, logscale='logx', style='-')

		# Loglog scale of the number filaments larger than a given length
		FilLen_massfunc_loglog = self.Call_plot_sameX(lengths, distribution, xlabel_len, '$\mathregular{N(>L)}$', self.Legends, self.Plot_colors_all, logscale='loglog')
		
		# Density of the number of filaments larger than a given length. 
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density = self.Call_plot_sameX(lengths, np.array(distribution)/(256.0**3), xlabel_len, '$\mathregular{N(>L)}/V$', self.Legends, self.Plot_colors_all, logscale='logx')
	
		# Density of the number of filaments larger than a given length. Loglogscale
		# Divided by box volume of particles (256^3), can use box size of disperse?
		FilLen_massfunc_density_loglog = self.Call_plot_sameX(lengths, np.array(distribution)/(256.0**3), xlabel_len, '$\mathregular{N(>L)}/V$', self.Legends, self.Plot_colors_all, logscale='loglog')

		# Relative difference of the lengths. Base model is LCDM.
		Deltas = []
		for i in range(1, len(distribution)):
			Deltas.append(self.relative_deviation(np.array(distribution), i))
		Deltas = np.array(Deltas)

		RelDiff_length = self.Plot_differences_sameX(lengths, distribution, xlabel_len, 'Relative difference of $N(>L)$', self.Legends[1:], self.Plot_colors_all[1:], logscale='logx')

		#num, bins = np.histogram(FilLengths[0], bins='fd')
		#binnum = len(bins)
		binnum = 30
		# Histogram lengths in bin-line form
		length_bins = OF.Get_common_bin(FilLengths, binnum=binnum)
		length_bin_values = []
		length_bin_error = []
		for fillens in FilLengths:
			binval_, binstd_ = OF.Bin_numbers_common(fillens, fillens, length_bins, std='poisson')
			length_bin_values.append(binval_)
			length_bin_error.append(binstd_)
		# With dots indicating bin location
		LengthHistComparison_digitized = self.Call_plot_sameX(length_bins, length_bin_values, xlabel_len, '$N$ filaments', self.Legends, self.Plot_colors_all, style='o-')

		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot = self.Call_plot_sameX(length_bins, length_bin_values, xlabel_len, '$N$ filaments', self.Legends, self.Plot_colors_all, style='-')

		# Histogram lengths in bin-line form - Logarithmic x-scale
		length_bins_logX = OF.Get_common_bin_logX(FilLengths, binnum=binnum)
		length_bin_values_logX = []
		length_bin_error_logX = []
		for fillens in FilLengths:
			binval_, binstd_ = OF.Bin_numbers_common(fillens, fillens, length_bins_logX, std='poisson')
			length_bin_values_logX.append(binval_)
			length_bin_error_logX.append(binstd_)
		# Plotting, with dots
		LengthHistComparison_digitized_logX = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, xlabel_len, '$N$ filaments', self.Legends, self.Plot_colors_all, style='o-')


		Bin_comparison = plt.figure()
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		plt.subplot(1,2,1)
		plt.plot(length_bins, length_bin_values[0], '-')
		plt.legend(['Linear scale'])
		plt.subplot(1,2,2)
		plt.plot(length_bins_logX, length_bin_values_logX[0])
		plt.legend(['Logarithmic scale'])
		plt.tight_layout()


		# Without dots indicating bin location
		LengthHistComparison_digitized_nodot_logX = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, xlabel_len, '$N$ filaments', self.Legends, self.Plot_colors_all, style='-')

		# Errorbar of the lengths
		LengthHistComparison_erorbar_logX = self.Plot_errobar_sameX(length_bins_logX, length_bin_values_logX, length_bin_error_logX,
																	 xlabel_len, '$N$ filaments', self.Legends, self.Plot_colors_all)
		
		# Without dots indicating bin location, LOGLOG scale
		LengthHistComparison_digitized_nodot_logX_loglog = self.Call_plot_sameX(length_bins_logX, length_bin_values_logX, xlabel_len, '$N$ filaments',
																		  self.Legends, self.Plot_colors_all, logscale='loglog', style='-')

		# Cumulative distribution of the lengths
		cumulative_list = []
		for i in range(len(length_bin_values)):
			Cumulative = np.cumsum(length_bin_values[i], dtype=np.float32)
			cumulative_list.append(Cumulative)

		Cumulative_plot = self.Call_plot_sameX(length_bins, cumulative_list, xlabel_len, 'Number of filaments', self.Legends, self.Plot_colors_all, title='Cumulative distribution')

		# Relative difference of N(>L) using the cumulative data
		Cumulative_distribution = np.array([len(FilLengths[i]) - cumulative_list[i] for i in range(len(FilLengths))])
		RelDiff_cumulative_reverse = self.Plot_differences_sameX(length_bins, Cumulative_distribution, xlabel_len, '$\mathregular{N(>L)}$', self.Legends[1:], self.Plot_colors_all[1:])

		# Computes propagated error on the relative difference of N(>L)
		Propagated_errors = np.array([self.Propagate_error_relNum(distribution[0], distribution[i], distribution_error[0], distribution_error[i]) for i in range(1,len(FilLengths))])
		Prop_error_reldiff_N = self.Plot_errobar_sameX(lengths, Deltas, Propagated_errors, xlabel_len, 'Relative difference of $\mathregular{N(>L)}$',
													 self.Legends[1:], self.Plot_colors_all[1:], fill_between=True, diff=True)

		RelDiff_num = [OF.relative_deviation_singular(length_bin_values_logX[0], length_bin_values_logX[i]) for i in range(1, len(FilLengths))]
		Prop_err_numLength = [OF.Propagate_error_reldiff(length_bin_values_logX[0], length_bin_values_logX[i], length_bin_error_logX[0], 
														length_bin_error_logX[i]) for i in range(1, len(FilLengths))]
		#### SEPARATING SYMMMETRON AND f(R) MODELS. Data includes LCDM as default
		Symm_legends = ['$\mathregular{\Lambda}$CDM', 'SymmA', 'SymmB', 'SymmC', 'SymmD']
		fofr_legends = ['$\mathregular{\Lambda}$CDM', 'fofr4', 'fofr5', 'fofr6']
		Symm_legends_only = ['SymmA', 'SymmB', 'SymmC', 'SymmD']
		fofr_legends_only = ['fofr4', 'fofr5', 'fofr6']
		
		Symm_length_values = length_bin_values_logX[0:5]
		fofr_length_values = np.array([length_bin_values_logX[i] for i in [0,5,6,7]])

		Symm_colors_only = self.Plot_colors_symm[1:]
		fofr_colors_only = self.Plot_colors_fofr[1:]
		# Length binned histogram, nodot
		xlim_len = (1, length_bins_logX[-1])
		Sep_FilLengths_S = self.Call_plot_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Number of filaments',
												Symm_legends, self.Plot_colors_symm, style='-')
		Sep_FilLengths_F = self.Call_plot_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Number of filaments',
												fofr_legends, self.Plot_colors_fofr, style='-')

		Sep_FilLengths_S_logx = pf.Call_plot_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Number of filaments', Symm_legends, self.Plot_colors_symm, 
													xscale='log', linestyles=self.Linestyles, legend_anchor=False, xlim=xlim_len)
		Sep_FilLengths_F_logx = pf.Call_plot_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Number of filaments', fofr_legends, self.Plot_colors_fofr, 
													xscale='log', linestyles=self.Linestyles, legend_anchor=False, xlim=xlim_len)
		
		Sep_FilLengths_S_loglog = self.Call_plot_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Number of filaments',
													 Symm_legends, self.Plot_colors_symm, style='-', logscale='loglog')
		Sep_FilLengths_F_loglog = self.Call_plot_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Number of filaments',
													 fofr_legends, self.Plot_colors_fofr, style='-', logscale='loglog')
		# Relative difference of N filaments, with propagated errors
		Set_FilLengths_properr_S_logx = pf.Call_plot_sameX(length_bins_logX, RelDiff_num[0:4], xlabel_len, '$(N_i - N_{\Lambda CDM})/N_{\Lambda CDM}$', 
															 Symm_legends_only, Symm_colors_only, error=Prop_err_numLength[0:4],
															 fillbetween=True, xscale='log', linestyles=self.Linestyles, reldiff=True)
		Set_FilLengths_properr_F_logx = pf.Call_plot_sameX(length_bins_logX, RelDiff_num[4:], xlabel_len, '$(N_i - N_{\Lambda CDM})/N_{\Lambda CDM}$',
															 fofr_legends_only, fofr_colors_only, error=Prop_err_numLength[4:],
															 fillbetween=True, xscale='log', linestyles=self.Linestyles, reldiff=True)
		
		# 'massfunction' of lengths
		Distribution_symm = distribution[0:5]
		Distribution_fofr = np.array([distribution[i] for i in [0,5,6,7]])
		RelDiff_distribution_symm = np.array([OF.relative_deviation_singular(distribution[0], Distribution_symm[i]) for i in range(1, len(Distribution_symm))])
		RelDiff_distribution_fofr = np.array([OF.relative_deviation_singular(distribution[0], Distribution_fofr[i]) for i in range(1, len(Distribution_fofr))])

		Lengths_symm = lengths[0:5]
		Lengths_fofr = np.array([lengths[i] for i in [0,5,6,7]])

		Distribution_lengths_S = self.Call_plot_sameX(lengths, Distribution_symm, xlabel_len, '$\mathregular{N(>L)}$', Symm_legends,
													 self.Plot_colors_symm, logscale='logx')
		Distribution_lengths_F = self.Call_plot_sameX(lengths, Distribution_fofr, xlabel_len, '$\mathregular{N(>L)}$', fofr_legends,
													 self.Plot_colors_fofr, logscale='logx')
		Distribution_lengths_S_loglog = self.Call_plot_sameX(lengths, Distribution_symm, xlabel_len, '$\mathregular{N(>L)}$', Symm_legends,
															 self.Plot_colors_symm, logscale='loglog')
		Distribution_lengths_F_loglog = self.Call_plot_sameX(lengths, Distribution_fofr, xlabel_len, '$\mathregular{N(>L)}$', fofr_legends,
															 self.Plot_colors_fofr, logscale='loglog')

		# Relative differences of N(>R)
		Sep_RelDiff_length_S = self.Plot_differences_sameX(lengths, Distribution_symm, xlabel_len, 'Relative difference of $N(>L)$', Symm_legends_only,
														 Symm_colors_only, logscale='logx')
		Sep_RelDiff_length_F = self.Plot_differences_sameX(lengths, Distribution_fofr, xlabel_len, 'Relative difference of $N(>L)$', fofr_legends_only,
														 fofr_colors_only, logscale='logx')
		# Relative difference of N(>R) with errorbar
		#Sep_RelDiff_length_error_S = self.Plot_errobar_sameX(lengths, RelDiff_distribution_symm, Propagated_errors[0:4], xlabel_len, 'Relative difference of $N(>L)$',
		#					 							Symm_legends_only, Symm_colors_only, fill_between=True, diff=True, logscale='logx')
		#Sep_RelDiff_length_error_F = self.Plot_errobar_sameX(lengths, RelDiff_distribution_fofr, Propagated_errors[4:], xlabel_len, 'Relative difference of $N(>L)$',
		#												 fofr_legends_only, fofr_colors_only, fill_between=True, diff=True, logscale='logx')
		Sep_RelDiff_length_error_S = pf.Call_plot_sameX(lengths, RelDiff_distribution_symm, xlabel_len, 'Relative difference of $N(>L)$', Symm_legends, 
														Symm_colors_only, error=Propagated_errors[:4], fillbetween=True, reldiff=True, xscale='log',
							 							linestyles=self.Linestyles, legend_anchor=False, ylim=(-1,1), xlim=xlim_len)
		Sep_RelDiff_length_error_F = pf.Call_plot_sameX(lengths, RelDiff_distribution_fofr, xlabel_len, 'Relative difference of $N(>L)$', fofr_legends, 
														fofr_colors_only, error=Propagated_errors[4:], fillbetween=True, reldiff=True, xscale='log',
														linestyles=self.Linestyles, legend_anchor=False, ylim=(-1,1), xlim=xlim_len)
		# Absolute Differences in Number of filaments and not N(>L)
		Sep_AbsDiff_Number_S = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, xlabel_len, '$\Delta N$ filaments', Symm_legends_only,
															Symm_colors_only, diff='abs')
		Sep_AbsDiff_Number_F = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, xlabel_len, '$\Delta N$ filaments', fofr_legends_only,
															fofr_colors_only, diff='abs')
		
		# Relative differences in Number of filaments and not N(>L) 
		Sep_RelDiff_Number_S = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
														 Symm_legends_only, Symm_colors_only)
		Sep_RelDiff_Number_F = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
														 fofr_legends_only, fofr_colors_only)

		# Absolute differences in Number of filaments and not N(>L). Logscale x
		Sep_AbsDiff_Number_S_logx = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, xlabel_len, '$\Delta N$ filaments', Symm_legends_only,
																Symm_colors_only, logscale='logx', diff='abs')
		Sep_AbsDiff_Number_F_logx = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, xlabel_len, '$\Delta N$ filaments', fofr_legends_only,
																fofr_colors_only, logscale='logx', diff='abs')
		
		# Relative differences in Number of filaments and not N(>L). Logscale x
		Sep_RelDiff_Number_S_logx = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
																Symm_legends_only, Symm_colors_only, logscale='logx')			
		Sep_RelDiff_Number_F_logx = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
																fofr_legends_only, fofr_colors_only, logscale='logx')

		# Relative differences in Number of filaments and not N(>L). Logscale x, including errorbars
		Sep_RelDiff_Number_error_S_logx = pf.Call_plot_sameX(length_bins_logX, RelDiff_num[:4], xlabel_len, 
															'$(N_i - N_{\Lambda\mathrm{CDM}})/N_{\Lambda\mathrm{CDM}}$',
															Symm_legends, Symm_colors_only, xscale='log', error=Prop_err_numLength[:4], 
															fillbetween=True, xlim=(1, np.max(length_bins_logX)), legend_anchor=False, reldiff=True,
															ylim=(-1,1), linestyles=self.Linestyles)
		Sep_RelDiff_Number_error_F_logx = pf.Call_plot_sameX(length_bins_logX, RelDiff_num[4:], xlabel_len, 
															'$(N_i - N_{\Lambda\mathrm{CDM}})/N_{\Lambda\mathrm{CDM}}$',
															fofr_legends, fofr_colors_only, xscale='log', error=Prop_err_numLength[4:],
															fillbetween=True, xlim=(1, np.max(length_bins_logX)), legend_anchor=False, reldiff=True,
															ylim=(-1,1), linestyles=self.Linestyles)

		# Relative differences in Number of filaments and not N(>L). LogLog scale
		Sep_RelDiff_Number_S_loglog = self.Plot_differences_sameX(length_bins_logX, Symm_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
																 Symm_legends_only, Symm_colors_only, logscale='loglog')			
		Sep_RelDiff_Number_F_loglog = self.Plot_differences_sameX(length_bins_logX, fofr_length_values, xlabel_len, 'Relative difference of $N$ (filaments)',
																 fofr_legends_only, fofr_colors_only, logscale='loglog')        

		# Relative differences of N(>L) with propagated errors
		Prop_err_symm = Propagated_errors[0:4]
		Prop_err_fofr = Propagated_errors[4:]
		Deltas_symm = Deltas[0:4]
		Deltas_fofr = Deltas[4:]
		Sep_RelDiff_length_PropErr_S = self.Plot_errobar_sameX(lengths, Deltas_symm, Prop_err_symm, xlabel_len, 'Relative difference of $N(>L)$',
															 Symm_legends_only, Symm_colors_only, fill_between=True, diff=True)
		Sep_RelDiff_length_PropErr_F = self.Plot_errobar_sameX(lengths, Deltas_fofr, Prop_err_fofr, xlabel_len, 'Relative difference of $N(>L)$',
															 fofr_legends_only, fofr_colors_only, fill_between=True, diff=True)
		
		# Absolute difference of N, with propagated errors
		Propagated_errors_lengths = np.array([self.Propagate_error_absNum(length_bin_error_logX[0], length_bin_error_logX[i]) for i in range(1, len(FilLengths))])
		Propp_err_lens_symm = Propagated_errors_lengths[0:4]
		Propp_err_lens_fofr = Propagated_errors_lengths[4:]
		AbsDiff_symm = np.array([Symm_length_values[i] - Symm_length_values[0] for i in range(1, len(Symm_length_values))])
		AbsDiff_fofr = np.array([fofr_length_values[i] - fofr_length_values[0] for i in range(1, len(fofr_length_values))])
		Sep_AbsDiff_Number_error_S = self.Plot_errobar_sameX(length_bins_logX, AbsDiff_symm, Propp_err_lens_symm, xlabel_len, '$\Delta N$ filaments',
																Symm_legends_only, Symm_colors_only, fill_between=True, diff=True, limit_yaxis=False)
		Sep_AbsDiff_Number_error_F = self.Plot_errobar_sameX(length_bins_logX, AbsDiff_fofr, Propp_err_lens_fofr, xlabel_len, '$\Delta N$ filaments',
																fofr_legends_only, fofr_colors_only, fill_between=True, diff=True, limit_yaxis=False)
		# Absolute difference of N, with propagated errors. LogX scale
		Sep_AbsDiff_Number_error_S_logx = self.Plot_errobar_sameX(length_bins_logX, AbsDiff_symm, Propp_err_lens_symm, xlabel_len, '$\Delta N$ filaments',
																Symm_legends_only, Symm_colors_only, logscale='logx', fill_between=True, diff=True, limit_yaxis=False)
		Sep_AbsDiff_Number_error_F_logx = self.Plot_errobar_sameX(length_bins_logX, AbsDiff_fofr, Propp_err_lens_fofr, xlabel_len, '$\Delta N$ filaments',
																fofr_legends_only, fofr_colors_only, logscale='logx', fill_between=True, diff=True, limit_yaxis=False)

		ConnectedHistComparison_subplot = plt.figure(figsize=(8,6))
		plt.gcf().set_size_inches((8*self.s, 6*self.s))
		ax = plt.subplot(1,2,1)
		connection_bins = OF.Get_common_bin(Nconnections, binnum=21)
		bin_val_connhist = []
		bin_std_connhist = []
		for i in range(0,5):
			binval, binstd = OF.Bin_numbers_common(Nconnections[i], Nconnections[i], connection_bins, std='poisson')
			bin_val_connhist.append(binval)
			bin_std_connhist.append(binstd)
			plt.plot(connection_bins, binval, color=self.Plot_colors_all[i])
			plt.fill_between(connection_bins, binval-binstd, binval+binstd, alpha=0.4, facecolor=self.Plot_colors_all[i])
		plt.legend(Symm_legends)
		plt.xscale('log')
		plt.yscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for i in [0,5,6,7]:
			binval, binstd = OF.Bin_numbers_common(Nconnections[i], Nconnections[i], connection_bins, std='poisson')
			if i != 0:
				bin_val_connhist.append(binval)
				bin_std_connhist.append(binstd)
			plt.plot(connection_bins, binval, color=self.Plot_colors_all[i])
			plt.fill_between(connection_bins, binval-binstd, binval+binstd, alpha=0.4, facecolor=self.Plot_colors_all[i])
		plt.legend(fofr_legends)
		plt.xscale('log')
		plt.yscale('log')
		ConnectedHistComparison_subplot.text(0.5, 0, 'Number connections', ha='center', fontsize=10)
		ConnectedHistComparison_subplot.text(0, 0.5, '$N$ filaments', ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()

		ConnectedHistComparison_subplot_reldiff = plt.figure(figsize=(5,4))
		ax = plt.subplot(1,2,1)
		for i in range(1,5):
			reldiff_chist = OF.relative_deviation_singular(bin_val_connhist[0], bin_val_connhist[i])
			reldiff_err_chist = OF.Propagate_error_reldiff(bin_val_connhist[0], bin_val_connhist[i], bin_std_connhist[0], bin_std_connhist[i])
			plt.plot(connection_bins, reldiff_chist, color=self.Plot_colors_all[i])
			plt.fill_between(connection_bins, reldiff_chist-reldiff_err_chist, reldiff_chist+reldiff_err_chist, alpha=0.4, facecolor=self.Plot_colors_all[i])
		plt.legend(Symm_legends_only)
		plt.xscale('log')
		ax.set_ylim(-2,2)
		#plt.yscale('log')
		ax2 = plt.subplot(1,2,2, sharey=ax)
		plt.setp(ax2.get_yticklabels(), visible=False)
		for i in range(5,8):
			reldiff_chist = OF.relative_deviation_singular(bin_val_connhist[0], bin_val_connhist[i])
			reldiff_err_chist = OF.Propagate_error_reldiff(bin_val_connhist[0], bin_val_connhist[i], bin_std_connhist[0], bin_std_connhist[i])
			plt.plot(connection_bins, reldiff_chist, color=self.Plot_colors_all[i])
			plt.fill_between(connection_bins, reldiff_chist-reldiff_err_chist, reldiff_chist+reldiff_err_chist, alpha=0.4, facecolor=self.Plot_colors_all[i])
		plt.legend(fofr_legends_only)
		plt.xscale('log')
		ax2.set_ylim(-2,2)
		ConnectedHistComparison_subplot_reldiff.text(0.5, 0, 'Number connections', ha='center', fontsize=10)
		ConnectedHistComparison_subplot_reldiff.text(0, 0.5, 'Relative difference of number connections', ha='center', va='center', rotation='vertical', fontsize=10)
		plt.tight_layout()	

		#### Using gridspec plotting
		NumLen_symm_logx_GRIDSPEC = pf.Do_gridspec_sameX(length_bins_logX, [np.array(Symm_length_values)/1000.0], [RelDiff_num[:4]],
													xlabel_len, r'$N \times 1000$', 'Relative difference', Symm_legends, self.Plot_colors_symm,
													Secerror=[Prop_err_numLength[:4]], xscale='log', linestyles=self.Linestyles, reldiff=True,
													fillbetween=True, xlim=xlim_len, ylim_diff=(-0.75, 0.75), rowcol=[2,1], legend_anchor=False)
		NumLen_fofr_logx_GRIDSPEC = pf.Do_gridspec_sameX(length_bins_logX, [np.array(fofr_length_values)/1000.0], [RelDiff_num[4:]],
													xlabel_len, r'$N \times 1000$', 'Relative difference', fofr_legends, self.Plot_colors_fofr,
													Secerror=[Prop_err_numLength[4:]], xscale='log', linestyles=self.Linestyles, reldiff=True,
													fillbetween=True, xlim=xlim_len, ylim_diff=(-0.3, 0.3), rowcol=[2,1], legend_anchor=False)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			self.savefigure(Bin_comparison, 'Binning_comparison')
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
			self.savefigure(RelDiff_cumulative_reverse, 'Filament_lengths_relative_difference_CumulativeData')
			self.savefigure(Prop_error_reldiff_N, 'Filament_lengths_relative_difference_propagatedError')
			## Seperating stuffers
			self.savefigure(Sep_FilLengths_S, 'Filament_lengths_cSymmetron')
			self.savefigure(Sep_FilLengths_F, 'Filament_lengths_cFofr')
			self.savefigure(Sep_FilLengths_S_logx, 'Filament_lengths_cSymmetron_logX')
			self.savefigure(Sep_FilLengths_F_logx, 'Filament_lengths_cFofr_logX')
			self.savefigure(Sep_FilLengths_S_loglog, 'Filament_lengths_cSymmetron_loglog')
			self.savefigure(Sep_FilLengths_F_loglog, 'Filament_lengths_cFofr_loglog')
			self.savefigure(Set_FilLengths_properr_S_logx, 'Filament_length_RelativeDiff_cSymmetron_logx')
			self.savefigure(Set_FilLengths_properr_F_logx, 'Filament_length_RelativeDiff_cFofr_logx')
			self.savefigure(Sep_RelDiff_length_S, 'Filament_lengths_relative_difference_cSymmetron')
			self.savefigure(Sep_RelDiff_length_F, 'Filament_lengths_relative_difference_cFofr')
			self.savefigure(Sep_RelDiff_length_error_S, 'Filament_lengths_cumulativeRev_reldiff_cSymmetron')
			self.savefigure(Sep_RelDiff_length_error_F, 'Filament_lengths_cumulativeRev_reldiff_cFofr')
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
			self.savefigure(Sep_RelDiff_Number_error_S_logx, 'Number_filaments_relDiff_ERRORBAR_LOGX_cSymmetron')
			self.savefigure(Sep_RelDiff_Number_error_F_logx, 'Number_filaments_relDiff_ERRORBAR_LOGX_cFofr')
			self.savefigure(Sep_RelDiff_Number_S_loglog, 'Number_filaments_relDiff_LOGLOG_cSymmetron')
			self.savefigure(Sep_RelDiff_Number_F_loglog, 'Number_filaments_relDiff_LOGLOG_cFofr')
			self.savefigure(Sep_RelDiff_length_PropErr_S, 'Filament_lengths_relative_difference_PropErr_cSymmetron')
			self.savefigure(Sep_RelDiff_length_PropErr_F, 'Filament_lengths_relative_difference_PropErr_cFofr')
			self.savefigure(Sep_AbsDiff_Number_error_S, 'Filament_lengths_absolute_difference_PropErr_cSymmetron')
			self.savefigure(Sep_AbsDiff_Number_error_F, 'Filament_lengths_absolute_difference_PropErr_cFofr')
			self.savefigure(Sep_AbsDiff_Number_error_S_logx, 'Filament_lengths_absolute_difference_PropErr_logX_cSymmetron')
			self.savefigure(Sep_AbsDiff_Number_error_F_logx, 'Filament_lengths_absolute_difference_PropErr_logX_cFofr')
			self.savefigure(ConnectedHistComparison_subplot, 'Number_Connected_Filaments_subplot')
			self.savefigure(ConnectedHistComparison_subplot_reldiff, 'Number_Connected_Filaments_subplot_RelDiff')
			## Gridspec stuff
			self.savefigure(NumLen_symm_logx_GRIDSPEC, 'NumberLength_logX_cSymmetron_gridspec')
			self.savefigure(NumLen_fofr_logx_GRIDSPEC, 'NumberLength_logX_cFofr_gridspec')
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