### Set comment on the two below to plot. Only use if running on euclid21 or other institute computers. 
#import matplotlib
#matplotlib.use('Agg')
# Basic stuff
import os
import time
import argparse
import cPickle as pickle
from collections import Counter
from functools import partial
import multiprocessing as mp

# Numpy, matplotlib and scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from scipy import spatial
from scipy.interpolate import griddata

# Disable warnings from matplotlib	
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# Own modules
import BoundaryChecker
import FilamentMasking
import ReadGadgetFile
import MaskCritPts
import FilamentsPerSigma
import Histogram_comparison as HComp
import ReadFilamentData
import InterpolateDensity
import ParticlesPerFilament

class Disperse_Plotter():
	"""
	This class does all the analysis from the data created by DisPerSE. This includes reading of filament points, critical points etc.
	If provided, it will also analyse the dark matter particles.
	Note that this class only solves for one model at a time. 
	"""
	def __init__(self, savefile, savefigDirectory, nPart, model, redshift, SigmaArg=False):
		self.savefile = savefile
		self.PlotChecker = 0
		self.nPart = nPart
		if self.nPart == 64:
			self.nPart_text = '$\mathregular{64^3}$'
		elif self.nPart == 128:
			self.nPart_text = '$\mathregular{128^3}$'
		elif self.nPart == 256:
			self.nPart_text = '$\mathregular{256^3}$'
		elif self.nPart == 512:
			self.nPart_text = '$\mathregular{512^3}$'
		else:
			raise ValueError('Argument nPart must be 64, 128, 256 or 512!')

		if type(savefigDirectory) != str:
			raise ValueError('Argument savefigDirectory must be a string!')

		# Creates folder when saving figure
		#script_dir = os.path.dirname(__file__) 	# Script running directory
		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir) and savefile == 1:
			os.makedirs(self.results_dir)

		if type(redshift) != int and type(redshift) != float:
			raise ValueError('Value of redshift must be an integer or float')

		if model == 'LCDM':
			self.ModelFilename = 'LCDM' + 'z' + str(redshift) + filetype
		elif model == 'SymmA':
			self.ModelFilename = 'SymmA' + 'z' + str(redshift) + filetype
		elif model == 'SymmB':
			self.ModelFilename = 'SymmB' + 'z' + str(redshift) + filetype
		elif model == 'SymmC':
			self.ModelFilename = 'SymmC' + 'z' + str(redshift) + filetype
		elif model == 'SymmD':
			self.ModelFilename = 'SymmD' + 'z' + str(redshift) + filetype
		elif model == 'fofr4':
			self.ModelFilename = 'fofr4' + 'z' + str(redshift) + filetype
		elif model == 'fofr5':
			self.ModelFilename = 'fofr5' + 'z' + str(redshift) + filetype
		elif model == 'fofr6':
			self.ModelFilename = 'fofr6' + 'z' + str(redshift) + filetype
		else:
			raise ValueError('Model name not properly set')
		self.model = model

		self.SigmaTitle = ''
		if SigmaComparison:
			if not SigmaArg:
				self.SigmaTitle = '$\mathregular{\sigma} = $3'
			else:
				self.SigmaTitle = '$\mathregular{\sigma} = $' + str(SigmaArg)
		if not SigmaArg:
			self.Alternative_sigmaTitle = '$\mathregular{\sigma} = $3'
		else:
			self.Alternative_sigmaTitle = '$\mathregular{\sigma} = $' + str(SigmaArg)


		if SigmaArg:
			self.SigmaArg = SigmaArg
		else:
			self.SigmaArg = False

		# A check that the funciton 'solve' has been called
		self.SolveRun = False

	def Filter_filaments(self, Sigma_threshold):
		CP_indices_to_be_deleted = np.where(self.Persistence_nsigmas <= Sigma_threshold)[0]
		
		CP_persistence_pairs = []
		for CP_index in CP_indices_to_be_deleted:
			Ppair = self.Neighbours_CP[CP_index]
			if self.Number_filaments_connecting_to_CP[Ppair] != 0:
				CP_persistence_pairs.append(Ppair)
		#CP_persistence_pairs = np.array(self.Neighbours_CP)[CP_indices_to_be_deleted]		
		CP_indices_to_be_deleted = np.unique(np.concatenate([CP_indices_to_be_deleted, CP_persistence_pairs]))
		Filaments_to_be_deleted = []
		
		for indices in CP_indices_to_be_deleted:
			CP_IDs = self.CP_id_of_connecting_filament[indices]
			for IDs in CP_IDs:
				if IDs == self.Neighbours_CP[indices]:
					Filament_ID_deleted = self.Critpts_filamentID[indices][np.where(self.CP_id_of_connecting_filament[indices] == IDs)[0]][0]
					Filaments_to_be_deleted.append(Filament_ID_deleted)
	
		"""
		for indices in CP_indices_to_be_deleted:
			Fil_IDs = self.Critpts_filamentID[indices]
			for IDs in Fil_IDs:
				Filaments_to_be_deleted.append(IDs)
		"""
		# Check duplicates
		#Duplicates = [k for k,v in Counter(Filaments_to_be_deleted).items() if v>1]
		#Duplicate_indices = np.where()
		#print Duplicate_indices
		#print len(Duplicates)
		Filaments_to_be_deleted = np.unique(np.array(Filaments_to_be_deleted))
		#print len(Filaments_to_be_deleted)

		self.xdimPos = np.delete(self.xdimPos, Filaments_to_be_deleted)
		self.ydimPos = np.delete(self.ydimPos, Filaments_to_be_deleted)
		self.zdimPos = np.delete(self.zdimPos, Filaments_to_be_deleted)
		self.FilamentPos = np.delete(self.FilamentPos, Filaments_to_be_deleted)

		self.NFils = len(self.xdimPos)

	def Filter_filaments2(self, Sigma_threshold):
		maximas = np.where(np.array(self.CP_type) == 3)[0]
		maxima_nsigmas = np.array(self.Persistence_nsigmas)[maximas]
		filaments_included = []
		for i in range(self.NFils-1):
			CP_maxima_index = self.PairIDS[i][1]
			if maxima_nsigmas[CP_maxima_index] >= Sigma_threshold:
				filaments_included.append(i)
		print np.array(filaments_included)
		print len(filaments_included)

	def Number_filament_connections(self):
		""" Computes the number of filament connections one filament has"""
		print 'Computing number connections'
		self.NumFilamentConnections = []
		for Connected_CP in self.PairIDS:
			counter = 0
			for CP_ids in Connected_CP:
				counter += self.Number_filaments_connecting_to_CP[CP_ids]
			self.NumFilamentConnections.append(counter)

		self.NumFilamentConnections = np.asarray(self.NumFilamentConnections)

	def Plot_Figures(self, filename, ndim=3):
		""" All plots done in this function	"""
		def savefigure(figure, name):
			""" Function that calls savefig based on figure instance and filename. """
			if type(name) != str:
				raise ValueError('filename not a string!')
			figure.savefig(self.results_dir + name + self.ModelFilename)

		print 'Plotting'
		if FilamentColors:
			ColorArray = np.linspace(0,1,110)
			ColorMap2D = np.array([np.mean(zvalues) for zvalues in self.zdimPos])
			ColorMapLength = np.array([lengths for lengths in self.LengthSplitFilament])
			if IncludeSlicing:
				ColorMap2DCutOff = np.array([np.mean(zvalues) for zvalues in self.CutOffzDim])
				ColorMap2DMasked = np.array([np.mean(zvalues) for zvalues in self.zdimMasked])
				ColorMapLengthCutOff = np.array([lengths for lengths in self.CutOffLengths])
				ColorMapLengthMasked = np.array([lengths for lengths in self.MaskedLengths])
		else:
			ColorArray = None
		
		if IncludeUnits:
			LegendText = ' - [Mpc/h]'
		else:
			LegendText = ' - [Percentage of box size]'

		if HistogramPlots:
			# Histogram of number of filament connections
			###### NEEDS FIXING!! Done.
			NumMin = min(self.NumFilamentConnections)
			NumMax = max(self.NumFilamentConnections)
			BinSize = (NumMax - NumMin)/(0.5) + 1
			Bin_list = np.linspace(NumMin, NumMax, BinSize)
			ConnectedHist = plt.figure()
			plt.hist(self.NumFilamentConnections, align='left', rwidth=1, bins=Bin_list)
			plt.xlabel('Number of connected filaments')
			plt.ylabel('Number of occurances')
			plt.title('Histogram of number connections between the critical points \n with %.2f critical points. ' \
						%self.NcritPts + '. Using '+ self.nPart_text+'$\mathregular{^3}$ particles.')
			plt.xlim([NumMin, NumMax])
			#plt.xticks(self.NumFilamentConnections)

			# Histogram of filament lengths
			LenMin = min(self.FilLengths)
			LenMax = max(self.FilLengths)
			BinSize_size = (LenMax - LenMin)/(0.5) + 1
			#BinList_FilLengths = np.linspace(LenMin, LenMax, BinSize_size)
			FilamentLengthsHist = plt.figure()
			plt.hist(self.FilLengths, align='left', rwidth=1, bins=600, histtype='step')#, bins=BinList_FilLengths)
			#plt.hold("on")
			#self.FilLengths.sort()
			#fit = stats.norm.pdf(self.FilLengths, np.mean(self.FilLengths), np.std(self.FilLengths))
			#plt.plot(self.FilLengths, fit)
			plt.xlabel('Length of filaments' + LegendText)
			plt.ylabel('Number of occurances')
			plt.title('Histogram of filament lengths with ' + self.nPart_text + '$\mathregular{^3}$ particles.')
			plt.xlim([LenMin, LenMax])

			# Histogram of number of points each filament has
			PtsMin = min(self.NFilamentPoints)
			PtsMax = max(self.NFilamentPoints)
			BinSize_FilPts = (PtsMax - PtsMin)/(0.5) + 1
			BinList_FilPts = np.linspace(PtsMin, PtsMax, BinSize_FilPts)
			FilamentPtsHis = plt.figure()
			plt.hist(self.NFilamentPoints, align='left', rwidth=1, bins=BinList_FilPts)
			plt.xlabel('Number of position points for a filament')
			plt.ylabel('Number of occurances')
			plt.title('Histogram of number of filament points with ' + self.nPart_text + '$\mathregular{^3}$ particles')
			#plt.xticks(self.NFilamentPoints)
			"""
			if IncludeDMParticles == 1:
				NPartFilBinMin = np.min(self.Particles_per_filament)
				NPartFilBinMax = np.max(self.Particles_per_filament)
				BinSize_NPartFil = (NPartFilBinMax - NPartFilBinMin)/(0.5) + 1
				BinList_NPartFil = np.linspace(NPartFilBinMin, NPartFilBinMax, BinSize_NPartFil)
				NumParticlesFilamentHist = plt.figure()
				plt.hist(self.Particles_per_filament, align='left', rwidth=1, bins=600)
				plt.xlabel('Number of particles per filament')
				plt.ylabel('Number of occurances')
				plt.title('Histogram of number of particles per filament with ' + self.nPart_text + '$\mathregular{^3}$ particles')
			"""
			
		if ndim == 2:
			# Plots for 2D dataset. Currently not used.
			if PlotFilaments:
				FilPositions = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection(line_segments)
				plt.xlabel('$\mathregular{x}$')
				plt.ylabel('$\mathregular{y}$')
				plt.title('Positions of the filaments with '+ self.nPart_text+ '$^3$ particles')
			if PlotFilamentsWCritPts:
				FilPositions_WCritPts = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection(line_segments)
				plt.hold(True)
				plt.plot(self.CritPointXpos, self.CritPointYpos, 'ro', alpha=0.7, markersize=3)
				#plt.plot(self.CritPointXposNOTCON, self.CritPointYposNOTCON, 'go', alpha=0.7, markersize=3)
				plt.xlabel('$\mathregular{x}$')
				plt.ylabel('$\mathregular{y}$')
				plt.title('Position of the filaments with critical points shown')
				plt.hold(False)
		if ndim == 3:
			# Plots of 3D datasets.
			if PlotFilaments:
				# Plots the filaments in a 3D box
				FilPositions = plt.figure()
				ax = FilPositions.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.rainbow)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				ax.set_xlabel('$\mathregular{x}$' + LegendText)
				ax.set_ylabel('$\mathregular{y}$' + LegendText)
				ax.set_zlabel('$\mathregular{z}$' + LegendText)
				plt.title('3D Position of the filaments with '+ self.nPart_text+ '$\mathregular{^3}$ particles.')
			if PlotFilamentsWCritPts:
				# Plots the filaments, including the critical points, in a 3D box
				FilPositions_WCritPts = plt.figure()
				ax = FilPositions_WCritPts.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				plt.hold(True)
				ax.plot(self.CritPointXpos, self.CritPointYpos, self.CritPointZpos, 'ro', alpha=0.7, markersize=3)
				ax.set_xlabel('$\mathregular{x}$' + LegendText)
				ax.set_ylabel('$\mathregular{y}$' + LegendText)
				ax.set_zlabel('$\mathregular{z}$' + LegendText)
				plt.title('3D Position of the filaments with critical points')
				plt.hold(False)
			if Projection2D:
				# Creates a 2D projection onto the X-Y plane of the 3D box
				FilPositions_2DProjection = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, array=ColorArray, cmap=plt.cm.rainbow)
				ax.add_collection(line_segments)
				ax.set_xlabel('$\mathregular{x}$' + LegendText)
				ax.set_ylabel('$\mathregular{y}$' + LegendText)
				plt.title('2D Position projection of the filaments')
				if ColorBarZDir:
					# 2D projection where the colorbar shows the value of the -average- z-value
					FilPositions_2DProjectionColorBarZDir = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.FilamentPos, array=ColorMap2D, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DProjectionColorBarZDir.colorbar(line_segmentsCbar)
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D Position projection of the filaments.\n Color based on average z-position')
				if ColorBarLength:
					# 2D projection where the colorbar shows the length of the filament
					FilPositions_2DProjectionColobarLength = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbarLen = LineCollection(self.FilamentPos, array=ColorMapLength, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbarLen)
					FilPositions_2DProjectionColobarLength.colorbar(line_segmentsCbarLen)
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D Position projection of the filaments.\n Color based on the lenght of the filament')

			if IncludeSlicing and PlotFilaments:
				# 3D box where filaments are masked over a slice of the box.
				FilamentSliced = plt.figure()
				ax = FilamentSliced.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.MaskedFilamentSegments, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimMasked, zdir='z')
				ax.set_xlabel('$\mathregular{x}$' + LegendText)
				ax.set_ylabel('$\mathregular{y}$' + LegendText)
				ax.set_zlabel('$\mathregular{z}$' + LegendText)
				plt.title('Sliced segment of the 3D box')
				# 3D box where filaments are masked. Filaments will be 'cut off' if it goes outside the masking boundary.
				FilamentCutOff = plt.figure()
				ax2 = FilamentCutOff.gca(projection='3d')
				ax2.set_xlim(self.xmin, self.xmax)
				ax2.set_ylim(self.ymin, self.ymax)
				ax2.set_zlim(self.zmin, self.zmax)
				line_segments_CO = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax2.add_collection3d(line_segments_CO, self.CutOffzDim, zdir='z')
				ax2.set_xlabel('$\mathregular{x}$' + LegendText)
				ax2.set_ylabel('$\mathregular{y}$' + LegendText)
				ax2.set_zlabel('$\mathregular{z}$' + LegendText)
				plt.title('Sliced segment of the 3D box, with filaments cut off outside of boundary')
				
				# 2D projection onto the X-Y plane with masked filaments. Segment in z-dir cut off outside of boundary. Colorbar based on filament sigma value
				FilPositions_2DProjection_sigmacolorbar, ax_sigmacbar = plt.subplots()
				ax_sigmacbar.set_xlim(self.xmin, self.ymax)
				ax_sigmacbar.set_ylim(self.ymin, self.ymax)
				line_segments_sigma = LineCollection(self.CutOffFilamentSegments, array=self.Filament_sigma)
				ax_sigmacbar.add_collection(line_segments_sigma)
				CbarSigma =FilPositions_2DProjection_sigmacolorbar.colorbar(line_segments_sigma)
				CbarSigma.set_clim(vmin=2, vmax=8)
				ax_sigmacbar.set_xlabel('$\mathregular{x}$' + LegendText)
				ax_sigmacbar.set_ylabel('$\mathregular{y}$' + LegendText)
				plt.title('2D Position projection sliced box. \n Colorbar based on filament sigma value')
				
				if Projection2D:
					# 2D projection onto the X-Y plane with masked filaments
					FilPositions_2DProjectionSliced = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segments2 = LineCollection(self.MaskedFilamentSegments, array=ColorArray, cmap=plt.cm.gist_ncar)
					ax.add_collection(line_segments2)
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D Position projection sliced box')

				if ColorBarZDir:
					# 2D projection onto the X-Y plane with masked filaments. Colorbar based on -average- z-value
					FilPositions_2DSlicedProjectionColorBarZDir = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.MaskedFilamentSegments, array=ColorMap2DMasked, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DSlicedProjectionColorBarZDir.colorbar(line_segmentsCbar)
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D projection of the filaments in a sliced segment of the box.\n Color based on average z-position')

				if ColorBarLength:
					# 2D projection onto the X-Y plane with masked filaments. Colorbar based on filament length
					FilPositions_2DSlicedProjectionColorBarLen = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.MaskedFilamentSegments, array=ColorMapLengthMasked, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DSlicedProjectionColorBarLen.colorbar(line_segmentsCbar)
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D projection of the filaments in a sliced segment of the box.\n Color based on filament length')
					# Same as above, but now including masked critical point over the same masked region
					FilPositions_2DSlicedProjection_Cbar_CritPts = plt.figure()
					plt.hold(True)
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.MaskedFilamentSegments, array=ColorMapLengthMasked, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DSlicedProjection_Cbar_CritPts.colorbar(line_segmentsCbar)
					plt.plot(self.MaskedXCP, self.MaskedYCP, 'ro')
					ax.set_xlabel('$\mathregular{x}$' + LegendText)
					ax.set_ylabel('$\mathregular{y}$' + LegendText)
					plt.title('2D projection of the filaments in a sliced segment of the box.\n Color based on filament length. Includes Masked Crit points.')

			if IncludeDMParticles:
				# Plots a 2D histogram of the masked dark matter particles
				self.PartPosX = PartPosX
				self.PartPosY = PartPosY
				self.PartPosZ = PartPosZ
				DMParticleHist = plt.figure()
				plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
				plt.xlabel('$\mathregular{x}$' + LegendText)
				plt.ylabel('$\mathregular{y}$' + LegendText)
				plt.title('Dark matter density field over a segment of the particle box.')

				# One dimensional histogram of x-position of the masked dark matter particles. Only for testing purposes
				ONEDHistX = plt.figure()
				plt.hist(self.PartPosX, bins=100)
				plt.xlabel('$\mathregular{x}$' + LegendText)

				# 2D histogram of masked dark matter particles overplotted with masked filaments in the X-Y plane.
				DMParticleHistwFilaments = plt.figure()
				plt.hold(True)
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segmentsDM = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMap2DCutOff, cmap=plt.cm.rainbow)
				ax.add_collection(line_segmentsDM)
				DMParticleHistwFilaments.colorbar(line_segmentsDM)
				plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
				plt.xlabel('$\mathregular{x}$' + LegendText)
				plt.ylabel('$\mathregular{y}$' + LegendText)
				plt.hold("off")
				plt.title('Dark matter density field over a segment of the particle box. \n Includes filaments with colorbar.'\
						+ 'Colors indicate average z-value.')
				plt.hold(False)

				# 2D histogram of masked dark matter particles overplotted with masked filaments in the X-Y plane. Colorbar based on filament length.
				DMParticleHistwFilamentsLengthCbar = plt.figure()
				plt.hold(True)
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segmentsDMlen = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMapLengthMasked, cmap=plt.cm.rainbow)
				ax.add_collection(line_segmentsDMlen)
				DMParticleHistwFilamentsLengthCbar.colorbar(line_segmentsDMlen)
				plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
				plt.xlabel('$\mathregular{x}$' + LegendText)
				plt.ylabel('$\mathregular{y}$' + LegendText)
				plt.hold(False)
				plt.title('Dark matter density field over a segment of the particle box. \n Includes filaments with colorbar.'\
						 +'Colors indicate length of a filament.')

				# Values in griddata = mass. Currently normalized to 1.
				# Compute mass as 0.23*rho_{crit,0}*Volume_box/Num_particles
				# See discussion with Max
				"""
				Mpc = 3.08568025e22
				G_grav = 6.67258e-11
				H_0 = 0.7*100*1e3/Mpc
				mass = 0.23*(3*H_0**2/(8*np.pi*G_grav))*(256.0*Mpc/0.7)**3/(512.0)**3griddata')
				"""
			if parsed_arguments.bwMethod and IncludeDMParticles:
				if parsed_arguments.bwMethod != 'Scott':
					if len(parsed_arguments.bwMethod) == 2:
						Column = 2
						Row = 1
					elif len(parsed_arguments.bwMethod) > 2 and len(parsed_arguments.bwMethod) <= 4:
						Column = 2
						Row = 2
					elif len(parsed_arguments.bwMethod) > 4 and len(parsed_arguments.bwMethod) <= 6:
						Column = 3
						Row = 2
					elif len(parsed_arguments.bwMethod) > 6 and len(parsed_arguments.bwMethod) <= 9:
						Column = 3
						Row = 3
				
				# Using gaussian kernel to plot density field of DM particle positions
				DMParticles_kernelPlot, ax_kernel = plt.subplots()
				if len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] == 'Scott':
					cax = ax_kernel.imshow(np.rot90(self.Interpolated_Z), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel.set_title('Bandwidth = Scott')
					cbar = DMParticles_kernelPlot.colorbar(cax)
				elif len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] != 'Scott':
					cax = ax_kernel.imshow(np.rot90(self.Interpolated_Z), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel.set_title('Bandwidth = ' + parsed_arguments.bwMethod[0]) 
					cbar = DMParticles_kernelPlot.colorbar(cax)
				else:
					for j in range(1,len(self.Interpolated_Z)+1):
						plt.subplot(Column, Row, j)
						cax = plt.imshow(np.rot90(self.Interpolated_Z[j-1][0]), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
						cbar = DMParticles_kernelPlot.colorbar(cax)
						plt.title('Bandwidth = ' + parsed_arguments.bwMethod[j-1])
				ax_kernel.set_xlim([self.xmin, self.xmax])
				ax_kernel.set_ylim([self.ymin, self.ymax])
				plt.xlabel('$\mathregular{x}$' + LegendText)
				plt.ylabel('$\mathregular{y}$' + LegendText)
				plt.tight_layout()

				# Plotting logarithmic value of the density
				DMParticles_kernelPlot_logarithmic, ax_kernel_log = plt.subplots()
				if len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] == 'Scott':
					cax_log = ax_kernel_log.imshow(np.rot90(self.Logarithmic_density), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel_log.set_title('Bandwidth = Scott. Logarithmic')
					cbar_log = DMParticles_kernelPlot_logarithmic.colorbar(cax_log)
				elif len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] != 'Scott':
					cax_log = ax_kernel_log.imshow(np.rot90(self.Logarithmic_density), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel_log.set_title('Bandwidth = ' + parsed_arguments.bwMethod[0] + '. Logarithmic') 
					cbar_log = DMParticles_kernelPlot_logarithmic.colorbar(cax_log)
				else:
					for j in range(1,len(self.Logarithmic_density)+1):
						plt.subplot(Column, Row, j)
						cax_log = plt.imshow(np.rot90(self.Logarithmic_density[j-1][0]), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
						cbar_log = DMParticles_kernelPlot_logarithmic.colorbar(cax_log)
						plt.title('Bandwidth = ' + parsed_arguments.bwMethod[j-1] + '. Logarithmic')
				ax_kernel_log.set_xlim([self.xmin, self.xmax])
				ax_kernel_log.set_ylim([self.ymin, self.ymax])
				plt.xlabel('$\mathregular{x}$' + LegendText)
				plt.ylabel('$\mathregular{y}$' + LegendText)
				plt.tight_layout()
					
				# Overplot with filaments on the gaussian kernel plots
				# Only done if input parameters for bw_method is none or only one scalar
				if parsed_arguments.bwMethod[0] == 'Scott' or len(parsed_arguments.bwMethod) == 1:
					# Filaments overplotted on the density field
					DMParticles_kernelPlot_wFilaments, ax_kernel_wfil = plt.subplots()
					ax_kernel_wfil.imshow(np.rot90(self.Interpolated_Z), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel_wfil.set_title('Bandwidth = ' + parsed_arguments.bwMethod[0] +'\n' 
											+ self.nPart_text + ' particle subsample. ' + self.Alternative_sigmaTitle) 
					ax_kernel_wfil.set_xlim([self.xmin, self.xmax])
					ax_kernel_wfil.set_ylim([self.ymin, self.ymax])
					line_segmentsDMlen_kernel = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMapLengthMasked, cmap=plt.cm.rainbow)
					ax_kernel_wfil.add_collection(line_segmentsDMlen_kernel)
					Cbar_kernelwFil = DMParticles_kernelPlot_wFilaments.colorbar(line_segmentsDMlen_kernel)
					Cbar_kernelwFil.set_clim(vmin=0, vmax=max(ColorMapLengthMasked))
					plt.xlabel('$\mathregular{x}$' + LegendText)
					plt.ylabel('$\mathregular{y}$' + LegendText)

					# Filaments overplotted on the logarithmic scaled density field
					DMParticles_kernelPlot_wFilaments_log, ax_kernel_wfil_log = plt.subplots()
					ax_kernel_wfil_log.imshow(np.rot90(self.Logarithmic_density), extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					ax_kernel_wfil_log.set_title('Bandwidth = ' + parsed_arguments.bwMethod[0] + \
										'. Logarithmic density. \n' + self.nPart_text + ' particle subsample. ' + self.Alternative_sigmaTitle) 
					ax_kernel_wfil_log.set_xlim([self.xmin, self.xmax])
					ax_kernel_wfil_log.set_ylim([self.ymin, self.ymax])
					line_segmentsDMlen_kernel_log = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMapLengthMasked, cmap=plt.cm.rainbow)
					ax_kernel_wfil_log.add_collection(line_segmentsDMlen_kernel_log)
					Cbar_kernelwFil_log = DMParticles_kernelPlot_wFilaments_log.colorbar(line_segmentsDMlen_kernel_log)
					Cbar_kernelwFil_log.set_clim(vmin=0, vmax=max(ColorMapLengthMasked))
					plt.xlabel('$\mathregular{x}$' + LegendText)
					plt.ylabel('$\mathregular{y}$' + LegendText)

					# Different zoom-ins on the density field with filaments
					DMParticles_kernelPlot_wFilaments_log_Zoomed, ax_kernel_wfil_log_zoom = plt.subplots(figsize=(8.5,8.5))
					DMParticles_kernelPlot_wFilaments_log_Zoomed.suptitle(
							'Zoomed in segments of the (logarithmic) density field with filaments \n Colorbar based on filament length'
							+'\n' + self.nPart_text + 'particle subsample. ' + self.Alternative_sigmaTitle)
					plt.subplot(2,2,1)
					plt.imshow(np.rot90(self.Log_zoomed_density[0]), extent=[70, 180, 140, 250])#extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					line_segmentsDMlen_kernel_log_zoom = LineCollection(
						self.CutOffFilamentSegments, linestyle='solid',
						array=ColorMapLengthMasked, cmap=plt.cm.autumn)
					plt.gca().add_collection(line_segmentsDMlen_kernel_log_zoom)
					plt.xlim(70, 180)
					plt.ylim(140,250)

					plt.subplot(2,2,2)
					plt.imshow(np.rot90(self.Log_zoomed_density[1]), extent=[100,150,175,225]) #extent=[self.xmin, self.xmax, self.ymin, self.ymax])
					line_segmentsDMlen_kernel_log_zoom = LineCollection(
						self.CutOffFilamentSegments, linestyle='solid',
						array=ColorMapLengthMasked, cmap=plt.cm.autumn)
					plt.gca().add_collection(line_segmentsDMlen_kernel_log_zoom)
					plt.xlim(100, 150)
					plt.ylim(175,225)
					
					plt.subplot(2,2,3)
					plt.imshow(np.rot90(self.Log_zoomed_density[2]), extent=[120, 145, 190, 215])
					line_segmentsDMlen_kernel_log_zoom = LineCollection(
						self.CutOffFilamentSegments, linestyle='solid',
						array=ColorMapLengthMasked, cmap=plt.cm.autumn)
					plt.gca().add_collection(line_segmentsDMlen_kernel_log_zoom)
					plt.xlim(120, 145)
					plt.ylim(190,215)
					
					plt.subplot(2,2,4)
					plt.imshow(np.rot90(self.Log_zoomed_density[2]), extent=[115, 125, 195, 210])
					line_segmentsDMlen_kernel_log_zoom = LineCollection(
						self.CutOffFilamentSegments, linestyle='solid',
						array=ColorMapLengthMasked, cmap=plt.cm.autumn)
					plt.gca().add_collection(line_segmentsDMlen_kernel_log_zoom)
					plt.xlim(115, 125)
					plt.ylim(195,210)

					DMParticles_kernelPlot_wFilaments_log_Zoomed.text(0.5, 0.04, '$\mathregular{x}$' + LegendText, ha='center', fontsize=16)
					DMParticles_kernelPlot_wFilaments_log_Zoomed.text(0.04, 0.5, '$\mathregular{y}$' + LegendText, va='center', rotation='vertical', fontsize=16)
					# Colorbar at the bottom of the figure. [left, bottom, width, height], % of box
					cax = DMParticles_kernelPlot_wFilaments_log_Zoomed.add_axes([0.92, 0.08, 0.03, 0.8])
					DMParticles_kernelPlot_wFilaments_log_Zoomed.colorbar(line_segmentsDMlen_kernel_log_zoom, cax=cax)#, orientation='horizontal')
					#plt.tight_layout()
					
		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			if HistogramPlots:
				savefigure(ConnectedHist, 'NumberFilamentConnectedHistogram')
				savefigure(FilamentLengthsHist, 'FilamentLengthsHistogram')
				savefigure(FilamentPtsHis, 'FilamentPointsHistogram')
				#if IncludeDMParticles:
				#	savefigure(NumParticlesFilamentHist, 'NumParticlesPerFilament')
			if PlotFilaments:
				savefigure(FilPositions, 'FilamentPositions')
			if PlotFilamentsWCritPts:
				savefigure(FilPositions_WCritPts, 'FilamentPositionsWithCritPts')
			if Projection2D:
				savefigure(FilPositions_2DProjection, '2DProjectedFilamentPosition')
				if ColorBarZDir:
					savefigure(FilPositions_2DProjectionColorBarZDir, '2DProjectionColorBarZDir')
				if ColorBarLength:
					savefigure(FilPositions_2DProjectionColobarLength, '2DProjectionColobarLength')
			if IncludeSlicing and PlotFilaments:
				savefigure(FilamentSliced, 'Sliced3dBox')
				savefigure(FilamentCutOff, 'CutOffFilaments')
				savefigure(FilPositions_2DProjection_sigmacolorbar, 'CutoffFilaments_sigmacolorbar')
				if Projection2D:
					savefigure(FilPositions_2DProjectionSliced, '2DProjectionSliced3dBox')
				if ColorBarZDir:
					savefigure(FilPositions_2DSlicedProjectionColorBarZDir, '2DProjectionSlicedColobarZDir')
				if ColorBarLength:
					savefigure(FilPositions_2DSlicedProjectionColorBarLen, '2DProjectionSlicedColobarLength')
					savefigure(FilPositions_2DSlicedProjection_Cbar_CritPts, '2DProjectionSliced_CbarLength_CritPts')
			if IncludeDMParticles:
				Masked_filename = 'Masked' + 'X'*MaskXdir + 'Y'*MaskYdir + 'Z'*MaskZdir
				savefigure(DMParticleHist, 'DMParticleHistogram_' + Masked_filename)
				savefigure(DMParticleHistwFilaments, 'DMParticleHistogramWFIlaments' + Masked_filename)
				savefigure(DMParticleHistwFilamentsLengthCbar, 'DMParticleHistogramWFilaments_LengthCbar_' + Masked_filename)
				savefigure(ONEDHistX, 'DMParticle1DHistogramXposition')
			if parsed_arguments.bwMethod and IncludeDMParticles:
				if len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] == 'Scott':
					savefigure(DMParticles_kernelPlot, 'DMParticles_kernelPlot_Scott')
					savefigure(DMParticles_kernelPlot_logarithmic, 'DMParticles_kernelPlot_Scott_logarithmic')
					savefigure(DMParticles_kernelPlot_wFilaments,  'DMParticles_kernelPlot_Scott_wFilaments')
					savefigure(DMParticles_kernelPlot_wFilaments_log, 'DMParticles_kernelPlot_Scott_wFilaments_logarithmic')
				elif len(parsed_arguments.bwMethod) == 1 and parsed_arguments.bwMethod[0] != 'Scott':
					savefigure(DMParticles_kernelPlot, 'DMParticles_kernelPlot_Scott' + parsed_arguments.bwMethod[0])
					savefigure(DMParticles_kernelPlot_logarithmic, 'DMParticles_kernelPlot_Scott_logarithmic' + parsed_arguments.bwMethod[0])
					savefigure(DMParticles_kernelPlot_wFilaments,  'DMParticles_kernelPlot_Scott_wFilaments' + parsed_arguments.bwMethod[0])
					savefigure(DMParticles_kernelPlot_wFilaments_log, 
							'DMParticles_kernelPlot_Scott_wFilaments_logarithmic' + parsed_arguments.bwMethod[0])
					savefigure(DMParticles_kernelPlot_wFilaments_log_Zoomed, 
							'DMParticles_kernelPlot_wFilaments_logarithmic_Zoomed' + parsed_arguments.bwMethod[0])
				else:
					savefigure(DMParticles_kernelPlot, 'DMParticles_kernelPlot_subplots')
					savefigure(DMParticles_kernelPlot_logarithmic, 'DMParticles_kernelPlot_subplots_logarithmic')
		else:
			plt.show()

		plt.close('all')

	def Unpack_filament_data(self, Box_info, CP_coord, Fil_coord, CP_data, Fil_data):
		""" Unpacks filament data from the read filament data module. """
		self.xmin = Box_info[0]
		self.xmax = Box_info[1]
		self.ymin = Box_info[2]
		self.ymax = Box_info[3]
		self.zmin = Box_info[4]
		self.zmax = Box_info[5]

		self.CritPointXpos = CP_coord[0]
		self.CritPointYpos = CP_coord[1]
		self.CritPointZpos = CP_coord[2]
		self.CP_type = CP_coord[3]
		self.CP_persistent_pair = CP_coord[4] 
		self.Critpts_filamentID = CP_coord[5] 
		self.CP_id_of_connecting_filament = CP_coord[6] 
		self.Number_filaments_connecting_to_CP = CP_coord[7]

		self.NFils = Fil_coord[0]
		self.FilamentPos = Fil_coord[1]
		self.xdimPos = Fil_coord[2]
		self.ydimPos = Fil_coord[3]
		self.zdimPos = Fil_coord[4]
		self.NFilamentPoints = Fil_coord[5]
		self.FilID = Fil_coord[6]
		self.PairIDS = Fil_coord[7]

		self.Persistence_ratio = CP_data[0]
		self.Persistence_nsigmas = CP_data[1]
		self.Persistence = CP_data[2]
		self.Persistence_pair = CP_data[3]
		self.Parent_index = CP_data[4]
		self.Parent_log_index = CP_data[5]
		self.log_CritPoint_field_value = CP_data[6]
		self.CritPoint_field_value = CP_data[7]
		self.CritPoint_cell = CP_data[8]
		
		self.Filament_field_value = Fil_data[0]
		self.orientation = Fil_data[1]
		self.Filament_cell = Fil_data[2]
		self.log_Filament_field_value = Fil_data[3]
		self.Filament_type = Fil_data[4]
		self.Filament_sigma = Fil_data[5]
		
	def Solve(self, filename, ndim=3, Sigma_threshold = False, robustness=False):
		""" 
		Runs the whole thing.
		Creates a pickle file of certain data unless it already exist.
		Removes the old pickle files if any of the mask directions are changed.
		"""
		# Determines folder name of pickle files
		cachedir_foldername_extra = self.model + 'npart'+str(self.nPart)
		if self.SigmaArg:
			cachedir_foldername_extra += 'nsig'+str(self.SigmaArg)
		else:
			cachedir_foldername_extra += 'nsig3'

		if robustness:
			cachedir_foldername_extra += 'TRIM'

		if HOMEPC == 0:
			cachedir='/PythonCaches/Disperse_analysis/'+cachedir_foldername_extra+'/'
		else:
			cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/' + cachedir_foldername_extra + '/'
		if not os.path.isdir(cachedir):
			os.makedirs(cachedir)

		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'Conditions_check.p'
		Disperse_data_check_fn = cachedir + 'Disperse_data.p'

		Pickle_check_list = [Mask_direction_check, Mask_boundary_list, Sigma_threshold]
		if os.path.isfile(Pickle_check_fn):
			"""
			Check if pickle file exist. If it exists, it will read from the existing pickle file and get data from there.
			The pickle file will check if the masking conditions are changed or not. If it is changed, will recalculate stuff.
			If pickle file does not exist, then the program will calculate stuff. 
			"""
			Pickle_check = pickle.load(open(Pickle_check_fn, 'rb'))
			if Pickle_check != Pickle_check_list:
				print 'Some global parameters are changed. Recomputing everything and creating new pickle files'
				os.remove(Boundary_check_cachefn)
				os.remove(Mask_slice_cachefn)
				pickle.dump(Pickle_check_list, open(Pickle_check_fn, 'wb'))
		else:
			pickle.dump(Pickle_check_list, open(Pickle_check_fn, 'wb'))

		if os.path.isfile(Disperse_data_check_fn):
			## Read filament data from the skeleton file
			print "reading from filament data pickle file..."
			Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = pickle.load(open(Disperse_data_check_fn, 'rb'))
		else:
			Read_filament_data_instance = ReadFilamentData.read_disperse_output(directory=file_directory, UnitConverter=UnitConverter)
			Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data = Read_filament_data_instance.get_data(filename=filename)
			pickle.dump(
				[Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data],
				open(Disperse_data_check_fn, 'wb')
				)
		self.Unpack_filament_data(Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data)
		
		if os.path.isfile(Boundary_check_cachefn):
			## Computing filaments crossing boundaries
			print "reading from boundary check pickle file..."
			BC_instance_variables = pickle.load(open(Boundary_check_cachefn, 'rb'))
		else:
			BC_instance = BoundaryChecker.BoundaryChecker(
					self.xmin, self.xmax, self.xdimPos,
					self.ydimPos, self.zdimPos, self.FilID, self.NFils
					)
			BC_instance_variables = BC_instance.Get_periodic_boundary()
			pickle.dump(BC_instance_variables, open(Boundary_check_cachefn, 'wb'))

		self.FilamentIDs, self.FilamentPos, self.xdimPos, self.ydimPos, self.zdimPos, self.LengthSplitFilament, self.FilLengths = BC_instance_variables

		if IncludeSlicing:
			## Masking filament and critical points
			if os.path.isfile(Mask_slice_cachefn):
				print "reading from mask_slice pickle file..."
				Mask_instance_variables, Masked_critpts = pickle.load(open(Mask_slice_cachefn, 'rb'))
			else:
				Mask_instance = FilamentMasking.FilamentMasking(
						self.FilamentPos, self.xdimPos, self.ydimPos, self.zdimPos,
						self.FilLengths , self.NFils, Mask_direction_check, Mask_boundary_list
						)
				Mask_instance_variables = Mask_instance.Mask_slices()
				Masked_critpts = MaskCritPts.Mask_CPs(
							self.CritPointXpos, self.CritPointYpos, self.CritPointZpos,
							Mask_boundary_list, Mask_direction_check
							)
				pickle.dump([Mask_instance_variables, Masked_critpts], open(Mask_slice_cachefn, 'wb'))

			self.MaskedFilamentSegments = Mask_instance_variables[0]
			self.MaskedLengths = Mask_instance_variables[1]
			self.zdimMasked = Mask_instance_variables[2]
			self.CutOffFilamentSegments = Mask_instance_variables[3]
			self.CutOffLengths = Mask_instance_variables[4]
			self.CutOffzDim = Mask_instance_variables[5]
			self.MaskedXCP, self.MaskedYCP, self.MaskedZCP = Masked_critpts

		# Computing other properties
		self.Number_filament_connections()
		if Sigma_threshold:
			# Filters filament based on a given sigma threshold. Still work in progress
			self.Filter_filaments2(Sigma_threshold)
		
		if HOMEPC == 1 and IncludeDMParticles and parsed_arguments.bwMethod:
			# Pickle file for interpolated densities using Gaussian KDE
			Density_dir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/Interpolated_Densities/'
			if not os.path.isdir(Density_dir):
				os.makedirs(Density_dir)

			Masked_density_dir = Density_dir
			if MaskXdir:
				Masked_density_dir += 'Xmin' + str(LowerBoundaryXDir/UnitConverter) + 'max' + str(UpperBoundaryXDir/UnitConverter) + '/'
			if MaskYdir:
				Masked_density_dir += 'Ymin' + str(LowerBoundaryYDir/UnitConverter) + 'max' + str(UpperBoundaryYDir/UnitConverter) + '/'
			if MaskZdir:
				Masked_density_dir += 'Zmin' + str(LowerBoundaryZDir/UnitConverter) + 'max' + str(UpperBoundaryZDir/UnitConverter) + '/'

			if not os.path.isdir(Masked_density_dir):
				os.makedirs(Masked_density_dir)
			
			Particle_positions = [PartPosX, PartPosY, PartPosZ]
			Box = [self.xmin, self.xmax, self.ymin, self.ymax]
			self.Zoom_areas = [[70, 180, 140, 250], [100, 150, 170, 230], [120, 145, 190, 215], [115, 125, 195, 210]]
			if len(parsed_arguments.bwMethod) == 1:
				# If there is only one argument
				Interpolated_density_cachefn = Masked_density_dir + 'InterpolatedDensities_bandwidth_' + parsed_arguments.bwMethod[0] + '.p'
				Interpolated_density_zoomed_cachefn = Masked_density_dir + 'InterpolatedDensities_Zoomed_bandwidth_' + parsed_arguments.bwMethod[0] + '.p'
				Interpolation_instance = InterpolateDensity.InterpolateDensity(parsed_arguments.bwMethod[0], Particle_positions, Box, npoints=200j)
				if os.path.isfile(Interpolated_density_cachefn):
					print "reading from interpolated density pickle file, with bandwidth = " + parsed_arguments.bwMethod[0] + "..."
					self.Interpolated_Z, self.Logarithmic_density = pickle.load(open(Interpolated_density_cachefn, 'rb'))
				else:
					self.Interpolated_Z, self.Logarithmic_density = Interpolation_instance.Interpolate_DM_particles()
					pickle.dump([self.Interpolated_Z, self.Logarithmic_density],
						open(Interpolated_density_cachefn ,'wb'))

				if os.path.isfile(Interpolated_density_zoomed_cachefn):
					self.Zoomed_density, self.Log_zoomed_density = pickle.load(open(Interpolated_density_zoomed_cachefn, 'rb'))
				else:
					self.Zoomed_density, self.Log_zoomed_density = Interpolation_instance.compute_zoomed_density(self.Zoom_areas)
					pickle.dump([self.Zoomed_density, self.Log_zoomed_density],
					 	open(Interpolated_density_zoomed_cachefn, 'wb'))

			elif len(parsed_arguments.bwMethod) > 1:
				# If there are multiple arguments
				self.Interpolated_Z = []
				self.Logarithmic_density = []
				for bandwidths in parsed_arguments.bwMethod:
					Interpolated_density_cachefn = Masked_density_dir + "InterpolatedDensities_bandwidth_" + bandwidths + '.p'
					Interpolation_instance = InterpolateDensity.InterpolateDensity(bandwidths, Particle_positions, Box, npoints=200j)
					if os.path.isfile(Interpolated_density_cachefn):
						print "reading from interpolated density pickle file, with bandwidth = " + bandwidths  + "..."
						Density, Log_density = pickle.load(open(Interpolated_density_cachefn, 'rb'))
						self.Interpolated_Z.append(Density)
						self.Logarithmic_density.append(Log_density)
					else:
						Density, Log_density = Interpolation_instance.Interpolate_DM_particles()
						pickle.dump([Density, Log_density], open(Interpolated_density_cachefn ,'wb'))
						self.Interpolated_Z.append(Density)
						self.Logarithmic_density.append(Log_density)
			
			#Particle_positions = [PartPosX, PartPosY, PartPosZ]
			#Box = [self.xmin, self.xmax, self.ymin, self.ymax]
			#Interpolation_instance = InterpolateDensity.InterpolateDensity(parsed_arguments.bwMethod[0], Particle_positions, Box, npoints=20j)
			#self.Interpolated_Z, self.Logarithmic_density = Interpolation_instance.Interpolate_DM_particles()
			#self.Zoomed_density, self.Log_zoomed_density = Interpolation_instance.compute_zoomed_density([[70, 180, 140, 250], [100, 150, 175, 225], [120,145,190,215], [115, 125, 195, 210]])

		if not Comparison:
			if self.savefile == 2:
				print 'Done! No figures saved.'
			else:
				self.Plot_Figures(filename, ndim)
		
		self.SolveRun = True
		return self.NumFilamentConnections, sorted(self.FilLengths), self.NFilamentPoints

	def get_3D_pos(self):
		""" 
		Create 3D filament position. Must have called solve() function first. 
		Also returns filament IDs.
		"""
		if not self.SolveRun:
			raise ValueError('You have not called the solve() function yet! Cannot give 3D filament position.')
		Filament_3DPos = []
		for i in range(len(self.xdimPos)):
			Filament_3DPos.append(np.column_stack((self.xdimPos[i], self.ydimPos[i], self.zdimPos[i])))
		return np.array(Filament_3DPos), self.FilamentIDs

class Read_solve_files():
	def __init__(self):
		self.read_solvefile()
		self.Create_Mask()
		#self.Create_KDTree()

	def read_solvefile(self):
		""" 
		Reads the .solve file which contains dark matter particle positions and velocities
		Saves all the particles which is later masked
		"""
		time_start = time.clock()
		print 'Reading data for the file: ', solve_file_dir, solve_filename, '. May take a while...'
		self.PartPosX = []
		self.PartPosY = []
		self.PartPosZ = []
		
		PartVelX = []
		PartVelY = []
		PartVelZ = []

		self.ParticlePos = []
		with open(os.path.join(file_directory+solve_file_dir, solve_filename), 'r') as datafiles:
			next(datafiles)
			for line in datafiles:
				data_set = line.split()
				self.ParticlePos.append(np.array([float(data_set[0]),float(data_set[1]),float(data_set[2])])*UnitConverter)
				#self.PartPosX.append(float(data_set[0])*UnitConverter)
				#self.PartPosY.append(float(data_set[1])*UnitConverter)
				#self.PartPosZ.append(float(data_set[2])*UnitConverter)
				
			
		#self.PartPosX = np.asarray(self.PartPosX)
		#self.PartPosY = np.asarray(self.PartPosY)
		#self.PartPosZ = np.asarray(self.PartPosZ)
		self.ParticlePos = np.asarray(self.ParticlePos)
		#mask = np.logical_and(self.ParticlePos[:,2] < UpperBoundaryZDir, self.ParticlePos[:,2] > LowerBoundaryZDir)
		print 'Read solve file time: ', time.clock() - time_start, 's'

	def Create_Mask(self):
		""" Creates a mask for the dark matter particles """
		if not MaskXdir and not MaskYdir and MaskZdir:
			self.mask = np.logical_and(self.ParticlePos[:,2] < UpperBoundaryZDir, self.ParticlePos[:,2] > LowerBoundaryZDir)
		elif not MaskXdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(self.ParticlePos[:,0] < UpperBoundaryXDir, self.ParticlePos[:,0] > LowerBoundaryXDir)
		elif not MaskYdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(self.ParticlePos[:,1] < UpperBoundaryYDir, self.ParticlePos[:,1] > LowerBoundaryYDir)
		elif MaskXdir and not MaskYdir and MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.ParticlePos[:,2] < UpperBoundaryZDir, self.ParticlePos[:,2] > LowerBoundaryZDir),\
									   np.logical_and(self.ParticlePos[:,0] < UpperBoundaryXDir, self.ParticlePos[:,0] > LowerBoundaryXDir))
		elif MaskXdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.ParticlePos[:,1] < UpperBoundaryYDir, self.ParticlePos[:,1] > LowerBoundaryYDir),\
									   np.logical_and(self.ParticlePos[:,0] < UpperBoundaryXDir, self.ParticlePos[:,0] > LowerBoundaryXDir))
		elif not MaskXdir and MaskYdir and MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.ParticlePos[:,2] < UpperBoundaryZDir, self.ParticlePos[:,2] > LowerBoundaryZDir),\
									   np.logical_and(self.ParticlePos[:,1] < UpperBoundaryYDir, self.ParticlePos[:,1] > LowerBoundaryYDir))
		else:
			self.mask = False

	def Create_KDTree(self):
		""" Creates a KDTree of the masked dark matter particles """
		time_start = time.clock()
		if self.mask is not False:
			self.PartPosX = self.ParticlePos[self.mask,0]
			self.PartPosY = self.ParticlePos[self.mask,1]
			self.PartPosZ = self.ParticlePos[self.mask,2]
		else:
			self.PartPosX = self.ParticlePos[:,0]
			self.PartPosY = self.ParticlePos[:,1]
			self.PartPosZ = self.ParticlePos[:,1]

		DM_points = np.dstack((self.PartPosX.ravel(), self.PartPosY.ravel(), self.PartPosZ.ravel()))
		self.DM_tree = spatial.KDTree(DM_points[0])
		print 'Dark matter particle KDTRee creation time: ', time.clock() - time_start, 's'


def Argument_parser():
	""" Parses optional argument when program is run from the command line """
	print 'Run python code with -h argument to see extra optional arguments'

	parser = argparse.ArgumentParser()
	# Optional arguments
	parser.add_argument("-hp", "--HOMEPC", help="Determines if program is run in UiO or laptop. Set 1 if run in UiO. 0 by default", type=int, default=0)
	parser.add_argument("-kernel", "--KernelMethod", help="Selects different kernels for interpolation. Default = gaussian. tp = tophat, "\
						, type=str, default='gaussian')
	parser.add_argument("-bw_m", "--bwMethod", nargs='*', help="Sets bw_method argument of scipy.stats_gaussian_kde. " \
					+  "Input arguments can be float values, may be multiple. If using default scipy settings, set argument Scott")
	parser.add_argument("-Nproc", "--NumProcesses", help="Sets the number of processes for multiprocessing module. Default to 4", type=int, default=4)
	parser.add_argument("-comp", "--Comparisons", help="If set to 1, compares different particle subsamples and/or gravity models. Default to 0", type=int, default=0)
	parser.add_argument("-hpart", "--HigherPart", help="Includes particles of larger subsamples if set to 1."\
					+ "Aimed to include seperate simulations for larger number of particles. Defalult to 0 (only runs 64^3 particles)", type=int, default=0)
	parser.add_argument("-sigcomp", "--SigmaComp", help="Set to 1 to compare simulations of different sigmas. 0 by default.", type=int, default=0)
	parser.add_argument("-modelcomp", "--ModelCompare", help="Set to 1 to compare all the modified gravity models. 0 by default.", type=int, default=0)
	parser.add_argument("-NpartModel", "--NumPartModel", help="Computes number of particles per filament for an input model. Models may be:" \
					+ "lcdm, symmX (X=A,B,C,D) or fofrY (Y=4,5,6). Runs all by default.", type=str)

	# Parse arguments
	args = parser.parse_args()
	return args

def Multiprocess_filament_per_filament(Instance, distance_threshold, Filament, i):
	""" 
	Computes number of particles per filament in parallel.
	Input must be 3D filament position, as well as a given distance threshold from particle to filament.
	"""
	NumParticles = Instance.solve(Filament[i], distance_threshold)
	return NumParticles

def Multiprocess_masked_particle_ids(Instance, Filament, i):
	""" 
	Computes the masked particle indices.
	Input must be 3D filament position.
	"""
	Indices = Instance.get_masked_ids(Filament[i])
	return Indices

def Multiprocess_particle_distances(Instance, Masked_ids, Filament, i):
	""" 
	Computes the masked particle indices.
	Input must be 3D filament position.
	"""
	Distances = Instance.get_distance(Filament[i], Masked_ids[i])
	return Distances

def Save_NumPartPerFil(name, FilPos, FilID, npart, nsig):
	""" 
	Saves number of particle per filament data in pickle file.
	Input is the model name and 3D filament positions.
	Box expand expands the box which covers the given filament. Ensures that it includes all particles.
	Box_expand must be larger than distance_threshold
	"""
	if box_expand < distance_threshold:
		raise ValueError('Box expand must be larger than the distance threshold!')

	if type(name) != str:
		raise ValueError('Wrong type for name. Current type is %s' %type(name))

	toggle = 0
	Model_check = ['lcdm', 'symm_A', 'symm_B', 'symm_C', 'symm_D', 'fofr4', 'fofr5', 'fofr6']
	for models in Model_check:
		if name == models:
			toggle = True
	if not toggle:
		raise ValueError('Model input name %s not correctly set into the number particle per filament saver.' %name)

	cachedir_ppf_distances = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/Distances/'
	cachedir_ppf_ids = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticlesPerFilament/MaskedParticlesIDs/'
	cache_particledata = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/ParticleData/'
	if not os.path.isdir(cachedir_ppf_distances):
		os.makedirs(cachedir_ppf_distances)
	if not os.path.isdir(cachedir_ppf_ids):
		os.makedirs(cachedir_ppf_ids)
		
	cachefile_distances = cachedir_ppf_distances + name + '_' + str(npart) + 'part_nsig' + str(nsig) + '_BoxExpand' + str(box_expand) + '.p'
	cachefile_ids = cachedir_ppf_ids + name +  '_' + str(npart) + 'part_nsig' + str(nsig) + '_BoxExpand' + str(box_expand) + '.p' 
	"""
	if os.path.isfile(cachefile):
		print 'Reading number of particles pickle file for ' + name + '...'
		NumPartPerFil, NumPartPerFil_particleIDs = pickle.load(open(cachefile, 'rb'))
	else:
		print 'Computing number of particles per filament for ' + name + '. May take a while...'
		ppf_instance = ParticlesPerFilament.particles_per_filament(name, Mask_check_list2, BoundaryCheck_list2)
		Fixed_args = partial(Multiprocess_filament_per_filament, ppf_instance, distance_threshold, FilPos)
		Multiproc_output = proc.map(Fixed_args, range(len(FilPos)))
		NumPartPerFil = np.array(Multiproc_output)[:,0]
		NumPartPerFil_particleIDs = np.array(Multiproc_output)[:,1]
		proc.close()
		proc.join()
		pickle.dump([NumPartPerFil, NumPartPerFil_particleIDs], open(cachefile, 'wb'))
	"""
	Nfils_iterate = range(len(FilPos))
	if os.path.isfile(cachefile_ids):
		print 'Reading masked particle IDs for ' + name + '...'
		Toggle_distance_computing = 0
		Masked_ids = pickle.load(open(cachefile_ids, 'rb'))
	else:
		print 'Computing masked particle indices for ' + name + '. May take a while...'
		Toggle_distance_computing = 1
		ppf_instance = ParticlesPerFilament.particles_per_filament(name, Mask_check_list2, BoundaryCheck_list2, box_expand)
		Fixed_args_ids = partial(Multiprocess_masked_particle_ids, ppf_instance, FilPos)
		Masked_ids = np.array(proc.map(Fixed_args_ids, Nfils_iterate))
		pickle.dump(Masked_ids, open(cachefile_ids, 'wb'))
		#proc.close()
		#proc.join()
	Distances = []
	if os.path.isfile(cachefile_distances):
		print 'Reading number of particles pickle file for ' + name + '...'
		Distances = pickle.load(open(cachefile_distances, 'rb'))
	elif not os.path.isfile(cachefile_distances) and Toggle_distance_computing:
		print 'Computing number of particles per filament for ' + name + '. may take a while...'
		#Fixed_args_distances = partial(Multiprocess_particle_distances, ppf_instance, Masked_ids, FilPos)
		#Distances = np.array(proc.map(Fixed_args_distances, Nfils_iterate))
		for ids in Nfils_iterate:
			Distances.append(ppf_instance.get_distance(FilPos[ids], Masked_ids[ids]))
		#proc.close()
		#proc.join()
		pickle.dump(Distances, open(cachefile_distances, 'wb'))
        
	NumPartPerFil = []
	for dist in Distances:
		Accepted_distances = np.where(dist <= distance_threshold)[0]
		NumPartPerFil.append(len(Accepted_distances))
	#NumPartPerFil = len(Accepted_distances)

	cache_model = cache_particledata + name + '_particleIDs.p'
	if os.path.isfile(cache_model):
		Particle_IDs = pickle.load(open(cache_model, 'rb'))
	else:
		raise ValueError("Particle ID pickle file does not exist!")

	Filament_part_IDs = []
	for mask in Masked_ids:
		Filament_part_IDs.append(Particle_IDs[mask])

	Npart_per_fil = []
	#Npart_per_fil_partIDs = []
	#def append_ids(id_list):
	#	temp_id = []
	#	for id_list2 in id_list:
	#		for ids in id_list2:
	#			temp_id.append(ids)
	#	Npart_per_fil_partIDs.append(np.array(temp_id))

	total = 0
	#total_partids = []
	#id_append_check = 0
	for i in range(len(FilID)):	
		if i != len(FilID) - 1:
			if FilID[i+1] == FilID[i]:
				#total_partids.append(NumPartPerFil_particleIDs[i])
				total += NumPartPerFil[i]
				if i == len(FilID) - 2:
					#id_append_check = 1
					Npart_per_fil.append(total)
			else:
				if total != 0:
					Npart_per_fil.append(total)
					#append_ids(total_partids)
					total = 0
					#total_partids = []
				else:
					Npart_per_fil.append(NumPartPerFil[i])
					#Npart_per_fil_partIDs.append(NumPartPerFil_particleIDs[i])
		else:
			if FilID[i] == FilID[i-1]:
				Npart_per_fil[-1] += NumPartPerFil[i]
				#if id_append_check:
				#	total_partids.append(NumPartPerFil_particleIDs[i])
				#	append_ids(total_partids)
				#else:
				#	raise ValueError("Something went wrong with id_append_check if test!")
				
			else:
				Npart_per_fil.append(NumPartPerFil[i])
				#Npart_per_fil_partIDs.append(NumPartPerFil_particleIDs[i])
	print 'Done'
	return np.array(Npart_per_fil), np.array(Filament_part_IDs)

if __name__ == '__main__':
	parsed_arguments = Argument_parser()
	HOMEPC = parsed_arguments.HOMEPC	# Set 1 if working in UiO terminal

	# Filament and dark matter particle plotting
	FilamentLimit = 0			# Limits the number of lines read from file. Reads all if 0
	PlotFilaments = 0			# Set 1 to plot actual filaments
	PlotFilamentsWCritPts = 0	# Set to 1 to plot filaments with critical points
	Projection2D = 0			# Set to 1 to plot a 2D projection of the 3D case
	FilamentColors = 1 			# Set to 1 to get different colors for different filaments
	ColorBarZDir = 1 			# Set 1 to include colorbar for z-direction
	ColorBarLength = 1 			# Set 1 to include colorbars based on length of the filament
	IncludeDMParticles = 1		# Set to 1 to include dark matter particle plots
	IncludeSlicing = 1			# Set 1 to include slices of the box
	MaskXdir = 0 				# Set 1 to mask one or more directions.
	MaskYdir = 0
	MaskZdir = 1

	# Histogram plots
	HistogramPlots = 0			# Set to 1 to plot histograms
	Comparison = parsed_arguments.Comparisons	# Set 1 if you want to compare different number of particles. Usual plots will not be plotted!
	ModelCompare = parsed_arguments.ModelCompare 	# Set to 1 to compare histograms of different models. Particle comparisons will not be run.
	SigmaComparison = parsed_arguments.SigmaComp 	# Set to 1 to compare histograms and/or plots based on different sigma values by MSE.
								# Must also set Comparison=1 to compare histograms
	
	# Run simulation for different models. Set to 1 to run them. 
	LCDM_model = 1 
	SymmA_model = 0
	SymmB_model = 0
	SymmC_model = 0
	SymmD_model = 0
	fofr4_model = 0
	fofr5_model = 0
	fofr6_model = 0
		
	# Global properties to be set
	IncludeUnits = 1			# Set to 1 to include 'rockstar' units, i.e Mpc/h and km/s
	SaveAsPNG = 1				# Set 1 to save figures as PNG
	SaveAsPDF = 0 				# Set 1 to save figures as PDF

	print '=== INFORMATION ==='
	# Some if tests before the simulation runs
	if FilamentLimit == 1:
		print '   Running program with limited amount of filaments.'

	if IncludeDMParticles == 1 and HOMEPC == 1:
		print '   Dark matter particles will be included'
		IncludeSlicing = 1
	else:
		print '   Dark matter particles not included'

	if IncludeSlicing == 1 and MaskXdir == 0 and MaskYdir == 0 and MaskZdir == 0:
		raise ValueError('IncludeSlicing set to 1, but all mask direction set to 0. Set at least one of them to 1.')

	if ModelCompare == 1:
		print '   Will compare histograms for different models'
	elif Comparison == 1 and ModelCompare == 0:
		print '   Will compare histograms over different models or number of particles.'
	if ModelCompare == 1 and HOMEPC:
		LCDM_model = 0

	UnitConverter = 256.0 if IncludeUnits == 1 else 1
	### Slice selection can be chosen here
	if IncludeSlicing == 1:
		print '   Slicing included'
		LowerBoundaryXDir = 0
		LowerBoundaryYDir = 0
		LowerBoundaryZDir = 0
		UpperBoundaryXDir = 0
		UpperBoundaryYDir = 0
		UpperBoundaryZDir = 0
		if MaskXdir == 1:
			LowerBoundaryXDir = 0.45*UnitConverter
			UpperBoundaryXDir = 0.55*UnitConverter
			print '   --Masking X direction'
		if MaskYdir == 1:
			LowerBoundaryYDir = 0.45*UnitConverter
			UpperBoundaryYDir = 0.55*UnitConverter
			print '   --Masking Y direction'
		if MaskZdir == 1:
			LowerBoundaryZDir = 0.45*UnitConverter
			UpperBoundaryZDir = 0.55*UnitConverter
			print '   --Masking Z direction'

	Mask_boundary_list = [UpperBoundaryXDir, UpperBoundaryYDir, UpperBoundaryZDir, LowerBoundaryXDir, LowerBoundaryYDir, LowerBoundaryZDir]
	Mask_direction_check = [MaskXdir, MaskYdir, MaskZdir]

	if SaveAsPNG == 1 and SaveAsPDF	== 1:
		raise ValueError('Cannot save both PDF and PNG at the same time. Only allow one at a time.')
	elif SaveAsPDF == 1 and SaveAsPNG == 0:
		print '   Saving figures as PDF files'
		filetype = '.pdf'
	elif SaveAsPNG == 1 and SaveAsPDF == 0:
		print '   Saving figures as PNG files'
		filetype = '.png'
	else:
		raise ValueError('Figure filetype to save not selected.')


	if HOMEPC == 0:
		file_directory = 'C:/Users/Alex/Documents/Masters_project/Disperse'
		savefile_directory = file_directory
		IncludeDMParticles = 0
		if LCDM_model == 1:
			print '=== Running for the LCDM model ==='
			LCDM_z0_64_dir = 'lcdm_z0_testing/LCDM_z0_64PeriodicTesting/'
			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'PlotsTest/', nPart=64, model='LCDM', redshift=0)
			NConn_64PartLCDM, FilLen_64PartLCDM, NPts_64PartLCDM = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl', Sigma_threshold=4.0)

			#LCDM_z0_128_dir = 'lcdm_z0_testing/LCDM_z0_128PeriodicTesting/'
			#LCDM_z0_128Instance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128, model='LCDM', redshift=0)
			#NConn_128PartLCDM, FilLen_128PartLCDM, NPts_128PartLCDM = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl')
			
			#LCDM_nsig4Instance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_64_dir+'Sigma4PlotsTest/', nPart=64, model='LCDM', redshift=0, SigmaArg=4)
			#NConn_nsig4, FilLen_nsig4, NPts_nsig4 = LCDM_nsig4Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064_nsig4.a.NDskl')
				
			if SigmaComparison:
				LCDM_nsig4Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Sigma4Plots/', nPart=64, model='LCDM', redshift=0, SigmaArg=4)
				LCDM_nsig5Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Sigma5Plots/', nPart=64, model='LCDM', redshift=0, SigmaArg=5)
				NConn_nsig4, FilLen_nsig4, NPts_nsig4 = LCDM_nsig4Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064_nsig4.a.NDskl')
				NConn_nsig5, FilLen_nsig5, NPts_nsig5 = LCDM_nsig5Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064_nsig5.a.NDskl')
				
			Comparison_dir = 'lcdm_z0_testing/Comparison_plots/'
			if Comparison == 1 and ModelCompare == 0:
				if SigmaComparison:
					NumConnections_list = [NConn_64PartLCDM, NConn_nsig4, NConn_nsig5]
					FilLengths_list = [FilLen_64PartLCDM, FilLen_nsig4, FilLen_nsig5]
					FilPoints_list = [NPts_64PartLCDM, NPts_nsig4, NPts_nsig5]
					Histogram_instance = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir+'SigmaComparisons/',\
									 redshift=0, LCDM=1, nsigComparison=1)
					Histogram_instance.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=64)
				else:
					NumConnections_list = [NConn_64PartLCDM, NConn_128PartLCDM]
					FilLengths_list = [FilLen_64PartLCDM, FilLen_128PartLCDM]
					FilPoints_list = [NPts_64PartLCDM, NPts_128PartLCDM]
					Histogram_instance = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, redshift=0, LCDM=1)
					Histogram_instance.Run(NumConnections_list, FilLengths_list, FilPoints_list)

		if SymmA_model:
			print '=== Running for Symm_A model ==='
			SymmA_z064_directory = 'SymmA_data/SymmA_z0_64Particles/'
			SymmA_z064_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmA_z064_directory+'Plots/', nPart=64, model='SymmA', redshift=0)
			Nconn_64PartSymmA, FilLen_64PartSymmA, NPts_64PartSymmA = SymmA_z064_instance.Solve(SymmA_z064_directory+'SkelconvOutput_SymmAz064Part.a.NDskl')

		if SymmB_model:
			print '=== Running for Symm_B model ==='
			SymmB_z064_directory = 'SymmB_data/SymmB_z0_64Particles/'
			SymmB_z064_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmB_z064_directory+'Plots/', nPart=64, model='SymmB', redshift=0)
			Nconn_64PartSymmB, FilLen_64PartSymmB, NPts_64PartSymmB = SymmB_z064_instance.Solve(SymmB_z064_directory+'SkelconvOutput_SymmBz064Part.a.NDskl')

		if ModelCompare:
			if LCDM_model and not SymmA_model and not SymmB_model:
				raise ValueError('Only LCDM model is ticked on. Cannot compare models if SymmA and/or SymmB is ticked as well')

			elif not LCDM_model and SymmA_model and not SymmB_model:
				raise ValueError('Only SymmA model is ticked on. Cannot compare models if LCDM and/or SymmB is ticked as well')

			elif not LCDM_model and not SymmA_model and SymmB_model:
				raise ValueError('Only SymmB model is ticked on. Cannot compare models if LCDM and/or SymmA is ticked as well')

			ModelComparisonDir = 'Model_comparisons/'
			ModelsFolder = ''
			if LCDM_model:
				ModelsFolder += 'LCDM'
			if SymmA_model:
				ModelsFolder += 'SymmA'
			if SymmB_model:
				ModelsFolder += 'SymmB'

			NumConnectionsModels_list = [NConn_64PartLCDM, Nconn_64PartSymmA, Nconn_64PartSymmB]
			FilLengthsModels_list = [FilLen_64PartLCDM, FilLen_64PartSymmA, FilLen_64PartSymmB]
			FilPointsModels_list = [NPts_64PartLCDM, NPts_64PartSymmA, NPts_64PartSymmB]

			ModelCompareInstance = HComp.Histogram_Comparison(savefile=1, savefile_directory = savefile_directory, 
							 savefigDirectory=ModelComparisonDir+ModelsFolder+'Plots/', filetype=filetype, redshift=0,\
							LCDM=LCDM_model, SymmA=SymmA_model, SymmB=SymmB_model)
			ModelCompareInstance.Run(NumConnectionsModels_list, FilLengthsModels_list, FilPointsModels_list, nPart=64)

	if HOMEPC == 1:
		file_directory = '/mn/stornext/d5/aleh'
		savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
		
		if LCDM_model == 1:
			print '=== Running for the LCDM model ==='
			solve_file_dir = '/lcdm_testing/'
			solve_filename = 'lcdm_z0_test.solve'
			
			if IncludeDMParticles == 1:
				Gadget_instance = ReadGadgetFile.Read_Gadget_file(Mask_direction_check, Mask_boundary_list)
				PartPosX, PartPosY, PartPosZ, PartIDs, DMHistogram, DMBinXedges, DMBinYedges, DM_KDTree = Gadget_instance.Get_particles('lcdm', includeKDTree=False)
				"""
				#SolveReadInstance = Read_solve_files()
				SolveReadInstance = Read_Gadget_file()
				#ParticlePos = SolveReadInstance.ParticlePos
				PartPosX = SolveReadInstance.PartPosX
				PartPosY = SolveReadInstance.PartPosY
				PartPosZ = SolveReadInstance.PartPosZ
				DM_KDTree = SolveReadInstance.DM_tree
				"""
			
			LCDM_z0_64_dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/'
			LCDM_z0_128_dir = 'lcdm_testing/LCDM_z0_128PeriodicTesting/'
			LCDM_z0_256_dir = 'lcdm_testing/LCDM_z0_256PeriodicTesting/'
			LCDM_z0_512_dir = 'lcdm_testing/LCDM_z0_512PeriodicTesting/'

			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plotstest2_periodic_png/', nPart=128, model='LCDM', redshift=0)
			NumConn_64LCDM, FilLen_64LCDM, NPts_64LCDM = LCDM_z0_64Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDMz0128.a.NDskl')
			
			if parsed_arguments.HigherPart:
				LCDM_z0_128Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plotstest2_png/', nPart=128, model='LCDM', redshift=0)
				NumConn_128LCDM, FilLen_128LCDM, NPts_128LCDM = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl')
				
				LCDM_z0_256Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/', nPart=256, model='LCDM', redshift=0)
				NumConn_256LCDM, FilLen_256LCDM, NPts_256LCDM = LCDM_z0_256Instance.Solve(LCDM_z0_256_dir+'SkelconvOutput_LCDMz0256.a.NDskl')
				
				LCDM_z0_512Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/', nPart=512, model='LCDM', redshift=0)
				NumConn_512LCDM, FilLen_512LCDM, NPts_512LCDM = LCDM_z0_512Instance.Solve(LCDM_z0_512_dir+'SkelconvOutput_LCDMz0512.a.NDskl')
				

			#LCDM_z0_64_Robustness_dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/Robustness_argument/'
			#LCDM_z0_64_RobustnessInstance = Disperse_Plotter(savefile=1, \
			#	savefigDirectory=LCDM_z0_64_Robustness_dir+'Plots_Robustness/', nPart=64, model='LCDM', redshift=0)
			#NumConn_64LCDM_rob, FilLen_64LCDM_rob, NPts_64LCDM_rob = \
			#  LCDM_z0_64_RobustnessInstance.Solve(LCDM_z0_64_Robustness_dir+'SkelconvOutput_LCDMz064_robustness3.TRIM.a.NDskl', robustness=True)

			if SigmaComparison:
				LCDM64_instance_nsig4 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plotstest2_png/Sigma4/', nPart=64, model='LCDM', redshift=0, SigmaArg=4)
				LCDM64_instance_nsig5 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plotstest2_png/Sigma5/', nPart=64, model='LCDM', redshift=0, SigmaArg=5)
				NConn_64nsig4, FilLen_64nsig4, NPts_64nsig4 = LCDM64_instance_nsig4.Solve(LCDM_z0_64_dir+'Sigma4/SkelconvOutput_LCDMz064_nsig4.a.NDskl')
				NConn_64nsig5, FilLen_64nsig5, NPts_64nsig5 = LCDM64_instance_nsig5.Solve(LCDM_z0_64_dir+'Sigma5/SkelconvOutput_LCDMz064_nsig5.a.NDskl')
				if parsed_arguments.HigherPart:
					LCDM128_instance_nsig4 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plotstest2_png/Sigma4/', nPart=128, model='LCDM', redshift=0, SigmaArg=4)
					LCDM128_instance_nsig5 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plotstest2_png/Sigma5/', nPart=128, model='LCDM', redshift=0, SigmaArg=5)
					NConn_128nsig4, FilLen_128nsig4, NPts_128nsig4 = LCDM128_instance_nsig4.Solve(LCDM_z0_128_dir+'Sigma4/SkelconvOutput_LCDMz0128_nsig4.a.NDskl')
					NConn_128nsig5, FilLen_128nsig5, NPts_128nsig5 = LCDM128_instance_nsig5.Solve(LCDM_z0_128_dir+'Sigma5/SkelconvOutput_LCDMz0128_nsig5.a.NDskl')

					LCDM256_instance_nsig4 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/Sigma4/', nPart=256, model='LCDM', redshift=0, SigmaArg=4)
					LCDM256_instance_nsig5 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/Sigma5/', nPart=256, model='LCDM', redshift=0, SigmaArg=5)
					LCDM256_instance_nsig6 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/Sigma6/', nPart=256, model='LCDM', redshift=0, SigmaArg=6)
					LCDM256_instance_nsig7 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/Sigma7/', nPart=256, model='LCDM', redshift=0, SigmaArg=7)
					LCDM256_instance_nsig8 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plotstest2_png/Sigma8/', nPart=256, model='LCDM', redshift=0, SigmaArg=8)
					NConn_256nsig4, FilLen_256nsig4, NPts_256nsig4 = LCDM256_instance_nsig4.Solve(LCDM_z0_256_dir+'Sigma4/SkelconvOutput_LCDMz0256_nsig4.a.NDskl')
					NConn_256nsig5, FilLen_256nsig5, NPts_256nsig5 = LCDM256_instance_nsig5.Solve(LCDM_z0_256_dir+'Sigma5/SkelconvOutput_LCDMz0256_nsig5.a.NDskl')
					NConn_256nsig6, FilLen_256nsig6, NPts_256nsig6 = LCDM256_instance_nsig6.Solve(LCDM_z0_256_dir+'Sigma6/SkelconvOutput_LCDMz0256_nsig6.a.NDskl')
					NConn_256nsig7, FilLen_256nsig7, NPts_256nsig7 = LCDM256_instance_nsig7.Solve(LCDM_z0_256_dir+'Sigma7/SkelconvOutput_LCDMz0256_nsig7.a.NDskl')
					NConn_256nsig8, FilLen_256nsig8, NPts_256nsig8 = LCDM256_instance_nsig8.Solve(LCDM_z0_256_dir+'Sigma8/SkelconvOutput_LCDMz0256_nsig8.a.NDskl')

					LCDM512_instance_nsig4 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/Sigma4/', nPart=512, model='LCDM', redshift=0, SigmaArg=4)
					LCDM512_instance_nsig5 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/Sigma5/', nPart=512, model='LCDM', redshift=0, SigmaArg=5)
					LCDM512_instance_nsig6 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/Sigma6/', nPart=512, model='LCDM', redshift=0, SigmaArg=6)
					LCDM512_instance_nsig7 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/Sigma7/', nPart=512, model='LCDM', redshift=0, SigmaArg=7)
					LCDM512_instance_nsig8 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plotstest2_png/Sigma8/', nPart=512, model='LCDM', redshift=0, SigmaArg=8)
					NConn_512nsig4, FilLen_512nsig4, NPts_512nsig4 = LCDM512_instance_nsig4.Solve(LCDM_z0_512_dir+'Sigma4/SkelconvOutput_LCDMz0512_nsig4.a.NDskl')
					NConn_512nsig5, FilLen_512nsig5, NPts_512nsig5 = LCDM512_instance_nsig5.Solve(LCDM_z0_512_dir+'Sigma5/SkelconvOutput_LCDMz0512_nsig5.a.NDskl')
					NConn_512nsig6, FilLen_512nsig6, NPts_512nsig6 = LCDM512_instance_nsig6.Solve(LCDM_z0_512_dir+'Sigma6/SkelconvOutput_LCDMz0512_nsig6.a.NDskl')
					NConn_512nsig7, FilLen_512nsig7, NPts_512nsig7 = LCDM512_instance_nsig7.Solve(LCDM_z0_512_dir+'Sigma7/SkelconvOutput_LCDMz0512_nsig7.a.NDskl')
					NConn_512nsig8, FilLen_512nsig8, NPts_512nsig8 = LCDM512_instance_nsig8.Solve(LCDM_z0_512_dir+'Sigma8/SkelconvOutput_LCDMz0512_nsig8.a.NDskl')

			Comparison_dir = 'lcdm_testing/Comparison_plots/'
			if Comparison == 1 and ModelCompare == 0:
				if SigmaComparison:
					NumConnections_list = [NumConn_512LCDM, NConn_512nsig4, NConn_512nsig5]
					FilLengths_list = [FilLen_512LCDM, FilLen_512nsig4, FilLen_512nsig5]
					FilPoints_list = [NPts_512LCDM, NPts_512nsig4, NPts_512nsig5]

					NumConnections_list_expanded = [NConn_64nsig5, NConn_64nsig4, NumConn_64LCDM ,NConn_512nsig5, NConn_512nsig4, NumConn_512LCDM]
					FilLengths_list_expanded = [FilLen_64LCDM, FilLen_64nsig4, FilLen_64nsig5, FilLen_512LCDM, FilLen_512nsig4, FilLen_512nsig5]
					FilPoints_list_expanded = [NPts_64LCDM, NPts_64nsig4, NPts_64nsig5, NPts_512LCDM, NPts_512nsig4, NPts_512nsig5]

					Nconnections_list_nsig3 = [NumConn_64LCDM, NumConn_128LCDM, NumConn_256LCDM, NumConn_512LCDM]
					Nconnections_list_nsig4 = [NConn_64nsig4, NConn_128nsig4, NConn_256nsig4, NConn_512nsig4]
					Nconnections_list_nsig5 = [NConn_64nsig5, NConn_128nsig5, NConn_256nsig5, NConn_512nsig5]
					FilLengths_list_nsig3 = [FilLen_64LCDM, FilLen_128LCDM, FilLen_256LCDM, FilLen_512LCDM]
					FilLengths_list_nsig4 = [FilLen_64nsig4, FilLen_128nsig4, FilLen_256nsig4, FilLen_512nsig4]
					FilLengths_list_nsig5 = [FilLen_64nsig5, FilLen_128nsig5, FilLen_256nsig5, FilLen_512nsig5]
					FilPoints_list_nsig3 = [NPts_64LCDM, NPts_128LCDM, NPts_256LCDM, NPts_512LCDM]
					FilPoints_list_nsig4 = [NPts_64nsig4, NPts_128nsig4, NPts_256nsig4, NPts_512nsig4]
					FilPoints_list_nsig5 = [NPts_64nsig5, NPts_128nsig5, NPts_256nsig5, NPts_512nsig5]					

					ComparisonInstance_LCDM = HComp.Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir+'SigmaComparisons',\
							savefile_directory = savefile_directory, filetype=filetype, redshift=0, LCDM=1, nsigComparison=1)
					ComparisonInstance_LCDM.Convergence_tests(NumConnections_list_expanded, FilLengths_list_expanded, FilPoints_list_expanded)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig3, FilLengths_list_nsig3, FilPoints_list_nsig3, 3)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig4, FilLengths_list_nsig4, FilPoints_list_nsig5, 4)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig5, FilLengths_list_nsig5, FilPoints_list_nsig4, 5)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=512)
				elif parsed_arguments.HigherPart:
					NumConnections_list = [NumConn_64LCDM, NumConn_128LCDM , NumConn_256LCDM, NumConn_512LCDM]
					FilLengths_list = [FilLen_64LCDM, FilLen_128LCDM, FilLen_256LCDM, FilLen_512LCDM]
					FilPoints_list = [NPts_64LCDM, NPts_128LCDM, NPts_256LCDM, NPts_512LCDM]
					ComparisonInstance_LCDM = Hcomp.Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, \
							savefile_directory=savefile_directory, filetype=filetype, redshift=0)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list)
				else:
					SkeletonFiles = [LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl', LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl',\
									LCDM_z0_256_dir+'SkelconvOutput_LCDMz0256.a.NDskl', LCDM_z0_512_dir+'SkelconvOutput_LCDMz0512.a.NDskl']
					
					results_dir_filpersig = os.path.join(savefile_directory, Comparison_dir+'Plotstest/')
					if not os.path.isdir(results_dir_filpersig):
						os.makedirs(results_dir_filpersig)
					
					N_sigmas = 100
					sigma_values = np.linspace(3, 10, N_sigmas)
					Results = []
					for filenames in SkeletonFiles:
						Instance = FilamentsPerSigma.FilamentsPerSigma(filenames)
						Instance_output = Instance.Filaments_per_sigma_Maxima(sigma_values)
						Results.append(Instance_output)

					fig_filpersig, ax_filpersig = plt.subplots()
					fig_filpersig_log, ax_filpersig_log = plt.subplots()
					for data in Results:
						ax_filpersig.plot(sigma_values, data)
						ax_filpersig_log.semilogy(sigma_values, data)
					ax_filpersig.legend(['$\mathregular{64^3}$ subsample','$\mathregular{128^3}$ subsample',\
								'$\mathregular{256^3}$ subsample','$\mathregular{512^3}$ particles'])
					ax_filpersig_log.legend(['$\mathregular{64^3}$ subsample','$\mathregular{128^3}$ subsample',\
								'$\mathregular{256^3}$ subsample','$\mathregular{512^3}$ particles'])
					ax_filpersig.set_xlabel('Sigma')
					ax_filpersig_log.set_xlabel('Sigma')
					ax_filpersig.set_ylabel('Number of filaments')
					ax_filpersig_log.set_ylabel('Number of filaments')

					print 'saving in', results_dir_filpersig
					fig_filpersig.savefig(results_dir_filpersig + 'Filaments_per_sigma_Maxima' + str(N_sigmas) + '.png')
					fig_filpersig_log.savefig(results_dir_filpersig + 'Filaments_per_sigma_Maxima'+ str(N_sigmas) + '_logarithmic.png')

	if ModelCompare and HOMEPC == 1:
		""" 
		Runs for all modified gravity simulations.
		All runs are on 256^3 particle subsample at sigma=5.
		"""
		lcdm_dir = 'lcdm_testing/LCDM_z0_128Particles/Sigma3/'
		SymmA_dir = 'SymmA_data/SymmA_z0_128Particles/Sigma3/'
		SymmB_dir = 'SymmB_data/SymmB_z0_128Particles/Sigma3/'
		SymmC_dir = 'SymmC_data/SymmC_z0_128Particles/Sigma3/'
		SymmD_dir = 'SymmD_data/SymmD_z0_128Particles/Sigma3/'
		fofr4_dir = 'fofr4_data/fofr4_z0_128Particles/Sigma3/'
		fofr5_dir = 'fofr5_data/fofr5_z0_128Particles/Sigma3/'
		fofr6_dir = 'fofr6_data/fofr6_z0_128Particles/Sigma3/'
		
		LCDM_instance = Disperse_Plotter(savefile=2, savefigDirectory=lcdm_dir+'Plots/', nPart=128, model='LCDM', redshift=0, SigmaArg=3)
		NumConn_LCDM, FilLen_LCDM, NPts_LCDM = LCDM_instance.Solve(lcdm_dir+'SkelconvOutput_LCDMz0128_nsig3.a.NDskl')
		Fil3DPos_LCDM, FilID_LCDM = LCDM_instance.get_3D_pos()
		"""
		SymmA_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmA_dir+'Plots/', nPart=128, model='SymmA', redshift=0, SigmaArg=3)
		NummConn_SymmA, FilLen_SymmA, NPts_SymmA = SymmA_instance.Solve(SymmA_dir+'SkelconvOutput_SymmAz0128_nsig3.a.NDskl')
		Fil3DPos_SymmA, FilID_SymmA = SymmA_instance.get_3D_pos()

		SymmB_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmB_dir+'Plots/', nPart=128, model='SymmB', redshift=0, SigmaArg=3)
		NummConn_SymmB, FilLen_SymmB, NPts_SymmB = SymmB_instance.Solve(SymmB_dir+'SkelconvOutput_SymmBz0128_nsig3.a.NDskl')
		Fil3DPos_SymmB, FilID_SymmB = SymmB_instance.get_3D_pos()
		
		SymmC_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmC_dir+'Plots/', nPart=128, model='SymmC', redshift=0, SigmaArg=3)
		NummConn_SymmC, FilLen_SymmC, NPts_SymmC = SymmC_instance.Solve(SymmC_dir+'SkelconvOutput_SymmCz0128_nsig3.a.NDskl')
		Fil3DPos_SymmC, FilID_SymmC = SymmC_instance.get_3D_pos()
		
		SymmD_instance = Disperse_Plotter(savefile=2, savefigDirectory=SymmD_dir+'Plots/', nPart=128, model='SymmD', redshift=0, SigmaArg=3)
		NummConn_SymmD, FilLen_SymmD, NPts_SymmD = SymmD_instance.Solve(SymmD_dir+'SkelconvOutput_SymmDz0128_nsig3.a.NDskl')
		Fil3DPos_SymmD, FilID_SymmD = SymmD_instance.get_3D_pos()
		
		fofr4_instance = Disperse_Plotter(savefile=2, savefigDirectory=fofr4_dir+'Plots/', nPart=128, model='fofr4', redshift=0, SigmaArg=3)
		NummConn_fofr4, FilLen_fofr4, NPts_fofr4 = fofr4_instance.Solve(fofr4_dir+'SkelconvOutput_fofr4z0128_nsig3.a.NDskl')
		Fil3dPos_fofr4, FilID_fofr4 = fofr4_instance.get_3D_pos()

		fofr5_instance = Disperse_Plotter(savefile=2, savefigDirectory=fofr5_dir+'Plots/', nPart=128, model='fofr5', redshift=0, SigmaArg=3)
		NummConn_fofr5, FilLen_fofr5, NPts_fofr5 = fofr5_instance.Solve(fofr5_dir+'SkelconvOutput_fofr5z0128_nsig3.a.NDskl')
		Fil3dPos_fofr5, FilID_fofr5 = fofr5_instance.get_3D_pos()
		
		fofr6_instance = Disperse_Plotter(savefile=2, savefigDirectory=fofr6_dir+'Plots/', nPart=128, model='fofr6', redshift=0, SigmaArg=3)
		NummConn_fofr6, FilLen_fofr6, NPts_fofr6 = fofr6_instance.Solve(fofr6_dir+'SkelconvOutput_fofr6z0128_nsig3.a.NDskl')
		Fil3dPos_fofr6, FilID_fofr6 = fofr6_instance.get_3D_pos()
		"""
		"""
		lcdm_64dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/'
		symma_64dir = 'SymmA_data/SymmA_z0_64Particles/'
		LCDM_test = Disperse_Plotter(savefile=2, savefigDirectory=lcdm_64dir+'Plots/', nPart=64, model='LCDM', redshift=0)
		NNLCDM, FLLCDM, NPLCDM = LCDM_test.Solve(lcdm_64dir+'SkelconvOutput_LCDMz064.a.NDskl')

		syma_test = Disperse_Plotter(savefile=2, savefigDirectory=symma_64dir+'Plots/', nPart=64, model='SymmA', redshift=0)
		NNSymmA, FLSymmA, NPSymmA = syma_test.Solve(symma_64dir+'SkelconvOutput_SymmAz064Part.a.NDskl')

		NumConnections_list = [NNLCDM, NNSymmA]
		FilLengths_list = [FLLCDM, FLSymmA]
		FilPoints_list = [NPLCDM, NPSymmA]

		CompI = HComp.Histogram_Comparison(savefile=0, savefigDirectory='LOL/', savefile_directory='WHAT/', filetype=filetype, redshift=0, LCDM=1)
		CompI.Compare_mg_models(NumConnections_list, FilLengths_list, FilPoints_list)
		"""
		# Compute number of particles per filament.
		# Units in Mpc/h for distance threshold and box_expand
		time_proc_start = time.clock()
		proc = mp.Pool(parsed_arguments.NumProcesses)
		print 'Processor preparation: ', time.clock() - time_proc_start, 's'
		BoundaryCheck_list2 = [0,0,0,0,0,0]
		Mask_check_list2 = [0,0,0]
		distance_threshold = 0.3
		box_expand = 1
		if parsed_arguments.NumPartModel == 'lcdm':
			NumPartPerFil_LCDM, NPPF_ids_LCDM = Save_NumPartPerFil('lcdm', Fil3DPos_LCDM[0:4], FilID_LCDM[0:4], 128, 3)
		elif parsed_arguments.NumPartModel == 'symmA':
			NumPartPerFil_SymmA, NPPF_ids_SymmA = Save_NumPartPerFil('symm_A', Fil3DPos_SymmA, FilID_SymmA, 128, 3)
		elif parsed_arguments.NumPartModel == 'symmB':
			NumPartPerFil_SymmB, NPPF_ids_SymmB = Save_NumPartPerFil('symm_B', Fil3DPos_SymmB, FilID_SymmB, 128, 3)
		elif parsed_arguments.NumPartModel == 'symmC':
			NumPartPerFil_SymmC, NPPF_ids_SymmC = Save_NumPartPerFil('symm_C', Fil3DPos_SymmC, FilID_SymmC, 128, 3)	
		elif parsed_arguments.NumPartModel == 'symmD':
			NumPartPerFil_SymmD, NPPF_ids_SymmD = Save_NumPartPerFil('symm_D', Fil3DPos_SymmD, FilID_SymmD, 128, 3)
		elif parsed_arguments.NumPartModel == 'fofr4':
			NumPartPerFil_fofr4, NPPF_ids_fofr4 = Save_NumPartPerFil('fofr4', Fil3DPos_fofr4, FilID_fofr4, 128, 3)
		elif parsed_arguments.NumPartModel == 'fofr5':
			NumPartPerFil_fofr5, NPPF_ids_fofr5 = Save_NumPartPerFil('fofr5', Fil3DPos_fofr5, FilID_fofr5, 128, 3)
		elif parsed_arguments.NumPartModel == 'fofr6':
			NumPartPerFil_fofr6, NPPF_ids_fofr6 = Save_NumPartPerFil('fofr6', Fil3DPos_fofr6, FilID_fofr6, 128, 3)
		else:
			NumPartPerFil_LCDM, NPPF_ids_LCDM = Save_NumPartPerFil('lcdm', Fil3DPos_LCDM, FilID_LCDM, 128, 3)
			NumPartPerFil_SymmA, NPPF_ids_SymmA = Save_NumPartPerFil('symm_A', Fil3DPos_SymmA, FilID_SymmA, 128, 3)
			NumPartPerFil_SymmB, NPPF_ids_SymmB = Save_NumPartPerFil('symm_B', Fil3DPos_SymmB, FilID_SymmB, 128, 3)
			NumPartPerFil_SymmC, NPPF_ids_SymmC = Save_NumPartPerFil('symm_C', Fil3DPos_SymmC, FilID_SymmC, 128, 3)
			NumPartPerFil_SymmD, NPPF_ids_SymmD = Save_NumPartPerFil('symm_D', Fil3DPos_SymmD, FilID_SymmD, 128, 3)
			NumPartPerFil_fofr4, NPPF_ids_fofr4 = Save_NumPartPerFil('fofr4', Fil3DPos_fofr4, FilID_fofr4, 128, 3)
			NumPartPerFil_fofr5, NPPF_ids_fofr5 = Save_NumPartPerFil('fofr5', Fil3DPos_fofr5, FilID_fofr5, 128, 3)
			NumPartPerFil_fofr6, NPPF_ids_fofr6 = Save_NumPartPerFil('fofr6', Fil3DPos_fofr6, FilID_fofr6, 128, 3)
		proc.close()
		proc.join()