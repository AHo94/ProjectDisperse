### Set comment on the two below to plot. Only use if running on euclid21 or other institute computers. 
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
#import scipy.stats as stats
from scipy import spatial
import os
import time
import cPickle as pickle
from collections import Counter

# Disable warnings from matplotlib	
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# Own modules
import BoundaryChecker
import FilamentMasking
import ReadGadgetFile
import MaskCritPts

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
			self.nPart_text = '$\mathregular{64}$'
		elif self.nPart == 128:
			self.nPart_text = '$\mathregular{128}$'
		elif self.nPart == 256:
			self.nPart_text = '$\mathregular{256}$'
		elif self.nPart == 512:
			self.nPart_text = '$\mathregular{512}$'
		else:
			raise ValueError('Argument nPart must be 64, 128, 256 or 512!')

		if type(savefigDirectory) != str:
			raise ValueError('Argument savefigDirectory must be a string!')

		# Creates folder when saving figure
		#script_dir = os.path.dirname(__file__) 	# Script running directory
		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir):
			os.makedirs(self.results_dir)

		if type(redshift) != int and type(redshift) != float:
			raise ValueError('Value of redshift must be an integer or float')

		if model == 'LCDM':
			self.ModelFilename = 'LCDM' + 'z' + str(redshift) + filetype
		elif model == 'SymmA':
			self.ModelFilename = 'SymmA' + 'z' + str(redshift) + filetype
		elif model == 'SymmB':
			self.ModelFilename = 'SymmB' + 'z' + str(redshift) + filetype
		else:
			raise ValueError('Model name not properly set')

		self.SigmaTitle = ''
		if SigmaComparison:
			if not SigmaArg:
				self.SigmaTitle = '\mathregular{\sigma} = 3'
			else:
				self.SigmaTitle = '\mathregular{\sigma} = ' + str(SigmaArg)


		if SigmaArg:
			self.SigmaArg = SigmaArg
		else:
			self.SigmaArg = False

	def ReadFile(self, filename, dimensions):
		""" Reads the data from the skeleton file from Disperse """
		time_start = time.clock()
		print 'Reading data for the file: ', filename, '...' 
		datafiles = open(os.path.join(file_directory, filename), 'r')
		self.CriticalPoints = []
		self.Filaments = []
		self.CriticalPointsData = []
		self.FilamentsData = []

		for line in datafiles:
			data_set = line.split()
			if 'BBOX' in line:
				BoxSize = line.split()
				BoxMin = BoxSize[1]
				BoxMax = BoxSize[2]
				MaxValues = BoxMax[BoxMax.index("[") + 1:BoxMax.rindex("]")].replace(",", " ").split()
				MinValues = BoxMin[BoxMin.index("[") + 1:BoxMin.rindex("]")].replace(",", " ").split()
				self.xmin = float(MinValues[0])*UnitConverter
				self.xmax = float(MaxValues[0])*UnitConverter
				self.ymin = float(MinValues[1])*UnitConverter
				self.ymax = float(MaxValues[1])*UnitConverter
				if dimensions == 3:
					self.zmin = float(MinValues[2])*UnitConverter
					self.zmax = float(MaxValues[2])*UnitConverter
				else:
					continue
				
			if '[CRITICAL POINTS]' in line:
				SETLIMIT = 0
				for lineCrit in datafiles:
					dataCritpts = lineCrit.split()
					if '[FILAMENTS]' in lineCrit:
						line = lineCrit
						break
					else:
						self.CriticalPoints.append(dataCritpts)
			if '[FILAMENTS]' in line:
				SETLIMIT = 0
				for lineFil in datafiles:
					dataFil = lineFil.split()
					if '[CRITICAL POINTS DATA]' in lineFil:
						line = lineFil
						break
					else:
						self.Filaments.append(dataFil)
					if FilamentLimit != 0:
						if SETLIMIT == FilamentLimit+1:
							break
						SETLIMIT += 1

			if '[CRITICAL POINTS DATA]' in line:
				for lineCritPtsData in datafiles:
					dataCritPtsData = lineCritPtsData.split()
					if '[FILAMENTS DATA]' in lineCritPtsData:
						line = lineCritPtsData
						break
					else:
						self.CriticalPointsData.append(dataCritPtsData)
			if '[FILAMENTS DATA]' in line:
				for lineFilData in datafiles:
					dataFilData = lineFilData.split()
					self.FilamentsData.append(dataFilData)

		datafiles.close()
		self.CriticalPoints = np.asarray(self.CriticalPoints)
		self.Filaments = np.asarray(self.Filaments)
		self.CriticalPointsData = np.asarray(self.CriticalPointsData)
		self.FilamentsData = np.asarray(self.FilamentsData)
		print 'Time elapsed for filament data reading: ', time.clock() - time_start, 's'

	def Sort_filament_coordinates(self, dimensions):
		""" 
		Sorts the coordinates of the filaments and critical points to their respective arrays 
		Data to be sorted: Critical points coordinate, ID of filament and filament coordinates
		"""
		time_start = time.clock()
		print 'Sorting filament and critical point positions ...'
		if FilamentLimit == 0:
			self.NFils = int(self.Filaments[0][0])
		else:
			self.NFils = FilamentLimit

		# General data
		# Critical points
		self.CritPointXpos = []
		self.CritPointYpos = []
		self.CritPointZpos = []
		self.NcritPts = int(self.CriticalPoints[0][0])
		self.Critpts_filamentID = []
		self.Neighbours_CP = []
		self.CP_id_of_connecting_filament = []
		self.Number_filaments_connecting_to_CP = []
		i = 1
		counter = 0
		while i < len(self.CriticalPoints):
			Info = self.CriticalPoints[i]	# Type, value and neighbouring CP not saved
			self.CritPointXpos.append(float(Info[1])*UnitConverter)
			self.CritPointYpos.append(float(Info[2])*UnitConverter)
			self.CritPointZpos.append(float(Info[3])*UnitConverter)
			self.Neighbours_CP.append(int(Info[5]))
			Filament_connections = int(self.CriticalPoints[i+1][0])
			self.Number_filaments_connecting_to_CP.append(Filament_connections)
			Temp_filID = []
			Temp_CPID = []
			if Filament_connections != 0:
				for j in range(1,Filament_connections+1):
					Temp_CPID.append(int(self.CriticalPoints[i+j+1][0]))
					Temp_filID.append(int(self.CriticalPoints[i+j+1][1]))
				self.Critpts_filamentID.append(np.array(Temp_filID))
				self.CP_id_of_connecting_filament.append(np.array(Temp_CPID))
			#else:
			#	break
			counter += 1
			i += 2 + Filament_connections

		"""
		self.NumFilamentConnections = []
		self.IDFilamentConnected = []
		self.CritPointInfo = []
		############ UNIT CONVERSION NOT FIXED HERE!!!! 
		self.NcritPts = int(self.CriticalPoints[0][0])
		for i in range(1, len(self.CriticalPoints)-1):
			stuff = self.CriticalPoints[i]
			if len(stuff) == 1:
				self.NumFilamentConnections.append(int(stuff[0]))
				IDconnections = []
				if int(stuff[0]) == 1:
					IDS = self.CriticalPoints[i+1]
					IDconnections.append(np.array([float(IDS[0]), float(IDS[1])]))
				else:	
					for j in range(1, int(stuff[0])+1):
						IDS = self.CriticalPoints[i+j]
						IDconnections.append(np.array([float(IDS[0]),float(IDS[1])]))
				self.IDFilamentConnected.append(np.asarray(IDconnections))
			elif len(stuff) > 2:
				InfoListTemp = []
				for k in range(len(stuff)):
					InfoListTemp.append(float(stuff[k]))
				self.CritPointInfo.append(np.asarray(InfoListTemp))

		self.NumFilamentConnections = np.asarray(self.NumFilamentConnections)
		self.IDFilamentConnected = np.asarray(self.IDFilamentConnected)
		self.CritPointInfo = np.asarray(self.CritPointInfo)
		"""

		# Filament positions etc
		self.FilamentPos = []
		self.FilID = []
		self.PairIDS = []
		self.xdimPos = []
		self.ydimPos = []
		self.NFilamentPoints = []
		k = 1
		NewID = 1
		if dimensions == 2:
			for i in range(1, self.NFils):
				Filstuff = self.Filaments[k]	# Contains info on filament ID, num crit pts etc.
				TempPositions = []
				xtemp = []
				ytemp = []
				self.FilID.append(NewID)
				self.PairIDS.append(np.array([Filstuff[0], Filstuff[1]]))
				for j in range(1, int(Filstuff[-1])+1):
					xPos = float(self.Filaments[k+j][0])*UnitConverter
					yPos = float(self.Filaments[k+j][1])*UnitConverter
					TempPositions.append(np.array([xPos, yPos]))
					xtemp.append(xPos)
					ytemp.append(yPos)
				self.FilamentPos.append(np.asarray(TempPositions))
				self.xdimPos.append(np.asarray(xtemp))
				self.ydimPos.append(np.asarray(ytemp))
				k += int(Filstuff[-1])+1
				NewID += 1
				if FilamentLimit != 0:
					if k > FilamentLimit:
						break
		elif dimensions == 3:
			self.zdimPos = []
			for i in range(1, self.NFils):
				Filstuff = self.Filaments[k]
				TempPositions = []
				xtemp = []
				ytemp = []
				ztemp = []
				self.FilID.append(NewID)
				self.PairIDS.append(np.array([int(Filstuff[0]), int(Filstuff[1])]))
				for j in range(1, int(Filstuff[-1])+1):
					if k+j >= len(self.Filaments):
						break
					xPos = float(self.Filaments[k+j][0])*UnitConverter
					yPos = float(self.Filaments[k+j][1])*UnitConverter
					zPos = float(self.Filaments[k+j][2])*UnitConverter
					TempPositions.append(np.array([xPos, yPos]))
					xtemp.append(xPos)
					ytemp.append(yPos)
					ztemp.append(zPos)		
				self.NFilamentPoints.append(int(Filstuff[-1]))
				self.FilamentPos.append(np.array(TempPositions))
				self.xdimPos.append(np.array(xtemp))
				self.ydimPos.append(np.array(ytemp))
				self.zdimPos.append(np.array(ztemp))
				k += int(Filstuff[-1])+1
				NewID += 1
				if FilamentLimit != 0:
					if k >= FilamentLimit:
						self.NFils = len(self.xdimPos)
						break


		self.FilamentPos = np.array(self.FilamentPos)
		self.NFilamentPoints = np.array(self.NFilamentPoints)
		self.xdimPos = np.array(self.xdimPos)
		self.ydimPos = np.array(self.ydimPos)
		self.Critpts_filamentID = np.array(self.Critpts_filamentID)
		self.CritPointXpos = np.asarray(self.CritPointXpos)
		self.CritPointYpos = np.asarray(self.CritPointYpos)
		self.CritPointZpos = np.asarray(self.CritPointZpos)
		print self.NFils, 'NUMBER OF FILAMENTS'
		if dimensions == 3:
			self.zdimPos = np.array(self.zdimPos)

		print 'Array sorting time: ', time.clock() - time_start, 's'

	def Sort_filament_data(self):
		""" 
		Sorts the data of the filaments and critical points.
		Data to be sorted includes: Persistence, persistence pairs, persistence ratios, field values, etc.
		"""
		print 'Sorting data of filaments and critical points'
		header_count_critpoint = int(self.CriticalPointsData[0][0])
		self.Persistence_ratio = []
		self.Persistence_nsigmas = []
		self.Persistence = []
		self.Persistence_pair = []
		self.Parent_index = []
		self.Parent_log_index = []
		self.log_CritPoint_field_value = []
		self.CritPoint_field_value = []
		self.CritPoint_cell = []

		for i in range(1, len(self.CriticalPointsData)):
			if i <= header_count_critpoint:
				continue
			else:
				if i <= len(self.Critpts_filamentID)+header_count_critpoint:
					# Do not include critical points that are not connected to filaments
					self.Persistence_ratio.append(float(self.CriticalPointsData[i][0]))
					self.Persistence_nsigmas.append(float(self.CriticalPointsData[i][1]))
					self.Persistence.append(float(self.CriticalPointsData[i][2]))
					self.Persistence_pair.append(int(self.CriticalPointsData[i][3]))
					self.Parent_index.append(int(self.CriticalPointsData[i][4]))
					self.Parent_log_index.append(int(self.CriticalPointsData[i][5]))
					self.log_CritPoint_field_value.append(float(self.CriticalPointsData[i][6]))
					self.CritPoint_field_value.append(float(self.CriticalPointsData[i][7]))
					self.CritPoint_cell.append(float(self.CriticalPointsData[i][8]))

		header_count_filament = int(self.FilamentsData[0][0])
		self.Filament_field_value = []
		self.orientation = []
		self.Filament_cell = []
		self.log_Filament_field_value = []
		self.Filament_type = []

		for j in range(1, len(self.FilamentsData)):
			if j <= header_count_filament:
				continue
			else:
				self.Filament_field_value.append(float(self.FilamentsData[i][0]))
				self.orientation.append(int(self.FilamentsData[i][1]))
				self.Filament_cell.append(float(self.FilamentsData[i][2]))
				self.log_Filament_field_value.append(float(self.FilamentsData[i][3]))
				self.Filament_type.append(int(self.FilamentsData[i][4]))

		self.Persistence_nsigmas = np.asarray(self.Persistence_nsigmas)

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
		print len(Filaments_to_be_deleted)

		self.xdimPos = np.delete(self.xdimPos, Filaments_to_be_deleted)
		self.ydimPos = np.delete(self.ydimPos, Filaments_to_be_deleted)
		self.zdimPos = np.delete(self.zdimPos, Filaments_to_be_deleted)
		self.FilamentPos = np.delete(self.FilamentPos, Filaments_to_be_deleted)

		self.NFils = len(self.xdimPos)

	def NumParticles_per_filament(self):
		""" 
		Checks the number of dark matter filaments per filament.
		The allowed distance between a particle and filament is arbitrarily chosen.
		"""
		time_start = time.clock()
		print 'Checking number of particles within each filament'
		self.Particles_per_filament = []

		TempPartPosX = self.PartPosX
		TempPartPosY = self.PartPosY
		TempPartPosZ = self.PartPosZ

		DistanceThreshold = 0.001*(self.xmax - self.xmin)
		DM_points = np.dstack((self.PartPosX.ravel(), self.PartPosY.ravel(), self.PartPosZ.ravel()))
		DM_tree = spatial.KDTree(DM_points[0])

		for i in range(len(self.zdimMasked)):
			PartCount = 0
			Fil_pts = np.dstack((np.array(self.MaskedFilamentSegments[i])[:,0].ravel(), np.array(self.MaskedFilamentSegments[i])[:,1].ravel(), self.zdimMasked[i]))
			Fil_tree = spatial.KDTree(Fil_pts[0])
			nearest_pts = DM_tree.query_ball_tree(Fil_tree, DistanceThreshold)
			for j in nearest_pts:
				if not len(j) == 0:
					PartCount += len(j)
			self.Particles_per_filament.append(PartCount)

		"""
		for i in range(len(self.zdimMasked)):
			PartCount = 0
			for j in range(len(self.zdimMasked[i])):
				index = []
				for k in range(len(TempPartPosX)):
					PointDistance = np.sqrt((self.MaskedFilamentSegments[i][j][0] - TempPartPosX[k])**2 + (self.MaskedFilamentSegments[i][j][1] - TempPartPosY[k])**2 \
									+ (self.zdimMasked[i][j] - TempPartPosZ[k])**2)
					if PointDistance < DistanceThreshold:
						PartCount += 1
						index.append(k)
						#TempPartPosX = np.delete(TempPartPosX, k)
						#TempPartPosY = np.delete(TempPartPosY, k)
						#TempPartPosZ = np.delete(TempPartPosZ, k)
					else:
						continue
				#TempPartPosX = np.delete(TempPartPosX, index)
				#TempPartPosY = np.delete(TempPartPosY, index)
				#TempPartPosZ = np.delete(TempPartPosZ, index)
			self.Particles_per_filament.append(PartCount)
		"""

		self.Particles_per_filament = np.asarray(self.Particles_per_filament)
		print 'Number particle per filament check time: ', time.clock() - time_start, 's'

	def NumParticles_per_filament_v2(self):
		""" 
		Checks the number of dark matter filaments per filament.
		The allowed distance between a particle and filament is arbitrarily chosen.
		"""
		time_start = time.clock()
		print 'Checking number of particles within each filament'
		self.Particles_per_filament = []

		def Find_distance(FilamentPos1, FilamentPos2, PartNeighbourIDs):
			""" 
			Computes the distances of the particle to the filament
			Order of the distance array corresponds to the ID of the particle
			"""
			Particle_distances = []
			particle_x = []
			particle_y = []
			particle_z = []
			for id1 in PartNeighbourID1:
				particle_x.append(PartPosX[id1])
				particle_y.append(PartPosY[id1])
				particle_z.append(PartPosZ[id1])
			
			ParticlePoints = np.dstack((np.array(particle_x).ravel(), np.array(particle_y).ravel(), np.array(particle_z).ravel()))[0]
			SegmentLength = np.linalg.norm(np.concatenate([FilamentPos1, FilamentPos2]))
			Distance1 = []
			Distance2 = []
			for PartCoord in ParticlePoints:
				Distance1.append(np.linalg.norm(np.concatenate([FilamentPos1, PartCoord])))
				Distance2.append(np.linalg.norm(np.concatenate([FilamentPos2, PartCoord])))
			
			Angle = np.arccos((SegmentLength**2 + np.array(Distance1)**2 - np.array(Distance2)**2)/(2*SegmentLength*np.array(Distance2)))
			Distance = np.array(Distance2)*np.sin(Angle)
			Particle_distances.append(Distance)
			return np.array(Particle_distances)

		BoxSize = self.xmax-self.xmin
		DistanceThreshold = 0.001*BoxSize
	
		DuplicateCount_array = np.histogram(self.FilamentIDs, bins=np.arange(1, self.NFils+1))[0]
		if not MaskXdir or not MaskYdir or MaskZdir:
			a = 1
			#while i < len(self.FilamentIDs)
		else:
			# Computes for all filaments
			while i < len(self.FilamentIDs):
				DuplicateCount = DuplicateCount_array[i]
				if DuplicateCount == 0:
					FilamentPoints = np.dstack((self.xdimPos[i].ravel(), self.ydimPos[i].ravel(), self.zdimPos[i].ravel()))
				else:
					xTemp = xdimPos[i]
					yTemp = ydimPos[i]
					zTemp = zdimPos[i]
					for k in range(1,DuplicateCount+1):
						xTemp = np.concatenate([xTemp, xdimPos[i+k]])
						yTemp = np.concatenate([yTemp, ydimPos[i+k]])
						zTemp = np.concatenate([zTemp, zdimPos[i+k]])
					FilamentPoints = np.dstack((self.xTemp.ravel(), self.yTemp.ravel(), self.zTemp.ravel()))

				Neighbours_indices = DM_KDTree.query_ball_point(FilamentPoints[0], BoxSize/2.0)
				Neighbours_indices = np.unique(np.concatenate(Neighbours_indices, axis=0))
				Accepted_IDs = np.array([])

				for j in range(len(FilamentPoints)-1):
					Distances = Find_distance(FilamentPoints[j], FilamentPoints[j+1], Neighbours_indices)
					Distance_masking = Distances >= DistanceThreshold
					Accepted_IDs = np.concatenate([Accepted_IDs, Neighbours_indices[Distance_masking]])

				self.Particles_per_filament.append(len(np.unique(Accepted_IDs)))
				i += DuplicateCount + 1
				
		print 'Number particle per filament check time: ', time.clock() - time_start, 's'

	def Number_filament_connections(self):
		""" Computes the number of filament connections one filament has"""
		self.NumFilamentConnections = []
		for Connected_CP in self.PairIDS:
			counter = 0
			for CP_ids in Connected_CP:
				counter += self.Number_filaments_connecting_to_CP[CP_ids]
			self.NumFilamentConnections.append(counter)

		self.NumFilamentConnections = np.asarray(self.NumFilamentConnections)

	def Plot_Figures(self, filename, ndim=3):
		""" All plots done in this function	"""
		print 'Plotting'
		if Comparison:
			#NormedArgument = True
			NormedArgument = False
		else:
			NormedArgument = False

		if FilamentColors:
			ColorArray = np.linspace(0,1,110)
			ColorMap2D = np.asarray([np.mean(zvalues) for zvalues in self.zdimPos])
			ColorMapLength = np.asarray([lengths for lengths in self.LengthSplitFilament])
			if IncludeSlicing:
				ColorMap2DCutOff = np.asarray([np.mean(zvalues) for zvalues in self.CutOffzDim])
				ColorMap2DMasked = np.asarray([np.mean(zvalues) for zvalues in self.zdimMasked])
				ColorMapLengthCutOff = np.asarray([lengths for lengths in self.CutOffLengths])
				ColorMapLengthMasked = np.asarray([lengths for lengths in self.MaskedLengths])
		else:
			ColorArray = None
		
		if IncludeUnits:
			LegendText = ' - [Mpc/h]'
		else:
			LegendText = ''

		if HistogramPlots:
			# Histogram of number of filament connections
			###### NEEDS FIXING!!
			NumMin = min(self.NumFilamentConnections)
			NumMax = max(self.NumFilamentConnections)
			BinSize = (NumMax - NumMin)/(0.5) + 1
			Bin_list = np.linspace(NumMin, NumMax, BinSize)
			ConnectedHist = plt.figure()
			plt.hist(self.NumFilamentConnections, align='left', rwidth=1, bins=Bin_list, normed=NormedArgument)
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
			plt.hist(self.FilLengths, align='left', rwidth=1, bins=600, normed=NormedArgument, histtype='step')#, bins=BinList_FilLengths)
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
			plt.hist(self.NFilamentPoints, align='left', rwidth=1, bins=BinList_FilPts, normed=NormedArgument)
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
				plt.hist(self.Particles_per_filament, align='left', rwidth=1, bins=600, normed=NormedArgument)
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
					line_segmentsDMlen = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMapLengthCutOff, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsDMlen)
					DMParticleHistwFilamentsLengthCbar.colorbar(line_segmentsDMlen)
					plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
					plt.xlabel('$\mathregular{x}$' + LegendText)
					plt.ylabel('$\mathregular{y}$' + LegendText)
					plt.hold(False)
					plt.title('Dark matter density field over a segment of the particle box. \n Includes filaments with colorbar.'\
							 +'Colors indicate length of a filament.')

					# Interpolated smoothed out 2D histogram of dark matter particles
					Interpolated_DM_particles_figure = plt.figure()
					ax = Interpolated_DM_particles_figure.add_subplot(111, title='Dark matter particle density field, interpolated', aspect='equal',\
									xlim=DMBinXedges[[0,-1]], ylim=DMBinYedges[[0,-1]])
					im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
					xcenters = (DMBinXedges[:-1] + DMBinYedges[1:])/2.0
					ycenters = (DMBinYedges[:-1] + DMBinYedges[1:])/2.0
					im.set_data(xcenters, ycenters, DMHistogram)
					ax.images.append(im)


		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			if HistogramPlots:
				ConnectedHist.savefig(self.results_dir + 'NumberFilamentConnectedHistogram' + self.ModelFilename)
				FilamentLengthsHist.savefig(self.results_dir + 'FilamentLengthsHistogram' + self.ModelFilename)
				FilamentPtsHis.savefig(self.results_dir + 'FilamentPointsHistogram' + self.ModelFilename)
				#if IncludeDMParticles:
				#	NumParticlesFilamentHist.savefig(self.results_dir + 'NumParticlesPerFilament' + self.ModelFilename)
			if PlotFilaments:
				FilPositions.savefig(self.results_dir + 'FilamentPositions' + self.ModelFilename)
			if PlotFilamentsWCritPts:
				FilPositions_WCritPts.savefig(self.results_dir + 'FilamentPositionsWithCritPts' + self.ModelFilename)
			if Projection2D:
				FilPositions_2DProjection.savefig(self.results_dir + '2DProjectedFilamentPosition' + self.ModelFilename)
				if ColorBarZDir:
					FilPositions_2DProjectionColorBarZDir.savefig(self.results_dir + '2DProjectionColorBarZDir' + self.ModelFilename)
				if ColorBarLength:
					FilPositions_2DProjectionColobarLength.savefig(self.results_dir + '2DProjectionColobarLength' + self.ModelFilename)
			if IncludeSlicing and PlotFilaments:
				FilamentSliced.savefig(self.results_dir + 'Sliced3dBox' + self.ModelFilename)
				FilamentCutOff.savefig(self.results_dir + 'CutOffFilaments' + self.ModelFilename)
				if Projection2D:
					FilPositions_2DProjectionSliced.savefig(self.results_dir + '2DProjectionSliced3dBox' + self.ModelFilename)
				if ColorBarZDir:
					FilPositions_2DSlicedProjectionColorBarZDir.savefig(self.results_dir + '2DProjectionSlicedColobarZDir' + self.ModelFilename)
				if ColorBarLength:
					FilPositions_2DSlicedProjectionColorBarLen.savefig(self.results_dir + '2DProjectionSlicedColobarLength' + self.ModelFilename)
					FilPositions_2DSlicedProjection_Cbar_CritPts.savefig(self.results_dir + '2DProjectionSliced_CbarLength_CritPts' + self.ModelFilename)
				if IncludeDMParticles:
					if MaskXdir == 0 and MaskYdir == 0 and MaskZdir == 1:
						DMParticleHist.savefig(self.results_dir + 'DMParticleHistogram_ZMasked' + self.ModelFilename)
						DMParticleHistwFilaments.savefig(self.results_dir + 'DMParticleHistogramWFIlaments_ZMasked' + self.ModelFilename)
						ONEDHistX.savefig(self.results_dir + 'DMParticle1DHistogramXposition' + self.ModelFilename)
						DMParticleHistwFilamentsLengthCbar.savefig(self.results_dir + 'DMParticleHistogramWFilaments_LengthCbar_ZMasked' + self.ModelFilename)
						Interpolated_DM_particles_figure.savefig(self.results_dir + 'DMPartickeHistogram_interpolated' + self.ModelFilename)
					if MaskXdir == 1 and MaskYdir == 1 and MaskZdir == 1:
						DMParticleHist.savefig(self.results_dir + 'DMParticleHistogram_XYZMasked' + self.ModelFilename)
						DMParticleHistwFilaments.savefig(self.results_dir + 'DMParticleHistogramWFIlaments_XYZMasked' + self.ModelFilename)
						DMParticleHistwFilamentsLengthCbar.savefig(self.results_dir + 'DMParticleHistogramWFilaments_LengthCbar_XYZMasked' + self.ModelFilename)

		else:
			plt.show()

		plt.close('all')

	def Solve(self, filename, ndim=3, Sigma_threshold = False, robustness=False):
		""" 
		Runs the whole thing.
		Creates a pickle file of certain data unless it already exist.
		Removes the old pickle files if any of the mask directions are changed.
		"""
		# Determines folder name of pickle files
		cachedir_foldername_extra = 'npart'+str(self.nPart)
		if self.SigmaArg:
			cachedir_foldername_extra += 'nsig'+str(self.SigmaArg)
		else:
			cachedir_foldername_extra += 'nsig3'

		if robustness:
			cachedir_foldername_extra += 'TRIM'

		if HOMEPC == 0:
			cachedir='/PythonCaches/Disperse_analysis/'+cachedir_foldername_extra+'/'
		else:
			#cachedir='/mn/stornext/u3/aleh/Masters_project/PythonCaches/Disperse_analysis/' + cachedir_foldername_extra + '/'
			cachedir = '/mn/stornext/d13/euclid/aleh/PythonCaches/Disperse_analysis/' + cachedir_foldername_extra + '/'
		if not os.path.isdir(cachedir):
			os.makedirs(cachedir)
		# Pickle filenames and folder directory
		Boundary_check_cachefn = cachedir + "check_boundary_compact.p"
		Mask_slice_cachefn = cachedir + "mask_slice.p"
		Pickle_check_fn = cachedir + 'masking_check.p'
		
		self.ReadFile(filename, ndim)
		self.Sort_filament_coordinates(ndim)
		self.Sort_filament_data()
		self.Number_filament_connections()
		# Get masked critical points
		self.MaskedXCP, self.MaskedYCP, self.MaskedZCP = MaskCritPts.Mask_CPs(self.CritPointXpos, self.CritPointYpos, self.CritPointZpos,\
																 Mask_boundary_list, Mask_direction_check)

		if Sigma_threshold:
			# Filters filament based on a given sigma threshold. Still work in progress
			self.Filter_filaments(Sigma_threshold)

		if os.path.isfile(Pickle_check_fn):
			"""
			Check if pickle file exist. If it exists, it will read from the existing pickle file and get data from there.
			The pickle file will check if the masking conditions are changed or not. If it is changed, will recalculate stuff.
			If pickle file does not exist, then the program will calculate stuff. 
			"""
			Pickle_check = pickle.load(open(Pickle_check_fn, 'rb'))
			if not MaskXdir == Pickle_check[0] or not MaskYdir == Pickle_check[1] or not MaskZdir == Pickle_check[2]\
				 or not Sigma_threshold == Pickle_check[3]:	
				print 'Mask directions changed. Removing old pickle files'
				os.remove(Boundary_check_cachefn)
				os.remove(Mask_slice_cachefn)
				pickle.dump([MaskXdir, MaskYdir, MaskZdir, Sigma_threshold], open(Pickle_check_fn,'wb'))
		else:
			pickle.dump([MaskXdir, MaskYdir, MaskZdir, Sigma_threshold], open(Pickle_check_fn,'wb'))

		if os.path.isfile(Boundary_check_cachefn):
			print "reading from boundary check pickle file..."
			BC_instance_variables = pickle.load(open(Boundary_check_cachefn, 'rb'))
		else:
			BC_instance = BoundaryChecker.BoundaryChecker(self.xmin, self.xmax, self.xdimPos, self.ydimPos, self.zdimPos, self.FilID, self.NFils)
			BC_instance_variables = BC_instance.Get_periodic_boundary()
			pickle.dump(BC_instance_variables, open(Boundary_check_cachefn, 'wb'))
		self.FilamentIDs, self.FilamentPos, self.xdimPos, self.ydimPos, self.zdimPos, self.LengthSplitFilament, self.FilLengths = BC_instance_variables

		if IncludeSlicing:
			if os.path.isfile(Mask_slice_cachefn):
				print "reading from mask_slice pickle file..."
				Mask_instance_variables = pickle.load(open(Mask_slice_cachefn, 'rb'))
			else:
				Mask_instance = FilamentMasking.FilamentMasking(self.FilamentPos, self.xdimPos, self.ydimPos, self.zdimPos,\
												self.LengthSplitFilament , self.NFils, Mask_direction_check, Mask_boundary_list)
				Mask_instance_variables = Mask_instance.Mask_slices()
				pickle.dump(Mask_instance_variables, open(Mask_slice_cachefn, 'wb'))
			self.MaskedFilamentSegments, self.MaskedLengths, self.zdimMasked, self.CutOffFilamentSegments, self.CutOffLengths, self.CutOffzDim = Mask_instance_variables
		
		#if IncludeSlicing and IncludeDMParticles and not Comparison:
		#	self.NumParticles_per_filament_v2()
		"""
		if IncludeSlicing and IncludeDMParticles and not Comparison:
			self.Read_SolveFile()
			#self.Mask_DMParticles()
			#self.NumParticles_per_filament()
		"""
		
		if not Comparison:
			if self.savefile == 2:
				print 'Done! No files saved.'
			else:
				self.Plot_Figures(filename, ndim)
		
		return self.NumFilamentConnections, sorted(self.FilLengths), self.NFilamentPoints

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


class Histogram_Comparison():
	def __init__(self, savefile, savefigDirectory, redshift, LCDM=False, SymmA=False, SymmB=False, nsigComparison=False):
		self.savefile = savefile
		self.LCDM_check = LCDM
		self.SymmA_check = SymmA
		self.SymmB_check = SymmB

		self.ParticleComparison = False
		self.ModelComparison = False
		self.SigmaComparison = False

		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir):
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
			self.ModelFilename += filetype
		elif nPart == 64:
			self.ModelFilename += '64Part' + filetype
		elif nPart == 128:
			self.ModelFilename += '128Part' + filetype
		elif nPart == 256:
			self.ModelFilename += '256Part' + filetype
		elif nPart == 512:
			self.ModelFilename += '512Part' + filetype
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
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', '$\mathregular{256^3}$ particles', '$\mathregular{512^3}$ particles']
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
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of number filament connections \n with '+str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)		
		plt.legend(self.LegendText)
		plt.hold(False)
	
		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(self.N):
			plt.hist(self.FilamentLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of filament length \n with' +str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)
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
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of number of datapoints per filament \n with' +str(self.nParticles) + '$\mathregular{^3}$ particle'+subsample_text)
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
		Legends = ['$\mathregular{\sigma=5}, 64^3$ part subsample', '$\mathregular{\sigma=4}, 64^3$ part subsample', '$\mathregular{\sigma=3}, 64^3$ part subsample', \
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
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of number connections per filament')
		plt.legend(Legends)
		plt.hold(False)

		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			plt.hist(FilLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
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
		plt.ylabel('Number of occurances')
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
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of number connections per filament for $\mathregular{\sigma=}$' + str(nsigma))
		plt.legend(Legends)
		plt.hold(False)

		LengthHistComparison = plt.figure()
		plt.hold(True)
		for i in range(N):
			plt.hist(FilLengths[i], align='mid', rwidth=1, bins=400, normed=False, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
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
		plt.ylabel('Number of occurances')
		plt.title('Histogram comparison of number data points per filament for $\mathregular{\sigma=}$' + str(nsigma))
		plt.legend(Legends)
		plt.hold(False)

		if self.savefile == 1:
			print '--- SAVING IN: ', self.results_dir, ' ---'
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison_AllParticles_nsig' + str(nsigma) + filetype)
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison_AllParticles_nsig' + str(nsigma) + filetype)
			NPointsHistComparison.savefig(self.results_dir + 'HistNPointsComparison_AllParticles_nsig' + str(nsigma) + filetype)
		elif self.savefile == 2:
			print 'Done! No histograms to plot.'
		else:
			plt.show()
		plt.close('all')

if __name__ == '__main__':
	HOMEPC = 1					# Set 1 if working in UiO terminal

	# Filament and dark matter particle plotting
	FilamentLimit = 0			# Limits the number of lines read from file. Reads all if 0
	PlotFilaments = 1			# Set 1 to plot actual filaments
	PlotFilamentsWCritPts = 0	# Set to 1 to plot filaments with critical points
	Projection2D = 0			# Set to 1 to plot a 2D projection of the 3D case
	FilamentColors = 1 			# Set to 1 to get different colors for different filaments
	ColorBarZDir = 1 			# Set 1 to include colorbar for z-direction
	ColorBarLength = 1 			# Set 1 to include colorbars based on length of the filament
	IncludeDMParticles = 1 		# Set to 1 to include dark matter particle plots
	IncludeSlicing = 1			# Set 1 to include slices of the box
	MaskXdir = 0 				# Set 1 to mask one or more directions.
	MaskYdir = 0
	MaskZdir = 1

	# Histogram plots
	HistogramPlots = 1			# Set to 1 to plot histograms
	Comparison = 0				# Set 1 if you want to compare different number of particles. Usual plots will not be plotted!
	ModelCompare = 0 			# Set to 1 to compare histograms of different models. Particle comparisons will not be run.
	SigmaComparison = 0 		# Set to 1 to compare histograms and/or plots based on different sigma values by MSE.
								# Must also set Comparison=1 to compare histograms
	
	# Run simulation for different models. Set to 1 to run them. 
	LCDM_model = 1 
	SymmA_model = 0
	SymmB_model = 0
		
	# Global properties to be set
	IncludeUnits = 1			# Set to 1 to include 'rockstar' units, i.e Mpc/h and km/s
	SaveAsPNG = 1				# Set 1 to save figures as PNG
	SaveAsPDF = 0 				# Set 1 to save figures as PDF

	# Some if tests before the simulation runs
	if FilamentLimit == 1:
		print 'Running program with limited amount of filaments.'

	if IncludeDMParticles == 1 and HOMEPC == 1:
		print 'Dark matter particles will be included'
		IncludeSlicing = 1
	else:
		print 'Dark matter particles not included'

	if IncludeSlicing == 1 and MaskXdir == 0 and MaskYdir == 0 and MaskZdir == 0:
		raise ValueError('IncludeSlicing set to 1, but all mask direction set to 0. Set at least one of them to 1.')

	if ModelCompare == 1:
		print 'Will compare histograms for different models'
	elif Comparison == 1 and ModelCompare == 0:
		print 'Will compare histograms over different models or number of particles.'

	UnitConverter = 256.0 if IncludeUnits == 1 else 1
	### Slice selection can be chosen here
	if IncludeSlicing == 1:
		print 'Slicing included'
		LowerBoundaryXDir = 0
		LowerBoundaryYDir = 0
		LowerBoundaryZDir = 0
		UpperBoundaryXDir = 0
		UpperBoundaryYDir = 0
		UpperBoundaryZDir = 0
		if MaskXdir == 1:
			LowerBoundaryXDir = 0.45*UnitConverter
			UpperBoundaryXDir = 0.55*UnitConverter
			print 'Masking X direction'
		if MaskYdir == 1:
			LowerBoundaryYDir = 0.45*UnitConverter
			UpperBoundaryYDir = 0.55*UnitConverter
			print 'Masking Y direction'
		if MaskZdir == 1:
			LowerBoundaryZDir = 0.45*UnitConverter
			UpperBoundaryZDir = 0.55*UnitConverter
			print 'Masking Z direction'

	Mask_boundary_list = [UpperBoundaryXDir, UpperBoundaryYDir, UpperBoundaryZDir, LowerBoundaryXDir, LowerBoundaryYDir, LowerBoundaryZDir]
	Mask_direction_check = [MaskXdir, MaskYdir, MaskZdir]

	if SaveAsPNG == 1 and SaveAsPDF	== 1:
		raise ValueError('Cannot save both PDF and PNG at the same time. Only allow one at a time.')
	elif SaveAsPDF == 1 and SaveAsPNG == 0:
		print 'Saving figures as PDF files'
		filetype = '.pdf'
	elif SaveAsPNG == 1 and SaveAsPDF == 0:
		print 'Saving figures as PNG files'
		filetype = '.png'
	else:
		raise ValueError('Figure filetype to save not selected.')

	
	if HOMEPC == 0:
		file_directory = 'C:/Users/Alex/Documents/Masters_project/Disperse'
		savefile_directory = file_directory
		IncludeDMParticles = 0
		if LCDM_model == 1:
			print '=== Running for the LCDM model ==='
			#solveInstance1 = Disperse_Plotter(savefile=1, savefigDirectory='Plot_Disperse_Example/', nPart=64)
			#solveInstance1.Plot("simu_2D.ND.NDnet_s3.up.NDskl.a.NDskl", ndim=2)
			#solveInstance1.Plot("simu_32_id.gad.NDnet_s3.5.up.NDskl.a.NDskl", ndim=3)

			LCDM_z0_64_dir = 'lcdm_z0_testing/LCDM_z0_64PeriodicTesting/'
			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'PlotsTest/', nPart=64, model='LCDM', redshift=0)
			NConn_64PartLCDM, FilLen_64PartLCDM, NPts_64PartLCDM = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl', Sigma_threshold=4.0)

			#LCDM_z0_128_dir = 'lcdm_z0_testing/LCDM_z0_128PeriodicTesting/'
			#LCDM_z0_128Instance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128, model='LCDM', redshift=0)
			#NConn_128PartLCDM, FilLen_128PartLCDM, NPts_128PartLCDM = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl')
			
			#LCDM_nsig4Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Sigma4PlotsTest/', nPart=64, model='LCDM', redshift=0, SigmaArg=4)
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
			SymmA_z064_instance = Disperse_Plotter(savefile=1, savefigDirectory=SymmA_z064_directory+'Plots/', nPart=64, model='SymmA', redshift=0)
			Nconn_64PartSymmA, FilLen_64PartSymmA, NPts_64PartSymmA = SymmA_z064_instance.Solve(SymmA_z064_directory+'SkelconvOutput_SymmAz064Part.a.NDskl')

		if SymmB_model:
			print '=== Running for Symm_B model ==='
			SymmB_z064_directory = 'SymmB_data/SymmB_z0_64Particles/'
			SymmB_z064_instance = Disperse_Plotter(savefile=1, savefigDirectory=SymmB_z064_directory+'Plots/', nPart=64, model='SymmB', redshift=0)
			Nconn_64PartSymmB, FilLen_64PartSymmB, NPts_64PartSymmB = SymmB_z064_instance.Solve(SymmB_z064_directory+'SkelconvOutput_SymmBz064Part.a.NDskl')

		if ModelCompare:
			if LCDM and not SymmA and not SymmB:
				raise ValueError('Only LCDM model is ticked on. Cannot compare models if SymmA and/or SymmB is ticked as well')

			elif not LCDM and SymmA and not SymmB:
				raise ValueError('Only SymmA model is ticked on. Cannot compare models if LCDM and/or SymmB is ticked as well')

			elif not LCDM and not SymmA and SymmB:
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

			ModelCompareInstance = Histogram_Comparison(savefile=1, savefigDirectory=ModelComparisonDir+ModelsFolder+'Plots/', redshift=0,\
							LCDM=LCDM_model, SymmA=SymmA_model, SymmB=SymmB_model)
			ModelCompareInstance.Run(NumConnectionsModels_list, FilLengthsModels_list, FilPointsModels_list, 64)

	if HOMEPC == 1:
		file_directory = '/mn/stornext/d5/aleh'
		savefile_directory = '/mn/stornext/u3/aleh/Masters_project/disperse_results'
		
		if LCDM_model == 1:
			print '=== Running for the LCDM model ==='
			solve_file_dir = '/lcdm_testing/'
			solve_filename = 'lcdm_z0_test.solve'
			
			if IncludeDMParticles == 1:
				Gadget_instance = ReadGadgetFile.Read_Gadget_file(Mask_direction_check, Mask_boundary_list)
				PartPosX, PartPosY, PartPosZ, DMHistogram, DMBinXedges, DMBinYedges, DM_KDTree = Gadget_instance.Get_particles(includeKDTree=False)
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
			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plotstest2/', nPart=64, model='LCDM', redshift=0)
			NumConn_64LCDM, FilLen_64LCDM, NPts_64LCDM = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl')
			"""
			LCDM_z0_128_dir = 'lcdm_testing/LCDM_z0_128PeriodicTesting/'
			LCDM_z0_128Instance = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128, model='LCDM', redshift=0)
			NumConn_128LCDM, FilLen_128LCDM, NPts_128LCDM = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl')
			
			LCDM_z0_256_dir = 'lcdm_testing/LCDM_z0_256PeriodicTesting/'
			LCDM_z0_256Instance = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_256_dir+'Plots/', nPart=256, model='LCDM', redshift=0)
			NumConn_256LCDM, FilLen_256LCDM, NPts_256LCDM = LCDM_z0_256Instance.Solve(LCDM_z0_256_dir+'SkelconvOutput_LCDMz0256.a.NDskl')
			
			LCDM_z0_512_dir = 'lcdm_testing/LCDM_z0_512PeriodicTesting/'
			LCDM_z0_512Instance = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_512_dir+'Plots/', nPart=512, model='LCDM', redshift=0)
			NumConn_512LCDM, FilLen_512LCDM, NPts_512LCDM = LCDM_z0_512Instance.Solve(LCDM_z0_512_dir+'SkelconvOutput_LCDMz0512.a.NDskl')
			"""

			#LCDM_z0_64_Robustness_dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/Robustness_argument/'
			#LCDM_z0_64_RobustnessInstance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_Robustness_dir+'Plots_Robustness/', nPart=64, model='LCDM', redshift=0)
			#NumConn_64LCDM_rob, FilLen_64LCDM_rob, NPts_64LCDM_rob = LCDM_z0_64_RobustnessInstance.Solve(LCDM_z0_64_Robustness_dir+'SkelconvOutput_LCDMz064_robustness3.TRIM.a.NDskl', robustness=True)

			if SigmaComparison:
				LCDM64_instance_nsig4 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_64_dir+'Sigma4/', nPart=64, model='LCDM', redshift=0, SigmaArg=4)
				LCDM64_instance_nsig5 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_64_dir+'Sigma5/', nPart=64, model='LCDM', redshift=0, SigmaArg=5)
				NConn_64nsig4, FilLen_64nsig4, NPts_64nsig4 = LCDM64_instance_nsig4.Solve(LCDM_z0_64_dir+'Sigma4/SkelconvOutput_LCDMz064_nsig4.a.NDskl')
				NConn_64nsig5, FilLen_64nsig5, NPts_64nsig5 = LCDM64_instance_nsig5.Solve(LCDM_z0_64_dir+'Sigma5/SkelconvOutput_LCDMz064_nsig5.a.NDskl')

				LCDM128_instance_nsig4 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_128_dir+'Sigma4/', nPart=128, model='LCDM', redshift=0, SigmaArg=4)
				LCDM128_instance_nsig5 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_128_dir+'Sigma5/', nPart=128, model='LCDM', redshift=0, SigmaArg=5)
				NConn_128nsig4, FilLen_128nsig4, NPts_128nsig4 = LCDM128_instance_nsig4.Solve(LCDM_z0_128_dir+'Sigma4/SkelconvOutput_LCDMz0128_nsig4.a.NDskl')
				NConn_128nsig5, FilLen_128nsig5, NPts_128nsig5 = LCDM128_instance_nsig5.Solve(LCDM_z0_128_dir+'Sigma5/SkelconvOutput_LCDMz0128_nsig5.a.NDskl')

				LCDM256_instance_nsig4 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_256_dir+'Sigma4/', nPart=256, model='LCDM', redshift=0, SigmaArg=4)
				LCDM256_instance_nsig5 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_256_dir+'Sigma5/', nPart=256, model='LCDM', redshift=0, SigmaArg=5)
				NConn_256nsig4, FilLen_256nsig4, NPts_256nsig4 = LCDM256_instance_nsig4.Solve(LCDM_z0_256_dir+'Sigma4/SkelconvOutput_LCDMz0256_nsig4.a.NDskl')
				NConn_256nsig5, FilLen_256nsig5, NPts_256nsig5 = LCDM256_instance_nsig5.Solve(LCDM_z0_256_dir+'Sigma5/SkelconvOutput_LCDMz0256_nsig5.a.NDskl')

				LCDM512_instance_nsig4 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_512_dir+'Sigma4/', nPart=512, model='LCDM', redshift=0, SigmaArg=4)
				LCDM512_instance_nsig5 = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_512_dir+'Sigma5/', nPart=512, model='LCDM', redshift=0, SigmaArg=5)
				NConn_512nsig4, FilLen_512nsig4, NPts_512nsig4 = LCDM512_instance_nsig4.Solve(LCDM_z0_512_dir+'Sigma4/SkelconvOutput_LCDMz0512_nsig4.a.NDskl')
				NConn_512nsig5, FilLen_512nsig5, NPts_512nsig5 = LCDM512_instance_nsig5.Solve(LCDM_z0_512_dir+'Sigma5/SkelconvOutput_LCDMz0512_nsig5.a.NDskl')
				
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

					ComparisonInstance_LCDM = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir+'SigmaComparisons/',\
										 redshift=0, LCDM=1, nsigComparison=1)
					ComparisonInstance_LCDM.Convergence_tests(NumConnections_list_expanded, FilLengths_list_expanded, FilPoints_list_expanded)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig3, FilLengths_list_nsig3, FilPoints_list_nsig3, 3)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig4, FilLengths_list_nsig4, FilPoints_list_nsig5, 4)
					ComparisonInstance_LCDM.Sigma_plot_comparison(Nconnections_list_nsig5, FilLengths_list_nsig5, FilPoints_list_nsig4, 5)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=512)
				else:
					NumConnections_list = [NumConn_64LCDM, NumConn_128LCDM , NumConn_256LCDM, NumConn_512LCDM]
					FilLengths_list = [FilLen_64LCDM, FilLen_128LCDM, FilLen_256LCDM, FilLen_512LCDM]
					FilPoints_list = [NPts_64LCDM, NPts_128LCDM, NPts_256LCDM, NPts_512LCDM]
					ComparisonInstance_LCDM = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, redshift=0, LCDM=1)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list)