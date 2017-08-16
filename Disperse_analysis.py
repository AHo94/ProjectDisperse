### Set comment on the two below to plot. Only use if running on papsukal, nekkar etc. 
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import os
#import scipy.stats as stats
from matplotlib import colors as mcolors
#from scipy import interpolate
import time
from scipy import spatial

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

	def ReadFile(self, filename, dimensions):
		""" Reads the data from the skeleton file from Disperse """
		time_start = time.clock()
		print 'Reading data for the file: ', filename, '...' 
		datafiles = open(os.path.join(file_directory, filename), 'r')
		self.CriticalPoints = []
		self.Filaments = []

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
				break

		datafiles.close()
		self.CriticalPoints = np.asarray(self.CriticalPoints)
		self.Filaments = np.asarray(self.Filaments)
		print 'Time elapsed for filament data reading: ', time.clock() - time_start, 's'

	def Read_SolveFile(self):
		time_start = time.clock()
		print 'Reading data for the file: ', solve_file_dir, solve_filename, '. May take a while...'
		#datafiles = open(os.path.join(file_directory+solve_file_dir, solve_filename), 'r')
		self.PartPosX = []
		self.PartPosY = []
		self.PartPosZ = []
		
		PartVelX = []
		PartVelY = []
		PartVelZ = []
		self.XYPartPos = []
		with open(os.path.join(file_directory+solve_file_dir, solve_filename), 'r') as datafiles:
			next(datafiles)
			for line in datafiles:
				data_set = line.split()
				if MaskZdir and not MaskYdir and not MaskXdir:
					if float(data_set[2])*UnitConverter > LowerBoundaryZDir and float(data_set[2])*UnitConverter < UpperBoundaryZDir:
						self.PartPosX.append(float(data_set[0])*UnitConverter)
						self.PartPosY.append(float(data_set[1])*UnitConverter)
						self.PartPosZ.append(float(data_set[2])*UnitConverter)
				elif MaskXdir == 1 and MaskYdir == 1 and MaskZdir == 1:
					if float(data_set[2])*UnitConverter > LowerBoundaryZDir and float(data_set[2])*UnitConverter < UpperBoundaryZDir\
					and float(data_set[1])*UnitConverter > LowerBoundaryYDir and float(data_set[1])*UnitConverter < UpperBoundaryYDir\
					and float(data_set[0])*UnitConverter > LowerBoundaryXDir and float(data_set[0])*UnitConverter < UpperBoundaryXDir:
						self.PartPosX.append(float(data_set[0])*UnitConverter)
						self.PartPosY.append(float(data_set[1])*UnitConverter)
						self.PartPosZ.append(float(data_set[2])*UnitConverter)

		"""
		SkipFirstLine = 0
		for line in datafiles:
			data_set = line.split()
			if IncludeSlicing == 1:
				if SkipFirstLine == 0:
					SkipFirstLine = 1
				else:
					if float(data_set[2])*UnitConverter > LowerBoundaryZDir and float(data_set[2])*UnitConverter < UpperBoundaryZDir:
						#self.XYPartPos.append([float(data_set[0]), float(data_set[1])])
						self.PartPosX.append(float(data_set[0])*UnitConverter)
						self.PartPosY.append(float(data_set[1])*UnitConverter)
						#PartPosZ.append(float(data_set[2]))
						#PartVelX.append(float(data_set[3]))
						#PartVelY.append(float(data_set[4]))
						#PartVelZ.append(float(data_set[5]))
		"""
		self.PartPosX = np.asarray(self.PartPosX)
		self.PartPosY = np.asarray(self.PartPosY)
		self.PartPosZ = np.asarray(self.PartPosZ)
		datafiles.close()
		print 'Read solve file time: ', time.clock() - time_start, 's'

	def Sort_arrays(self, dimensions):
		""" 
		Sorts data to their respective arrays 
		Data to be sorted: Critical points, ID of filament and filament points
		"""
		time_start = time.clock()
		print 'Sorting data ...'
		if FilamentLimit == 0:
			self.NFils = int(self.Filaments[0][0])
		else:
			self.NFils = FilamentLimit

		# General data
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
		# Data for non-connected critical points
		self.CriticalPointPosition = []
		self.CritPointXpos = []
		self.CritPointYpos = []
		self.CritPointZpos = []
		self.CritPointXposNOTCON = []
		self.CritPointYposNOTCON = []
		self.CritPointZposNOTCON = []
		if dimensions == 2:
			for i in range(self.NcritPts):
				if int(self.CritPointInfo[i][0]) == 0:
					self.CritPointXposNOTCON.append(float(self.CritPointInfo[i][1]))
					self.CritPointYposNOTCON.append(float(self.CritPointInfo[i][2]))
				else:
					self.CritPointXpos.append(float(self.CritPointInfo[i][1]))
					self.CritPointYpos.append(float(self.CritPointInfo[i][2]))
		elif dimensions == 3:
			for i in range(self.NcritPts):
				if int(self.CritPointInfo[i][0]) == 0:
					self.CritPointXposNOTCON.append(float(self.CritPointInfo[i][1]))
					self.CritPointYposNOTCON.append(float(self.CritPointInfo[i][2]))
					self.CritPointZposNOTCON.append(float(self.CritPointInfo[i][3]))
				else:
					self.CritPointXpos.append(float(self.CritPointInfo[i][1]))
					self.CritPointYpos.append(float(self.CritPointInfo[i][2]))
					self.CritPointZpos.append(float(self.CritPointInfo[i][3]))
		"""

		# Filament positions etc
		self.FilamentPos = []
		self.FilID = []
		self.xdimPos = []
		self.ydimPos = []
		self.NFilamentPoints = []
		k = 1
		if dimensions == 2:
			for i in range(1, self.NFils):
				Filstuff = self.Filaments[k]	# Contains info on filament ID, num crit pts etc.
				TempPositions = []
				xtemp = []
				ytemp = []
				self.FilID.append(float(Filstuff[0]))
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
				self.FilID.append(float(Filstuff[0]))
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

				if FilamentLimit != 0:
					if k >= FilamentLimit:
						self.NFils = len(self.xdimPos)
						break

		self.FilamentPos = np.array(self.FilamentPos)
		self.NFilamentPoints = np.array(self.NFilamentPoints)
		self.xdimPos = np.array(self.xdimPos)
		self.ydimPos = np.array(self.ydimPos)
		if dimensions == 3:
			self.zdimPos = np.array(self.zdimPos)

		print 'Array sorting time: ', time.clock() - time_start, 's'

	def Mask_slices(self):
		"""
		Creates a mask to plot a slice of the filament box. Boundary of slice chosen arbitrarily.
		The masking includes filaments that are within the given boundary.
		Also includes masking where the filament are cut off outside the boundary. 
		"""
		time_start = time.clock()
		print 'Computing masks'
		self.zdimMasked = []
		self.MaskedFilamentSegments = []
		self.MaskedLengths = []

		self.CutOffFilamentSegments = []
		self.CutOffzDim = []
		self.CutOffLengths = []
		# Mask Z direction only
		if MaskZdir == 1 and MaskXdir == 0 and MaskYdir == 0:
			for i in range(len(self.zdimPos)):
				if np.any(np.array(self.zdimPos[i]) > LowerBoundaryZDir) and np.any(np.array(self.zdimPos[i]) < UpperBoundaryZDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.zdimPos)):
				Indices = np.where(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append(np.array([self.xdimPos[i][idx], self.ydimPos[i][idx]]))
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)

		# Mask X direction only	
		if MaskXdir == 1 and MaskYdir == 0 and MaskZdir == 0:
			for i in range(len(self.xdimPos)):
				if np.any(np.array(self.xdimPos[i]) > LowerBoundaryXDir) and np.any(np.array(self.xdimPos[i]) < UpperBoundaryXDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.xdimPos)):
				Indices = np.where(np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)
		# Mask Y direction only		
		if MaskYdir == 1 and MaskXdir == 0 and MaskZdir == 0:
			for i in range(len(self.ydimPos)):
				if np.any(np.array(self.ydimPos[i]) > LowerBoundaryYDir) and np.any(np.array(self.ydimPos[i]) < UpperBoundaryYDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.ydimPos)):
				Indices = np.where(np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir)))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)

		# Mask X and Y direction
		if MaskXdir == 1 and MaskYdir == 1 and MaskZdir == 0:
			for i in range(len(self.xdimPos)):
				if np.any(np.array(self.xdimPos[i]) > LowerBoundaryXDir) and np.any(np.array(self.xdimPos[i]) < UpperBoundaryXDir)\
				 and np.any(np.array(self.ydimPos[i]) > LowerBoundaryYDir) and np.any(np.array(self.ydimPos[i]) < UpperBoundaryYDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.xdimPos)):
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)),\
					np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir))))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)
		# Mask Y and Z direction
		if MaskXdir == 0 and MaskYdir == 1 and MaskZdir == 1:
			for i in range(len(self.zdimPos)):
				if np.any(np.array(self.zdimPos[i]) > LowerBoundaryZDir) and np.any(np.array(self.zdimPos[i]) < UpperBoundaryZDir)\
				 and np.any(np.array(self.ydimPos[i]) > LowerBoundaryYDir) and np.any(np.array(self.ydimPos[i]) < UpperBoundaryYDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.xdimPos)):
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)),\
					np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir))))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)
		# Mask X and Z direction
		if MaskXdir == 1 and MaskYdir == 0 and MaskZdir == 1:
			for i in range(len(self.zdimPos)):
				if np.any(np.array(self.zdimPos[i]) > LowerBoundaryZDir) and np.any(np.array(self.zdimPos[i]) < UpperBoundaryZDir)\
				 and np.any(np.array(self.xdimPos[i]) > LowerBoundaryXDir) and np.any(np.array(self.xdimPos[i]) < UpperBoundaryXDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.xdimPos)):
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)),\
					np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir))))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)
		# Mask all directions
		if MaskXdir == 1 and MaskYdir == 1 and MaskZdir == 1:
			for i in range(len(self.zdimPos)):
				if np.any(np.array(self.zdimPos[i]) > LowerBoundaryZDir) and np.any(np.array(self.zdimPos[i]) < UpperBoundaryZDir)\
				and np.any(np.array(self.ydimPos[i]) > LowerBoundaryYDir) and np.any(np.array(self.ydimPos[i]) < UpperBoundaryYDir)\
				and np.any(np.array(self.xdimPos[i]) > LowerBoundaryXDir) and np.any(np.array(self.xdimPos[i]) < UpperBoundaryXDir):
					self.zdimMasked.append(self.zdimPos[i])
					self.MaskedFilamentSegments.append(self.FilamentPos[i])
					self.MaskedLengths.append(self.LengthSplitFilament[i])

			for i in range(len(self.xdimPos)):
				Indices = np.where(
					np.logical_and(np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir)),
					np.logical_and(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)),\
					np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)))))[0]
				if not len(Indices) == 0:
					FilSegmentTemp = []
					zDimTemp = []
					for idx in Indices:
						zDimTemp.append(self.zdimPos[i][idx])
						FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
					self.CutOffFilamentSegments.append(FilSegmentTemp)
					self.CutOffzDim.append(zDimTemp)
					TempLen = 0
					for j in range(len(zDimTemp)-1):
						TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2+ (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
								 + (zDimTemp[j+1] - zDimTemp[j])**2)
					self.CutOffLengths.append(TempLen)

		self.MaskedFilamentSegments = np.asarray(self.MaskedFilamentSegments)
		self.MaskedLengths = np.asarray(self.MaskedLengths)
		self.zdimMasked = np.asarray(self.zdimMasked)
		self.CutOffFilamentSegments = np.asarray(self.CutOffFilamentSegments)
		self.CutOffLengths = np.asarray(self.CutOffLengths)
		self.CutOffzDim = np.asarray(self.CutOffzDim)
		print 'Masking time: ', time.clock() - time_start, 's'

	def Mask_DMParticles(self):
		""" Computes a mask for the dark matter particles """
		if MaskZdir and not MaskXdir and not MaskYdir:
			mask = np.logical_and(ParticlePos[:,2] < UpperBoundaryZDir, ParticlePos[:,2] > LowerBoundaryZDir)
		elif MaskXdir and not MaskYdir and not MaskZdir:
			mask = np.logical_and(ParticlePos[:,0] < UpperBoundaryZDir, ParticlePos[:,0] > LowerBoundaryZDir)
		elif MaskYdir and not MaskXdir and not MaskZdir:
			mask = np.logical_and(ParticlePos[:,1] < UpperBoundaryZDir, ParticlePos[:,1] > LowerBoundaryZDir)
		
		self.PartPosX = ParticlePos[mask,0]
		self.PartPosY = ParticlePos[mask,1]
		self.PartPosZ = ParticlePos[mask,2]
			
	def Check_boundary(self):
		"""
		This function checks whether a filament crosses the boundary or not. 
		If a filament crosses the boundary, the filament will be split into two/three filaments, i.e different lists/arrays.
		More details of the algorithm, see paper.
		Function also computes the length of the filament.
		"""
		time_start = time.clock()
		print 'Checking boundaries'
		self.BoxSize = self.xmax - self.xmin
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
		self.LengthSplitFilament = []

		def New_points(P1, P2, boundary_dir, boundary):
			""" Computes the points of the two non-boundary crossing coordinates at the boundary. """
			if boundary == 0:
				DifferenceAdd = -self.BoxSize
				BoxBoundary = self.xmin
			elif boundary == 1:
				DifferenceAdd = self.BoxSize
				BoxBoundary = self.xmax

			if boundary_dir == 'x':
				t_variable = (BoxBoundary - P1[0])/(P2[0] - P1[0] - self.BoxSize)
				y_coord = P1[1] + (P2[1] - P1[1])*t_variable
				z_coord = P1[2] + (P2[2] - P1[2])*t_variable
				return y_coord, z_coord
			elif boundary_dir == 'y':
				t_variable = (BoxBoundary - P1[1])/(P2[1] - P1[1] - self.BoxSize)
				x_coord = P1[0] + (P2[0] - P1[0])*t_variable
				z_coord = P1[2] + (P2[2] - P1[2])*t_variable
				return x_coord, z_coord
			elif boundary_dir == 'z':
				t_variable = (BoxBoundary - P1[2])/(P2[2] - P1[2] - self.BoxSize)
				x_coord = P1[0] + (P2[0] - P1[0])*t_variable
				y_coord = P1[1] + (P2[1] - P1[1])*t_variable
				return x_coord, y_coord
			else:
				raise ValueError('boundary_dir argument not set properly')

		for i in range(len(self.xdimPos)):
			SplitFilament = 0
			xBoundary = 0
			yBoundary = 0
			zBoundary = 0
			xyTemp = []
			xTemp = []
			yTemp = []
			zTemp = []
			xyNewTemp = []
			xNewTemp = []
			yNewTemp = []
			zNewTemp = []
			if SplitFilament == 2:
				xyNewTemp2 = []
				xNewTemp2 = []
				yNewTemp2 = []
				zNewTemp2 = []
			for j in range(len(self.xdimPos[i])-1):
				xDiff = np.abs(self.xdimPos[i][j+1] - self.xdimPos[i][j])
				yDiff = np.abs(self.ydimPos[i][j+1] - self.ydimPos[i][j])
				zDiff = np.abs(self.zdimPos[i][j+1] - self.zdimPos[i][j])

				if SplitFilament == 0:
					xyTemp.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xTemp.append(self.xdimPos[i][j])
					yTemp.append(self.ydimPos[i][j])
					zTemp.append(self.zdimPos[i][j])

				elif SplitFilament == 1:
					if xBoundary == 1 or yBoundary == 1 or zBoundary == 1:
						Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
						Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
					if xBoundary == 1:
						if np.abs(self.xdimPos[i][j-1] - self.xmin) < 0.5*self.BoxSize:
							Ypoint, Zpoint = New_points(Point1, Point2, 'x', 0)
							xyTemp.append(np.array([self.xmin, Ypoint]))
							xTemp.append(self.xmin)
							yTemp.append(Ypoint)
							zTemp.append(Zpoint)
							xyNewTemp.append(np.array([self.xmax, Ypoint]))
							xNewTemp.append(self.xmax)
							yNewTemp.append(Ypoint)
							zNewTemp.append(Zpoint)			
							xBoundary = 2
						elif np.abs(self.xdimPos[i][j-1] - self.xmax) < 0.5*self.BoxSize:
							Ypoint, Zpoint = New_points(Point1, Point2, 'x', 1)
							xyTemp.append(np.array([self.xmax, Ypoint]))
							xTemp.append(self.xmax)
							yTemp.append(Ypoint)
							zTemp.append(Zpoint)
							xyNewTemp.append(np.array([self.xmin, Ypoint]))
							xNewTemp.append(self.xmin)
							yNewTemp.append(Ypoint)
							zNewTemp.append(Zpoint)			
							xBoundary = 2
					elif yBoundary == 1:
						if np.abs(self.ydimPos[i][j-1] - self.ymin) < 0.5*self.BoxSize:
							Xpoint, Zpoint = New_points(Point1, Point2, 'y', 0)
							xyTemp.append(np.array([Xpoint, self.ymin]))
							xTemp.append(Xpoint)
							yTemp.append(self.ymin)
							zTemp.append(Zpoint)
							xyNewTemp.append(np.array([Xpoint, self.ymax]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(self.ymax)
							zNewTemp.append(Zpoint)			
							yBoundary = 2
						elif np.abs(self.ydimPos[i][j-1] - self.ymax) < 0.5*self.BoxSize:
							Xpoint, Zpoint = New_points(Point1, Point2, 'y', 1)
							xyTemp.append(np.array([Xpoint, self.ymax]))
							xTemp.append(Xpoint)
							yTemp.append(self.ymax)
							zTemp.append(Zpoint)
							xyNewTemp.append(np.array([Xpoint, self.ymin]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(self.ymin)
							zNewTemp.append(Zpoint)			
							yBoundary = 2
					elif zBoundary == 1:
						if np.abs(self.zdimPos[i][j-1] - self.zmin) < 0.5*self.BoxSize:
							Xpoint, Ypoint = New_points(Point1, Point2, 'z', 0)
							xyTemp.append(np.array([Xpoint, Ypoint]))
							xTemp.append(Xpoint)
							yTemp.append(Ypoint)
							zTemp.append(self.zmin)
							xyNewTemp.append(np.array([Xpoint, Ypoint]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(Ypoint)
							zNewTemp.append(self.zmax)			
							zBoundary = 2
						elif np.abs(self.zdimPos[i][j-1] - self.zmax) < 0.5*self.BoxSize:
							Xpoint, Ypoint = New_points(Point1, Point2, 'z', 1)
							xyTemp.append(np.array([Xpoint, Ypoint]))
							xTemp.append(Xpoint)
							yTemp.append(Ypoint)
							zTemp.append(self.zmax)
							xyNewTemp.append(np.array([Xpoint, Ypoint]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(Ypoint)
							zNewTemp.append(self.zmin)			
							zBoundary = 2
					xyNewTemp.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp.append(self.xdimPos[i][j])
					yNewTemp.append(self.ydimPos[i][j])
					zNewTemp.append(self.zdimPos[i][j])

				elif SplitFilament == 2:
					if xBoundary == 1 or yBoundary == 1 or zBoundary == 1:
						Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
						Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
					if xBoundary == 1:
						if np.abs(self.xdimPos[i][j-1] - self.xmin) < 0.5*self.BoxSize:
							Ypoint, Zpoint = New_points(Point1, Point2, 'x', 0)
							xyNewTemp.append(np.array([self.xmin, Ypoint]))
							xNewTemp.append(self.xmin)
							yNewTemp.append(Ypoint)
							zNewTemp.append(Zpoint)
							xyNewTemp2.append(np.array([self.xmax, Ypoint]))
							xNewTemp2.append(self.xmax)
							yNewTemp2.append(Ypoint)
							zNewTemp2.append(Zpoint)			
							xBoundary = 2
						elif np.abs(self.xdimPos[i][j-1] - self.xmax) < 0.5*self.BoxSize:
							Ypoint, Zpoint = New_points(Point1, Point2, 'x', 1)
							xyNewTemp.append(np.array([self.xmax, Ypoint]))
							xNewTemp.append(self.xmax)
							yNewTemp.append(Ypoint)
							zNewTemp.append(Zpoint)
							xyNewTemp2.append(np.array([self.xmin, Ypoint]))
							xNewTemp2.append(self.xmin)
							yNewTemp2.append(Ypoint)
							zNewTemp2.append(Zpoint)			
							xBoundary = 2
					elif yBoundary == 1:
						if np.abs(self.ydimPos[i][j-1] - self.ymin) < 0.5*self.BoxSize:
							Xpoint, Zpoint = New_points(Point1, Point2, 'y', 0)
							xyNewTemp.append(np.array([Xpoint, self.ymin]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(self.ymin)
							zNewTemp.append(Zpoint)
							xyNewTemp2.append(np.array([Xpoint, self.ymax]))
							xNewTemp2.append(Xpoint)
							yNewTemp2.append(self.ymax)
							zNewTemp2.append(Zpoint)			
							yBoundary = 2
						elif np.abs(self.ydimPos[i][j-1] - self.ymax) < 0.5*self.BoxSize:
							Xpoint, Zpoint = New_points(Point1, Point2, 'y', 1)
							xyNewTemp.append(np.array([Xpoint, self.ymax]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(self.ymax)
							zNewTemp.append(Zpoint)
							xyNewTemp2.append(np.array([Xpoint, self.ymin]))
							xNewTemp2.append(Xpoint)
							yNewTemp2.append(self.ymin)
							zNewTemp2.append(Zpoint)			
							yBoundary = 2
					elif zBoundary == 1:
						if np.abs(self.zdimPos[i][j-1] - self.zmin) < 0.5*self.BoxSize:
							Xpoint, Ypoint = New_points(Point1, Point2, 'z', 0)
							xyNewTemp.append(np.array([Xpoint, Ypoint]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(Ypoint)
							zNewTemp.append(self.zmin)
							xyNewTemp2.append(np.array([Xpoint, Ypoint]))
							xNewTemp2.append(Xpoint)
							yNewTemp2.append(Ypoint)
							zNewTemp2.append(self.zmax)			
							zBoundary = 2
						elif np.abs(self.zdimPos[i][j-1] - self.zmax) < 0.5*self.BoxSize:
							Xpoint, Ypoint = New_points(Point1, Point2, 'z', 1)
							xyNewTemp.append(np.array([Xpoint, Ypoint]))
							xNewTemp.append(Xpoint)
							yNewTemp.append(Ypoint)
							zNewTemp.append(self.zmax)
							xyNewTemp2.append(np.array([Xpoint, Ypoint]))
							xNewTemp2.append(Xpoint)
							yNewTemp2.append(Ypoint)
							zNewTemp2.append(self.zmin)			
							zBoundary = 2
					xyNewTemp2.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp2.append(self.xdimPos[i][j])
					yNewTemp2.append(self.ydimPos[i][j])
					zNewTemp2.append(self.zdimPos[i][j])

				# Check if boundary has been crossed
				if xDiff > 0.5*self.BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					xBoundary = 1
				if yDiff > 0.5*self.BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					yBoundary = 1
				if zDiff > 0.5*self.BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					zBoundary = 1

			# Adds final points to the positions
			if SplitFilament == 0:
				xyTemp.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xTemp.append(self.xdimPos[i][-1])
				yTemp.append(self.ydimPos[i][-1])
				zTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 1:
				xyNewTemp.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xNewTemp.append(self.xdimPos[i][-1])
				yNewTemp.append(self.ydimPos[i][-1])
				zNewTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 2:
				xyNewTemp2.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xNewTemp2.append(self.xdimPos[i][-1])
				yNewTemp2.append(self.ydimPos[i][-1])
				zNewTemp2.append(self.zdimPos[i][-1])

			# Adds positions to temporary arrays
			if len(zTemp) > 1:
				FilPosTemp.append(xyTemp)
				xPosTemp.append(xTemp)
				yPosTemp.append(yTemp)
				zPosTemp.append(zTemp)
			if SplitFilament == 1:
				if len(zNewTemp) > 1:
					FilPosTemp.append(xyNewTemp)
					xPosTemp.append(xNewTemp)
					yPosTemp.append(yNewTemp)
					zPosTemp.append(zNewTemp)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					FilPosTemp.append(xyNewTemp2)
					xPosTemp.append(xNewTemp2)
					yPosTemp.append(yNewTemp2)
					zPosTemp.append(zNewTemp2)

			TempLength = 0
			if len(zTemp) > 1:
				for k in range(len(zTemp)-1):
					TempLength += np.sqrt((xyTemp[k+1][0]-xyTemp[k][0])**2 + (xyTemp[k+1][1] - xyTemp[k][1])**2 + (zTemp[k+1]-zTemp[k])**2)
			if len(zNewTemp) > 1:
				for k in range(len(zNewTemp)-1):
					TempLength += np.sqrt((xyNewTemp[k+1][0]-xyNewTemp[k][0])**2 + (xyNewTemp[k+1][1] - xyNewTemp[k][1])**2 \
										+ (zNewTemp[k+1]-zNewTemp[k])**2)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					for k in range(len(zNewTemp2)-1):
						TempLength += np.sqrt((xyNewTemp2[k+1][0]-xyNewTemp2[k][0])**2 + (xyNewTemp2[k+1][1] - xyNewTemp2[k][1])**2 \
											+ (zNewTemp2[k+1]-zNewTemp2[k])**2)

			self.FilLengths.append(TempLength)
			if len(zTemp) > 1:
				self.LengthSplitFilament.append(TempLength)
			if len(zNewTemp) > 1: 
				self.LengthSplitFilament.append(TempLength)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					self.LengthSplitFilament.append(TempLength)
				
		
		self.FilamentPos = np.asarray(FilPosTemp)
		self.xdimPos = np.asarray(xPosTemp)
		self.ydimPos = np.asarray(yPosTemp)
		self.zdimPos = np.asarray(zPosTemp)
		self.LengthSplitFilament = np.asarray(self.LengthSplitFilament)
		self.FilLengths = np.asarray(self.FilLengths)

		print 'Boundary check time:', time.clock() - time_start, 's'

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
		#DM_points = np.dstack((self.PartPosX.ravel(), self.PartPosY.ravel(), self.PartPosZ.ravel()))
		#DM_tree = spatial.KDTree(DM_points[0])
		print np.array(self.MaskedFilamentSegments[0])[:,0].ravel()
		print self.zdimMasked[0]
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

	def Filament_Length(self, dimensions):
		""" Computes the length of the filament """
		print 'Computing filament lengths'
		if dimensions == 3:
			self.FilLengths = []
			for i in range(self.NFils-1):
				TempLength = 0
				for j in range(len(self.xdimPos[i])-1):
					TempLength += np.sqrt((self.xdimPos[i][j+1] - self.xdimPos[i][j])**2.0 + (self.ydimPos[i][j+1] - self.ydimPos[i][j])**2.0\
								 + (self.zdimPos[i][j+1] - self.zdimPos[i][j])**2.0)
				self.FilLengths.append(TempLength)

	def Plot_Figures(self, filename, ndim=3):
		""" All plots done in this function	"""
		print 'Plotting'

		if Comparison:
			NormedArgument = True
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
			
			
		if ndim == 2:
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
				plt.hold("on")
				plt.plot(self.CritPointXpos, self.CritPointYpos, 'ro', alpha=0.7, markersize=3)
				#plt.plot(self.CritPointXposNOTCON, self.CritPointYposNOTCON, 'go', alpha=0.7, markersize=3)
				plt.xlabel('$\mathregular{x}$')
				plt.ylabel('$\mathregular{y}$')
				plt.title('Position of the filaments with critical points shown')
		if ndim == 3:
			if PlotFilaments:
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
				FilPositions_WCritPts = plt.figure()
				ax = FilPositions_WCritPts.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				plt.hold("on")
				ax.plot(self.CritPointXpos, self.CritPointYpos, self.CritPointZpos, 'ro', alpha=0.7, markersize=3)
				ax.set_xlabel('$\mathregular{x}$' + LegendText)
				ax.set_ylabel('$\mathregular{y}$' + LegendText)
				ax.set_zlabel('$\mathregular{z}$' + LegendText)
				plt.title('3D Position of the filaments with critical points')
			if Projection2D:
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

				if IncludeDMParticles:
					DMParticleHist = plt.figure()
					plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
					plt.xlabel('$\mathregular{x}$' + LegendText)
					plt.ylabel('$\mathregular{y}$' + LegendText)
					plt.title('Dark matter density field over a segment of the particle box.')

					DMParticleHistwFilaments = plt.figure()
					plt.hold("on")
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

					DMParticleHistwFilamentsLengthCbar = plt.figure()
					plt.hold("on")
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.xmax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsDMlen = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMapLengthCutOff, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsDMlen)
					DMParticleHistwFilamentsLengthCbar.colorbar(line_segmentsDMlen)
					plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
					plt.xlabel('$\mathregular{x}$' + LegendText)
					plt.ylabel('$\mathregular{y}$' + LegendText)
					plt.hold("off")
					plt.title('Dark matter density field over a segment of the particle box. \n Includes filaments with colorbar.'\
							 +'Colors indicate length of a filament.')


		if self.savefile == 1:
			if HistogramPlots:
				ConnectedHist.savefig(self.results_dir + 'NumberFilamentConnectedHistogram' + self.ModelFilename)
				FilamentLengthsHist.savefig(self.results_dir + 'FilamentLengthsHistogram' + self.ModelFilename)
				FilamentPtsHis.savefig(self.results_dir + 'FilamentPointsHistogram' + self.ModelFilename)
				if IncludeDMParticles:
					NumParticlesFilamentHist.savefig(self.results_dir + 'NumParticlesPerFilament' + self.ModelFilename)
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
				if IncludeDMParticles:
					if MaskXdir == 0 and MaskYdir == 0 and MaskZdir == 1:
						DMParticleHist.savefig(self.results_dir + 'DMParticleHistogram_ZMasked' + self.ModelFilename)
						DMParticleHistwFilaments.savefig(self.results_dir + 'DMParticleHistogramWFIlaments_ZMasked' + self.ModelFilename)
						DMParticleHistwFilamentsLengthCbar.savefig(self.results_dir + 'DMParticleHistogramWFilaments_LengthCbar_ZMasked' + self.ModelFilename)
					if MaskXdir == 1 and MaskYdir == 1 and MaskZdir == 1:
						DMParticleHist.savefig(self.results_dir + 'DMParticleHistogram_XYZMasked' + self.ModelFilename)
						DMParticleHistwFilaments.savefig(self.results_dir + 'DMParticleHistogramWFIlaments_XYZMasked' + self.ModelFilename)
						DMParticleHistwFilamentsLengthCbar.savefig(self.results_dir + 'DMParticleHistogramWFilaments_LengthCbar_XYZMasked' + self.ModelFilename)

		else:
			plt.show()

		plt.close('all')

	def Solve(self, filename, ndim=3):
		self.ReadFile(filename, ndim)
		self.Sort_arrays(ndim)
		#self.BoundaryStuff()
		self.Check_boundary()
		if IncludeSlicing:
			self.Mask_slices()
		if IncludeSlicing and IncludeDMParticles:
			self.Read_SolveFile()
			#self.Mask_DMParticles()
			self.NumParticles_per_filament()
		#if IncludeDMParticles == 1:
		#	self.NumParticles_per_filament()
		
		if Comparison == 0:
			if self.savefile == 2:
				print 'Done! No files saved.'
			else:
				self.Plot_Figures(filename, ndim)
		
		return self.NumFilamentConnections, sorted(self.FilLengths), self.NFilamentPoints

class Read_solve_files():
	def __init__(self):
		self.read_solvefile()
		self.Create_Mask()
		self.Create_KDTree()

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
		#if self.ParticleComparison:
		self.Plot_Histograms_particleComparison()


	def Check_Number_Comparisons(self):
		if self.ParticleComparison:
			if self.N == 2:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles']
			elif self.N == 3:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', '$\mathregular{256^3}$ particles']
			elif self.N == 4:
				self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', '$\mathregular{256^3}$ particles', '$\mathregular{256^3}$ particles']
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
		alphas = [0.3, 0.4, 0.5, 0.6]
		ConnectedHistComparison = plt.figure()
		plt.hold("on")
		for i in range(self.N):
			DataMin = min(self.NumberConnections[i])
			DataMax = max(self.NumberConnections[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(self.NumberConnections[i], align='mid', rwidth=1, bins=BinList, normed=True, alpha=alphas[i])
		plt.xlabel('Number of connected filaments')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)
		plt.hold("off")
	
		LengthHistComparison = plt.figure()
		plt.hold("on")
		for i in range(self.N):
			plt.hist(self.FilamentLengths[i], align='mid', rwidth=1, bins=600, normed=True, histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)
		plt.hold("off")

		NPointsHistComparison = plt.figure()
		plt.hold("on")
		for i in range(self.N):
			DataMin = min(self.NPointsPerFilament[i])
			DataMax = max(self.NPointsPerFilament[i])
			BinSize = (DataMax - DataMin)/(0.5) + 1
			BinList = np.linspace(DataMin, DataMax, BinSize)
			plt.hist(self.NPointsPerFilament[i], align='mid', rwidth=1, bins=BinList, normed=True, alpha=alphas[i])
		plt.xlabel('Number of points per filament')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)
		plt.hold("off")

		if self.savefile == 1:
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison' + self.ModelFilename)
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison' + self.ModelFilename)
			NPointsHistComparison.savefig(self.results_dir + 'HistNPointsComparison' + self.ModelFilename)
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
	IncludeSlicing = 1 			# Set 1 to include slices of the box
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
	SaveAsPNG = 0 				# Set 1 to save figures as PNG
	SaveAsPDF = 1 				# Set 1 to save figures as PDF

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
		if MaskXdir == 1:
			LowerBoundaryXDir = 1*UnitConverter
			UpperBoundaryXDir = 1*UnitConverter
			print 'Masking X direction'
		if MaskYdir == 1:
			LowerBoundaryYDir = 1*UnitConverter
			UpperBoundaryYDir = 1*UnitConverter
			print 'Masking Y direction'
		if MaskZdir == 1:
			LowerBoundaryZDir = 0.47*UnitConverter
			UpperBoundaryZDir = 0.53*UnitConverter
			print 'Masking Z direction'

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
			LCDM_z0_64Instance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_64_dir+'Plots/', nPart=64, model='LCDM', redshift=0)
			NConn_64PartLCDM, FilLen_64PartLCDM, NPts_64PartLCDM = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl')
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
					NumConnections_list = [NConnections_64Peri, NConnections_128Peri]
					FilLengths_list = [FilLengths_64Peri, FilLengths_128Peri]
					FilPoints_list = [NPoints_64Peri, NPoints_128Peri]
					Histogram_instance = Histogram_Comparison(savefile=0, savefigDirectory=Comparison_dir, redshift=0, LCDM=1)
					Histogram_instance.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=64)

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
			"""
			if IncludeDMParticles == 1:
				SolveReadInstance = Read_solve_files()
				#PartPosX = SolveReadInstance.ParticlePos[:,0]
				#PartPosY = SolveReadInstance.ParticlePos[:,1]
				#PartPosZ = SolveReadInstance.ParticlePos[:,2]
				#ParticlePos = SolveReadInstance.ParticlePos
				PartPosX = SolveReadInstance.PartPosX
				PartPosY = SolveReadInstance.PartPosY
				PartPosZ = SolveReadInstance.PartPosZ
			"""
			
			LCDM_z0_64_dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/'
			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plots/', nPart=64, model='LCDM', redshift=0)
			NumConn_64LCDM, FilLen_64LCDM, NPts_64LCDM = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDMz064.a.NDskl')
			"""
			LCDM_z0_128_dir = 'lcdm_testing/LCDM_z0_128PeriodicTesting/'
			LCDM_z0_128Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128, model='LCDM', redshift=0)
			NumConn_128LCDM, FilLen_128LCDM, NPts_128LCDM = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128.a.NDskl')
			
			LCDM_z0_256_dir = 'lcdm_testing/LCDM_z0_256PeriodicTesting/'
			LCDM_z0_256Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plots/', nPart=256, model='LCDM', redshift=0)
			NumConn_256LCDM, FilLen_256LCDM, NPts_256LCDM = LCDM_z0_256Instance.Solve(LCDM_z0_256_dir+'SkelconvOutput_LCDMz0256.a.NDskl')
			
			LCDM_z0_512_dir = 'lcdm_testing/LCDM_z0_512PeriodicTesting/'
			LCDM_z0_512Instance = Disperse_Plotter(savefile=2, savefigDirectory=LCDM_z0_512_dir+'Plots/', nPart=512, model='LCDM', redshift=0)
			NumConn_512LCDM, FilLen_512LCDM, NPts_512LCDM = LCDM_z0_512Instance.Solve(LCDM_z0_512_dir+'SkelconvOutput_LCDMz0512.a.NDskl')
			"""
			if SigmaComparison:
				LCDM512_instance_nsig4 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Sigma4/', nPart=512, model='LCDM', redshift=0, SigmaArg=4)
				LCDM512_instance_nsig5 = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Sigma5/', nPart=512, model='LCDM', redshift=0, SigmaArg=5)
				NConn_512nsig4, FilLen_512nsig4, NPts_512nsig4 =LCDM512_instance_nsig4.Solve(LCDM_z0_512_dir+'Sigma4/SkelconvOutput_LCDMz0512_nsig4.a.NDskl')
				NConn_512nsig5, FilLen_512nsig5, NPts_512nsig5 =LCDM512_instance_nsig5.Solve(LCDM_z0_512_dir+'Sigma5/SkelconvOutput_LCDMz0512_nsig5.a.NDskl')
				
			Comparison_dir = 'lcdm_testing/Comparison_plots'
			if Comparison == 1 and ModelCompare == 0:
				if SigmaComparison:
					NumConnections_list = [NumConn_512LCDM, NConn_512nsig4, NConn_512nsig5]
					FilLengths_list = [FilLen_512LCDM, FilLen_512nsig4, FilLen_512nsig5]
					FilPoints_list = [NPts_512LCDM, NPts_512nsig4, NPts_512nsig5]
					ComparisonInstance_LCDM = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir+'SigmaComparisons/',\
										 redshift=0, LCDM=1, nsigComparison=1)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=512)
				else:
					NumConnections_list = [NumConn_64LCDM, NumConn_128LCDM] # , NumConn_256LCDM, NumConn_512LCDM]
					FilLengths_list = [FilLen_64LCDM, FilLen_128LCDM] #, FilLen_256LCDM, FilLen_512LCDM]
					FilPoints_list = [NPts_64LCDM, NPts_128LCDM] #, NPts_256LCDM, NPts_512LCDM]
					ComparisonInstance_LCDM = Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, redshift=0, LCDM=1)
					ComparisonInstance_LCDM.Run(NumConnections_list, FilLengths_list, FilPoints_list, nPart=512)
					
				"""
				ComparisonInstance = Histogram_Comparison2(savefile=0, savefigDirectory=Comparison_dir, ndim=3, model='$\mathregular{\Lambda}$CDM',\
													 Nparticles=[64, 128])
				ComparisonInstance.Solve('Connections', NumConnections_list)
				ComparisonInstance.Solve('Lengths', FilLengths_list)
				ComparisonInstance.Solve('FilamentPoints', FilPoints_list)
				plt.show()
				"""