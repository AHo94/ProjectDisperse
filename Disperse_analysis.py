### Set comment on the two below to plot. Only use if running on papsukal, nekkar etc. 
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import os
import scipy.stats as stats
from matplotlib import colors as mcolors
from scipy import interpolate

class Disperse_Plotter():
	def __init__(self, savefile, savefigDirectory, nPart):
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

	def ReadFile(self, filename, dimensions):
		""" Reads the data from the skeleton file from Disperse """
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
				self.xmin = float(MinValues[0])
				self.xmax = float(MaxValues[0])
				self.ymin = float(MinValues[1])
				self.ymax = float(MaxValues[1])
				if dimensions == 3:
					self.zmin = float(MinValues[2])
					self.zmax = float(MaxValues[2])
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

	def Read_SolveFile(self):
		print 'Reading data for the file: ', solve_file_dir, solve_filename, '. May take a while...' 
		datafiles = open(os.path.join(file_directory+solve_file_dir, solve_filename), 'r')
		self.PartPosX = []
		self.PartPosY = []
		PartPosZ = []
		
		PartVelX = []
		PartVelY = []
		PartVelZ = []
		self.XYPartPos = []
		SkipFirstLine = 0
		for line in datafiles:
			data_set = line.split()
			if IncludeSlicing == 1:
				if SkipFirstLine == 0:
					SkipFirstLine = 1
				else:
					if float(data_set[2]) > LowerBoundary and float(data_set[2]) < UpperBoundary:
						#self.XYPartPos.append([float(data_set[0]), float(data_set[1])])
						self.PartPosX.append(float(data_set[0]))
						self.PartPosY.append(float(data_set[1]))
						#PartPosZ.append(float(data_set[2]))
						#PartVelX.append(float(data_set[3]))
						#PartVelY.append(float(data_set[4]))
						#PartVelZ.append(float(data_set[5]))

	def Sort_arrays(self, dimensions):
		""" 
		Sorts data to their respective arrays 
		Data to be sorted: Critical points, ID of filament and filament points
		"""
		print 'Sorting data ...'
		self.NcritPts = int(self.CriticalPoints[0][0])
		if FilamentLimit == 0:
			self.NFils = int(self.Filaments[0][0])
		else:
			self.NFils = FilamentLimit

		# General data
		self.NumFilamentConnections = []
		self.IDFilamentConnected = []
		self.CritPointInfo = []
		for i in range(1, len(self.CriticalPoints)-1):
			stuff = self.CriticalPoints[i]
			if len(stuff) == 1:
				self.NumFilamentConnections.append(int(stuff[0]))
				IDconnections = []
				if int(stuff[0]) == 1:
					IDS = self.CriticalPoints[i+1]
					IDconnections.append([float(IDS[0]), float(IDS[1])])
				else:	
					for j in range(1, int(stuff[0])+1):
						IDS = self.CriticalPoints[i+j]
						IDconnections.append([float(IDS[0]),float(IDS[1])])
				self.IDFilamentConnected.append(IDconnections)
			elif len(stuff) > 2:
				InfoListTemp = []
				for k in range(len(stuff)):
					InfoListTemp.append(float(stuff[k]))
				self.CritPointInfo.append(InfoListTemp)
		
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
					xPos = float(self.Filaments[k+j][0])
					yPos = float(self.Filaments[k+j][1])
					TempPositions.append([xPos, yPos])
					xtemp.append(xPos)
					ytemp.append(yPos)
				self.FilamentPos.append(TempPositions)
				self.xdimPos.append(xtemp)
				self.ydimPos.append(ytemp)
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
					xPos = float(self.Filaments[k+j][0])
					yPos = float(self.Filaments[k+j][1])
					zPos = float(self.Filaments[k+j][2])
					TempPositions.append([xPos, yPos])
					xtemp.append(xPos)
					ytemp.append(yPos)
					ztemp.append(zPos)		
				self.NFilamentPoints.append(int(Filstuff[-1]))
				self.FilamentPos.append(TempPositions)
				self.xdimPos.append(xtemp)
				self.ydimPos.append(ytemp)
				self.zdimPos.append(ztemp)
				k += int(Filstuff[-1])+1

				if FilamentLimit != 0:
					if k >= FilamentLimit:
						self.NFils = len(self.xdimPos)
						break

	def Mask_slices(self):
		"""
		Creates a mask to plot a slice of the filament box. Boundary of slice chosen arbitrarily.
		The masking includes filaments that are within the given boundary.
		Also includes masking where the filament are cut off outside the boundary. 
		"""
		print 'Computing masks'
		self.zdimMasked = []
		self.MaskedFilamentSegments = []
		for i in range(len(self.zdimPos)):
			if np.any(np.array(self.zdimPos[i]) > LowerBoundary) and np.any(np.array(self.zdimPos[i]) < UpperBoundary):
				self.zdimMasked.append(self.zdimPos[i])
				self.MaskedFilamentSegments.append(self.FilamentPos[i])

		self.CutOffFilamentSegments = []
		self.CutOffzDim = []

		for i in range(len(self.zdimPos)):
			Indices = np.where(np.logical_and(np.greater(self.zdimPos[i], LowerBoundary), np.less(self.zdimPos[i], UpperBoundary)))[0]
			if not len(Indices) == 0:
				FilSegmentTemp = []
				zDimTemp = []
				for idx in Indices:
					zDimTemp.append(self.zdimPos[i][idx])
					FilSegmentTemp.append([self.xdimPos[i][idx], self.ydimPos[i][idx]])
				self.CutOffFilamentSegments.append(FilSegmentTemp)
				self.CutOffzDim.append(zDimTemp)


	def Check_Boundary_and_compute_length(self):
		print 'Checking boundaries'
		BoxSize = self.xmax - self.xmin
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
		for i in range(len(self.xdimPos)):
			SplitFilament = 0
			xyTemp = []
			xTemp = []
			yTemp = []
			zTemp = []
			xyNewTemp = []
			xNewTemp = []
			yNewTemp = []
			zNewTemp = []
			for j in range(len(self.xdimPos[i])-1):
				xDiff = np.abs(self.xdimPos[i][j+1] - self.xdimPos[i][j])
				yDiff = np.abs(self.ydimPos[i][j+1] - self.ydimPos[i][j])
				zDiff = np.abs(self.zdimPos[i][j+1] - self.zdimPos[i][j])

				if SplitFilament == 0:
					xyTemp.append([self.xdimPos[i][j], self.ydimPos[i][j]])
					xTemp.append(self.xdimPos[i][j])
					yTemp.append(self.ydimPos[i][j])
					zTemp.append(self.zdimPos[i][j])
				elif SplitFilament == 1:
					xyNewTemp.append([self.xdimPos[i][j], self.ydimPos[i][j]])
					xNewTemp.append(self.xdimPos[i][j])
					yNewTemp.append(self.ydimPos[i][j])
					zNewTemp.append(self.zdimPos[i][j])
				if xDiff > 0.5*BoxSize or yDiff > 0.5*BoxSize or zDiff > 0.5*BoxSize:
					SplitFilament += 1			

			if SplitFilament == 0:
				xyTemp.append([self.xdimPos[i][-1], self.ydimPos[i][-1]])
				xTemp.append(self.xdimPos[i][-1])
				yTemp.append(self.ydimPos[i][-1])
				zTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 1:
				xyNewTemp.append([self.xdimPos[i][-1], self.ydimPos[i][-1]])
				xNewTemp.append(self.xdimPos[i][-1])
				yNewTemp.append(self.ydimPos[i][-1])
				zNewTemp.append(self.zdimPos[i][-1])

			if len(zTemp) > 1:
				FilPosTemp.append(xyTemp)
				xPosTemp.append(xTemp)
				yPosTemp.append(yTemp)
				zPosTemp.append(zTemp)
			if len(zNewTemp) > 1:
				FilPosTemp.append(xyNewTemp)
				xPosTemp.append(xNewTemp)
				yPosTemp.append(yNewTemp)
				zPosTemp.append(zNewTemp)

			TempLength = 0
			if len(zTemp) > 1:
				for k in range(len(zTemp)-1):
						TempLength += np.sqrt((xyTemp[k+1][0]-xyTemp[k][0])**2 + (xyTemp[k+1][1] - xyTemp[k][1])**2 + (zTemp[k+1]-zTemp[k])**2)
			if len(zNewTemp) > 1:
				for k in range(len(zNewTemp)-1):
						TempLength += np.sqrt((xyNewTemp[k+1][0]-xyNewTemp[k][0])**2 + (xyNewTemp[k+1][1] - xyNewTemp[k][1])**2 \
											+ (zNewTemp[k+1]-zNewTemp[k])**2)

			self.FilLengths.append(TempLength)
		
		self.FilamentPos = FilPosTemp
		self.xdimPos = xPosTemp
		self.ydimPos = yPosTemp
		self.zdimPos = zPosTemp

	def BoundaryStuff(self):
		"""
		This function checks whether a filament crosses the boundary or not. 
		If a filament crosses the boundary, the filament will be split into two/three filaments, i.e different lists/arrays.
		The algorithm of splitting the filament:
		1) Check which direction the filament crosses the boundary.
		2) For the point we are considering, check which boundary is closest to the point and add that to the old list. 
		3) Create new list and add other side of the boundary as the first point in the new list.
		4) For the other directions, add 'directionPoint*2/BoxSize' to the old list and subtract the same to the new list.
		5) Repeat if filament crosses boundary multiple times.
		Function also computes the length of the filament.
		"""
		BoxSize = self.xmax - self.xmin
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
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
					xyTemp.append([self.xdimPos[i][j], self.ydimPos[i][j]])
					xTemp.append(self.xdimPos[i][j])
					yTemp.append(self.ydimPos[i][j])
					zTemp.append(self.zdimPos[i][j])
				elif SplitFilament == 1:
					if xBoundary == 1:
						yBoundaryPoint = 2.0*(self.ydimPos[i][j] - self.ydimPos[i][j-1])/BoxSize
						zBoundaryPoint = 2.0*(self.zdimPos[i][j] - self.zdimPos[i][j-1])/BoxSize
						if np.abs(self.xdimPos[i][j-1] - self.xmin) > 0.5*BoxSize:
							xyTemp.append([self.xmax, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp.append([self.xmin, self.ydimPos[i][j-1]-yBoundaryPoint])
							xTemp.append(self.xmax)
							xNewTemp.append(self.xmin)
							yTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							xBoundary = 2
						elif np.abs(self.xdimPos[i][j-1] - self.xmax) > 0.5*BoxSize: 
							xyTemp.append([self.xmin, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp.append([self.xmax, self.ydimPos[i][j-1]-yBoundaryPoint])
							xTemp.append(self.xmin)
							xNewTemp.append(self.xmax)
							yTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							xBoundary = 2
					elif yBoundary == 1:
						xBoundaryPoint = 2.0*(self.xdimPos[i][j] - self.xdimPos[i][j-1])/BoxSize
						zBoundaryPoint = 2.0*(self.zdimPos[i][j] - self.zdimPos[i][j-1])/BoxSize
						if np.abs(self.ydimPos[i][j-1] - self.ymin) > 0.5*BoxSize:
							xyTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ymax])
							xyNewTemp.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ymin])
							xTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yTemp.append(self.ymax)
							yNewTemp.append(self.ymin)
							zTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							yBoundary = 2
						elif np.abs(self.ydimPos[i][j-1] - self.ymax) > 0.5*BoxSize: 
							xyTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ymin])
							xyNewTemp.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ymax])
							xTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yTemp.append(self.ymin)
							yNewTemp.append(self.ymax)
							zTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							yBoundary = 2
					elif zBoundary == 1:
						xBoundaryPoint = 2.0*(self.xdimPos[i][j] - self.xdimPos[i][j-1])/BoxSize
						yBoundaryPoint = 2.0*(self.ydimPos[i][j] - self.ydimPos[i][j-1])/BoxSize
						if np.abs(self.zdimPos[i][j-1] - self.zmin) > 0.5*BoxSize:
							xyTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ydimPos[i][j-1]-yBoundaryPoint])
							xTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zTemp.append(self.zmax)
							zNewTemp.append(self.zmin)
							zBoundary = 2
						elif np.abs(self.zdimPos[i][j-1] - self.zmax) > 0.5*BoxSize:
							xyTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ydimPos[i][j-1]-yBoundaryPoint])
							xTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zTemp.append(self.zmin)
							zNewTemp.append(self.zmax)
							zBoundary = 2
					xyNewTemp.append([self.xdimPos[i][j], self.ydimPos[i][j]])
					xNewTemp.append(self.xdimPos[i][j])
					yNewTemp.append(self.ydimPos[i][j])
					zNewTemp.append(self.zdimPos[i][j])
				elif SplitFilament == 2:
					if xBoundary == 1:
						yBoundaryPoint = 2.0*(self.ydimPos[i][j] - self.ydimPos[i][j-1])/BoxSize
						zBoundaryPoint = 2.0*(self.zdimPos[i][j] - self.zdimPos[i][j-1])/BoxSize
						if np.abs(self.xdimPos[i][j-1] - self.xmin) > 0.5*BoxSize:
							xyNewTemp.append([self.xmax, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp2.append([self.xmin, self.ydimPos[i][j-1]-yBoundaryPoint])
							xNewTemp.append(self.xmax)
							xNewTemp2.append(self.xmin)
							yNewTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp2.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp2.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							xBoundary = 2
						elif np.abs(self.xdimPos[i][j-1] - self.xmax) > 0.5*BoxSize: 
							xyNewTemp.append([self.xmin, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp2.append([self.xmax, self.ydimPos[i][j-1]-yBoundaryPoint])
							xNewTemp.append(self.xmin)
							xNewTemp2.append(self.xmax)
							yNewTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp2.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zNewTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp2.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							xBoundary = 2
					elif yBoundary == 1:
						xBoundaryPoint = 2.0*(self.xdimPos[i][j] - self.xdimPos[i][j-1])/BoxSize
						zBoundaryPoint = 2.0*(self.zdimPos[i][j] - self.zdimPos[i][j-1])/BoxSize
						if np.abs(self.ydimPos[i][j-1] - self.ymin) > 0.5*BoxSize:
							xyNewTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ymax])
							xyNewTemp2.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ymin])
							xNewTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp2.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yNewTemp.append(self.ymax)
							yNewTemp2.append(self.ymin)
							zNewTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp2.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							yBoundary = 2
						elif np.abs(self.ydimPos[i][j-1] - self.ymax) > 0.5*BoxSize: 
							xyNewTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ymin])
							xyNewTemp2.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ymax])
							xNewTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp2.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yNewTemp.append(self.ymin)
							yNewTemp2.append(self.ymax)
							zNewTemp.append(self.zdimPos[i][j-1]+zBoundaryPoint)
							zNewTemp2.append(self.zdimPos[i][j-1]-zBoundaryPoint)
							yBoundary = 2
					elif zBoundary == 1:
						xBoundaryPoint = 2.0*(self.xdimPos[i][j] - self.xdimPos[i][j-1])/BoxSize
						yBoundaryPoint = 2.0*(self.ydimPos[i][j] - self.ydimPos[i][j-1])/BoxSize
						if np.abs(self.zdimPos[i][j-1] - self.zmin) > 0.5*BoxSize:
							xyNewTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp2.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ydimPos[i][j-1]-yBoundaryPoint])
							xNewTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp2.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp2.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zNewTemp.append(self.zmax)
							zNewTemp2.append(self.zmin)
							zBoundary = 2
						elif np.abs(self.zdimPos[i][j-1] - self.zmax) > 0.5*BoxSize:
							xyNewTemp.append([self.xdimPos[i][j-1]+xBoundaryPoint, self.ydimPos[i][j-1]+yBoundaryPoint])
							xyNewTemp2.append([self.xdimPos[i][j-1]-xBoundaryPoint, self.ydimPos[i][j-1]-yBoundaryPoint])
							xNewTemp.append(self.xdimPos[i][j-1]+xBoundaryPoint)
							xNewTemp2.append(self.xdimPos[i][j-1]-xBoundaryPoint)
							yNewTemp.append(self.ydimPos[i][j-1]+yBoundaryPoint)
							yNewTemp2.append(self.ydimPos[i][j-1]-yBoundaryPoint)
							zNewTemp.append(self.zmin)
							zNewTemp2.append(self.zmax)
							zBoundary = 2
					xyNewTemp2.append([self.xdimPos[i][j], self.ydimPos[i][j]])
					xNewTemp2.append(self.xdimPos[i][j])
					yNewTemp2.append(self.ydimPos[i][j])
					zNewTemp2.append(self.zdimPos[i][j])

				if xDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					xBoundary = 1
				if yDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					yBoundary = 1
				if zDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					zBoundary = 1			

			if SplitFilament == 0:
				xyTemp.append([self.xdimPos[i][-1], self.ydimPos[i][-1]])
				xTemp.append(self.xdimPos[i][-1])
				yTemp.append(self.ydimPos[i][-1])
				zTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 1:
				xyNewTemp.append([self.xdimPos[i][-1], self.ydimPos[i][-1]])
				xNewTemp.append(self.xdimPos[i][-1])
				yNewTemp.append(self.ydimPos[i][-1])
				zNewTemp.append(self.zdimPos[i][-1])

			if len(zTemp) > 1:
				FilPosTemp.append(xyTemp)
				xPosTemp.append(xTemp)
				yPosTemp.append(yTemp)
				zPosTemp.append(zTemp)
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
		
		self.FilamentPos = FilPosTemp
		self.xdimPos = xPosTemp
		self.ydimPos = yPosTemp
		self.zdimPos = zPosTemp


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

	def Plot_Figures(self, filename, ndim=2):
		""" All plots done in this function	"""
		if ndim == 2:
			Filenamedimension = '2D.png'
		elif ndim == 3:
			Filenamedimension = '3D.png'
		else:
			raise ValueError('ndim less than 2 or greater than 3!')
		print 'Plotting'

		if Comparison == 1:
			NormedArgument = True
		else:
			NormedArgument = False

		if FilamentColors == 1:
			ColorArray = np.linspace(0,1,110)
			ColorMap2D = np.array([np.mean(zvalues) for zvalues in self.zdimPos])
		else:
			ColorArray = None
		

		if HistogramPlots == 1:
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
						%self.NcritPts + Filenamedimension.replace('.png', '')+ '. Using '+ self.nPart_text+'$\mathregular{^3}$ particles.')
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
			plt.xlabel('Length of filaments ')
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
			
		if ndim == 2:
			if PlotFilaments == 1:
				FilPositions = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection(line_segments)
				plt.xlabel('$\mathregular{x}$')
				plt.ylabel('$\mathregular{y}$')
				plt.title('Positions of the filaments with '+ self.nPart_text+ '$^3$ particles')
			if PlotFilamentsWCritPts == 1:
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
			if PlotFilaments == 1:
				FilPositions = plt.figure()
				ax = FilPositions.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.rainbow)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				ax.set_xlabel('$\mathregular{x}$')
				ax.set_ylabel('$\mathregular{y}$')
				ax.set_zlabel('$\mathregular{z}$')
				plt.title('3D Position of the filaments with '+ self.nPart_text+ '$^3$ particles.')
			if PlotFilamentsWCritPts == 1:
				FilPositions_WCritPts = plt.figure()
				ax = FilPositions_WCritPts.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				plt.hold("on")
				ax.plot(self.CritPointXpos, self.CritPointYpos, self.CritPointZpos, 'ro', alpha=0.7, markersize=3)
				ax.set_xlabel('$\mathregular{x}$')
				ax.set_ylabel('$\mathregular{y}$')
				ax.set_zlabel('$\mathregular{z}$')
				plt.title('3D Position of the filaments with critical points')
			if Projection2D == 1:
				FilPositions_2DProjection = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, array=ColorArray, cmap=plt.cm.rainbow)
				ax.add_collection(line_segments)
				ax.set_xlabel('$\mathregular{x}$')
				ax.set_ylabel('$\mathregular{y}$')
				plt.title('2D Position projection of the filaments')
				if ColorBarZDir == 1:
					FilPositions_2DProjectionColorBarZDir = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.FilamentPos, array=ColorMap2D, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DProjectionColorBarZDir.colorbar(line_segmentsCbar)
					ax.set_xlabel('$\mathregular{x}$')
					ax.set_ylabel('$\mathregular{y}$')
					plt.title('2D Position projection of the filaments.\n Color based on average z-position')

			if IncludeSlicing == 1:
				self.Mask_slices()
				FilamentSliced = plt.figure()
				ax = FilamentSliced.gca(projection='3d')
				ax.set_xlim(self.xmin, self.xmax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.MaskedFilamentSegments, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimMasked, zdir='z')
				ax.set_xlabel('$\mathregular{x}$')
				ax.set_ylabel('$\mathregular{y}$')
				ax.set_zlabel('$\mathregular{z}$')
				plt.title('Sliced segment of the 3D box')

				FilamentCutOff = plt.figure()
				ax2 = FilamentCutOff.gca(projection='3d')
				ax2.set_xlim(self.xmin, self.xmax)
				ax2.set_ylim(self.ymin, self.ymax)
				ax2.set_zlim(self.zmin, self.zmax)
				line_segments_CO = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax2.add_collection3d(line_segments_CO, self.CutOffzDim, zdir='z')
				ax2.set_xlabel('$\mathregular{x}$')
				ax2.set_ylabel('$\mathregular{y}$')
				ax2.set_zlabel('$\mathregular{z}$')
				plt.title('Sliced segment of the 3D box, with filaments cut off outside of boundary')
				
				if Projection2D == 1:
					FilPositions_2DProjectionSliced = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segments2 = LineCollection(self.MaskedFilamentSegments, array=ColorArray, cmap=plt.cm.gist_ncar)
					ax.add_collection(line_segments2)
					ax.set_xlabel('$\mathregular{x}$')
					ax.set_ylabel('$\mathregular{y}$')
					plt.title('2D Position projection sliced box')

				if ColorBarZDir == 1:
					FilPositions_2DSlicedProjectionColorBarZDir = plt.figure()
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.ymax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsCbar = LineCollection(self.FilamentPos, array=ColorMap2D, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsCbar)
					FilPositions_2DSlicedProjectionColorBarZDir.colorbar(line_segmentsCbar)
					ax.set_xlabel('$\mathregular{x}$')
					ax.set_ylabel('$\mathregular{y}$')
					plt.title('2D projection of the filaments in a sliced segment of the box.\n Color based on average z-position')

				if IncludeDMParticles == 1:
					DMParticleHist = plt.figure()
					plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
					plt.xlabel('$\mathregular{x}$')
					plt.ylabel('$\mathregular{y}$')
					plt.title('Dark matter density field over a segment of the particle box.')

					DMParticleHistwFilaments = plt.figure()
					plt.hold("on")
					ax = plt.axes()
					ax.set_xlim(self.xmin, self.xmax)
					ax.set_ylim(self.ymin, self.ymax)
					line_segmentsDM = LineCollection(self.CutOffFilamentSegments, linestyle='solid', array=ColorMap2D, cmap=plt.cm.rainbow)
					ax.add_collection(line_segmentsDM)
					DMParticleHistwFilaments.colorbar(line_segmentsDM)
					plt.hist2d(self.PartPosX, self.PartPosY, bins=100)
					plt.xlabel('$\mathregular{x}$')
					plt.ylabel('$\mathregular{y}$')
					plt.hold("off")
					plt.title('Dark matter density field over a segment of the particle box. \n Includes filaments with colorbar. Colors indicate average z-value.')

		if self.savefile == 1:
			if HistogramPlots == 1:
				ConnectedHist.savefig(self.results_dir + 'NumberFilamentConnectedHistogram' + Filenamedimension)
				FilamentLengthsHist.savefig(self.results_dir + 'FilamentLengthsHistogram' + Filenamedimension)
				FilamentPtsHis.savefig(self.results_dir + 'FilamentPointsHistogram' + Filenamedimension)
			if PlotFilaments == 1:
				FilPositions.savefig(self.results_dir + 'FilamentPositions' + Filenamedimension)
			if PlotFilamentsWCritPts == 1:
				FilPositions_WCritPts.savefig(self.results_dir + 'FilamentPositionsWithCritPts' + Filenamedimension)
			if Projection2D == 1:
				FilPositions_2DProjection.savefig(self.results_dir + '2DProjectedFilamentPosition' + Filenamedimension)
				if ColorBarZDir == 1:
					FilPositions_2DProjectionColorBarZDir.savefig(self.results_dir + '2DProjectionColorBarZDir' + Filenamedimension)
			if IncludeSlicing == 1:
				FilamentSliced.savefig(self.results_dir + 'Sliced3dBox' + Filenamedimension)
				FilamentCutOff.savefig(self.results_dir + 'CutOffFilaments' + Filenamedimension)
				if Projection2D == 1:
					FilPositions_2DProjectionSliced.savefig(self.results_dir + '2DProjectionSliced3dBox' + Filenamedimension)
				if IncludeDMParticles == 1:
					DMParticleHist.savefig(self.results_dir + 'DMParticleHistogram' + Filenamedimension)
					DMParticleHistwFilaments.savefig(self.results_dir + 'DMParticleHistogramWFIlaments' + Filenamedimension)

		else:
			plt.show()

	def Solve(self, filename, ndim=2):
		self.ReadFile(filename, ndim)
		if IncludeDMParticles == 1:
			self.Read_SolveFile()
		self.Sort_arrays(ndim)
		#self.Check_Boundary_and_compute_length()
		self.BoundaryStuff()
		#self.Filament_Length(ndim)
		
		#if ndim == 3:
		#	self.Check_Boundaries()
		
		if Comparison == 0:
			if self.savefile == 2:
				print 'Done! No files saved.'
			else:
				self.Plot_Figures(filename, ndim)
		
		return self.NumFilamentConnections, sorted(self.FilLengths), self.NFilamentPoints

class Histogram_Comparison():
	def __init__(self, savefile, savefigDirectory, ndim, NumberConnections, FilamentLengths):
		self.savefile = savefile
		self.ndim = ndim
		
		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir):
			os.makedirs(self.results_dir)

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

		for i in range(len(NumberConnections)-1):
			if len(NumberConnections[i+1]) < len(NumberConnections[i]):
				raise ValueError('List of NumberConnections must be in order to increasing number of particles')
			elif len(FilamentLengths[i+1]) < len(FilamentLengths[i]):
				raise ValueError('List of FilamentLengths must be in order to increasing number of particles')
		

		self.NumberConnections = NumberConnections
		self.FilamentLengths = FilamentLengths
		self.N = len(self.NumberConnections)

		self.Check_Number_Comparisons()
		self.Plot_Histograms()

	def Check_Number_Comparisons(self):
		if self.N == 2:
			self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles']
		elif self.N == 3:
			self.LegendText = ['$\mathregular{64^3}$ particles', '$\mathregular{128^3}$ particles', '$\mathregular{256^3}$ particles']
			
	def Plot_Histograms(self):
		if self.ndim == 2:
			Filenamedimension = '2D.png'
		elif self.ndim == 3:
			Filenamedimension = '3D.png'
		else:
			raise ValueError('ndim less than 2 or greater than 3!')
		print 'Plotting'

		alphas = [0.4, 0.5, 0.6, 0.7]
		
		ConnectedHistComparison = plt.figure()
		plt.hold("on")
		for i in range(self.N):
			PtsMin = min(self.NumberConnections[i])
			PtsMax = max(self.NumberConnections[i])
			BinSize_FilPts = (PtsMax - PtsMin)/(0.5) + 1
			BinList_FilPts = np.linspace(PtsMin, PtsMax, BinSize_FilPts)
			plt.hist(self.NumberConnections[i], align='mid', rwidth=1, bins=BinList_FilPts, normed=True, alpha=alphas[i])
		plt.xlabel('Number of connected filaments')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)
	
		LengthHistComparison = plt.figure()
		plt.hold("on")
		for i in range(self.N):
			plt.hist(self.FilamentLengths[i], align='mid', rwidth=1, bins=600, normed=True, alpha=alphas[i], histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)

		if self.savefile == 1:
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison' + Filenamedimension)
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison' + Filenamedimension)
		elif self.savefile == 2:
			print 'Done! No histograms to plot.'
		else:
			plt.show()


class Histogram_Comparison2():
	def __init__(self, savefile, savefigDirectory, ndim, model, Nparticles):
		self.savefile = savefile
		self.ndim = ndim
		self.model = model
		self.Nparticles = Nparticles
		self.alphas = [0.4, 0.5, 0.6, 0.7]
		
		if type(Nparticles) != list:
			raise ValueError('Argument Nparticles must be a list!')
		if type(model) != str:
			raise ValueError('Argument model must be a string!')

		if self.ndim == 2:
			self.Filenamedimension = '2D.png'
		elif self.ndim == 3:
			self.Filenamedimension = '3D.png'
		else:
			raise ValueError('ndim less than 2 or greater than 3!')
		print 'Plotting histogram comparisons'

		self.results_dir = os.path.join(savefile_directory, savefigDirectory)
		if not os.path.isdir(self.results_dir):
			os.makedirs(self.results_dir)

		self.Check_Number_Comparisons()

	def Check_Number_Comparisons(self):
		self.LegendText = []
		for particles in self.Nparticles:
			if particles == 64:
				self.LegendText.append('$\mathregular{64^3}$ particles')
			elif particles == 128:
				self.LegendText.append('$\mathregular{128^3}$ particles')
			elif particles == 256:
				self.LegendText.append('$\mathregular{256^3}$ particles')
			elif particles == 512:
				self.LegendText.append('$\mathregular{512^3}$ particles')
			else:
				raise ValueError('Invalid argument: particles = ' + str(particles))

	def Plot_Filament_Connections(self, histograms):
		ConnectedHistComparison = plt.figure()
		plt.hold("on")
		for i in range(len(histograms)):
			NumMin = min(histograms[i])
			NumMax = max(histograms[i])
			BinSize = (NumMax - NumMin)/(0.5) + 1
			Bin_list = np.linspace(NumMin, NumMax, BinSize)
			plt.hist(histograms[i], align='mid', rwidth=1, bins=50, normed=True, alpha=self.alphas[i])
		plt.xlabel('Number of connected filaments')
		plt.legend(self.LegendText)
		plt.title('Comparison of number of connected points for the model: ' + self.model)
		if self.savefile == 1:
			ConnectedHistComparison.savefig(self.results_dir + 'HistNumConnectedFilamentsComparison' + self.Filenamedimension)
		
	def Plot_Filament_Lengths(self, histograms):
		LengthHistComparison = plt.figure()
		plt.hold("on")
		for i in range(len(histograms)):
			plt.hist(histograms[i], align='mid', rwidth=1, bins=600, normed=True, alpha=self.alphas[i], histtype='step')
		plt.xlabel('Filament lengths')
		plt.ylabel('Number of occurances')
		plt.legend(self.LegendText)
		plt.title('Comparison of filament lengths for the model: ' + self.model)
		if self.savefile == 1:
			LengthHistComparison.savefig(self.results_dir + 'HistLengthComparison' + self.Filenamedimension)

	def Plot_Filament_Points(self, histograms):
		FilPointsHistComparison = plt.figure()
		plt.hold("on")
		for i in range(len(histograms)):
			NumMin = min(histograms[i])
			NumMax = max(histograms[i])
			BinSize = (NumMax - NumMin)/(0.5) + 1
			Bin_list = np.linspace(NumMin, NumMax, BinSize)
			#fit = stats.norm.pdf(sorted(histograms[i]), np.mean(sorted(histograms[i])), np.std(sorted(histograms[i])))
			#plt.plot(sorted(histograms[i]), fit)
			plt.hist(histograms[i], align='mid', rwidth=1, bins=Bin_list, normed=True, alpha=self.alphas[i])
		plt.xlabel('Number of points in a filament')
		plt.legend(self.LegendText)
		plt.title('Comparison of number of points in a filament for model: ' + self.model)
		if self.savefile == 1:
			FilPointsHistComparison.savefig(self.results_dir + 'HisFilamentPointsComparison' + self.Filenamedimension)

	def Solve(self, FilamentProperty, FilamentPropertyList):
		if len(FilamentPropertyList) < 2:
			raise ValueError('Less than two histograms to compare! Stopping program.')

		if FilamentProperty == 'Connections':
			self.Plot_Filament_Connections(FilamentPropertyList)
		elif FilamentProperty == 'Lengths':
			self.Plot_Filament_Lengths(FilamentPropertyList)
		elif FilamentProperty == 'FilamentPoints':
			self.Plot_Filament_Points(FilamentPropertyList)
		else:
			print 'Arguments not properly set! Try the following (arguments must be string):'
			print "FilamentProperty = 'Connections'"
			print "FilamentProperty = 'Lengths'"
			print "FilamentProperty = 'FilamentPoints'"
			return

		if self.savefile == 2:
			print 'Done! Nothing saved.'
		#else:
		#	plt.show()

if __name__ == '__main__':
	HOMEPC = 1					# Set 1 if working in UiO terminal
	FilamentLimit = 0			# Limits the number of lines read from file. Reads all if 0
	
	PlotFilaments = 0			# Set 1 to plot actual filaments
	PlotFilamentsWCritPts = 0	# Set to 1 to plot filaments with critical points
	Projection2D = 0			# Set to 1 to plot a 2D projection of the 3D case
	ColorBarZDir = 0
	FilamentColors = 1 			# Set to 1 to get different colors for different filaments
	IncludeDMParticles = 0 		# Set to 1 to include dark matter particle plots
	IncludeSlicing = 0 			# Set 1 to include slices of the box
	if IncludeSlicing == 1:
		LowerBoundary = 0.45
		UpperBoundary = 0.55

	HistogramPlots = 1 			# Set to 1 to plot histograms
	Comparison = 0				# Set 1 if you want to compare different number of particles. Usual plots will not be plotted!

	# Run simulation for different models. Set to 1 to run them. 
	LCDM_model = 1 				
	
	if HOMEPC == 0:
		file_directory = 'C:/Users/Alex/Documents/Masters_project/Disperse'
		savefile_directory = file_directory
		IncludeDMParticles = 0
		#solveInstance1 = Disperse_Plotter(savefile=1, savefigDirectory='Plot_Disperse_Example/', nPart=64)
		#solveInstance1.Plot("simu_2D.ND.NDnet_s3.up.NDskl.a.NDskl", ndim=2)
		#solveInstance1.Plot("simu_32_id.gad.NDnet_s3.5.up.NDskl.a.NDskl", ndim=3)
		"""
		LCDM_64Periodic_dir = 'lcdm_z0_testing/LCDM64_Periodic/'
		LCDM_z0_64Peri = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_64Periodic_dir + 'Plots/', nPart=64)
		NConnections_64Peri, FilLengths_64Peri, NPoints_64Peri = LCDM_z0_64Peri.Solve(LCDM_64Periodic_dir+'SkelconvOutput_LCDM64Periodic.a.NDskl',ndim=3)
		"""
		"""
		LCDM_128Periodic_dir = 'lcdm_z0_testing/LCDM128_Periodic/'
		LCDM_z0_128Peri = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_128Periodic_dir + 'Plots/', nPart=128)
		NConnections_128Peri, FilLengths_128Peri, NPoints_128Peri = LCDM_z0_128Peri.Solve(LCDM_128Periodic_dir+'SkelconvOutput_LCDM128Periodic.a.NDskl',ndim=3)
		"""
		
		LCDM_z0_64Test2_dir = 'lcdm_z0_testing/LCDM_z0_64PeriodicTesting/'
		LCDM_z0_64Test2Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64Test2_dir+'Plots/', nPart=64)
		NN, FF, FP = LCDM_z0_64Test2Instance.Solve(LCDM_z0_64Test2_dir+'SkelconvOutput_LCDMz064.a.NDskl', ndim=3)

		Comparison_dir = 'lcdm_z0_testing/Comparison_plots/'
		if Comparison == 1:
			NumConnections_list = [NConnections_64Peri, NConnections_128Peri]
			FilLengths_list = [FilLengths_64Peri, FilLengths_128Peri]
			FilPoints_list = [NPoints_64Peri, NPoints_128Peri]
			Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, ndim=3, NumberConnections=NumConnections_list, FilamentLengths=FilLengths_list)
	
	if HOMEPC == 1:
		file_directory = '/mn/stornext/d5/aleh'
		savefile_directory = '/uio/hume/student-u70/aleh/Masters_project/disperse_results'
		
		if LCDM_model == 1:
			solve_file_dir = '/lcdm_testing/'
			solve_filename = 'lcdm_z0_test.solve'
			"""
			LCDM_z0_64_dir = 'lcdm_testing/LCDM_z0_64Periodic/'
			LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plots/', nPart=64)
			NConnections_64, FilLengths_64, FilPoints_64 = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDM64Periodic.a.NDskl', ndim=3)
			
			LCDM_z0_128_dir = 'lcdm_testing/LCDM_z0_128Periodic/'
			LCDM_z0_128Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128)
			NConnections_128, FilLengths_128, FilPoints_128 = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128Periodic.a.NDskl', ndim=3)
			"""
			"""
			LCDM_z0_256_dir = 'lcdm_testing/LCDM_z0_256Periodic/'
			LCDM_z0_256Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plots/', nPart=256)
			NConnections_256, FilLengths_256, FilPoints_256 = LCDM_z0_256Instance.Solve(LCDM_z0_256_dir+'SkelconvOutput_LCDM256Periodic.a.NDskl', ndim=3)
			
			LCDM_z0_512_dir = 'lcdm_testing/LCDM_z0_512Periodic/'
			LCDM_z0_512Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plots/', nPart=512)
			NConnections_512, FilLengths_512, FilPoints_512 = LCDM_z0_512Instance.Solve(LCDM_z0_512_dir+'SkelconvOutput_LCDM512Periodic.a.NDskl', ndim=3)
			"""
			LCDM_z0_64Test2_dir = 'lcdm_testing/LCDM_z0_64PeriodicTesting/'
			LCDM_z0_64Test2Instance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_64Test2_dir+'Plots/', nPart=64)
			NN, FF, FP = LCDM_z0_64Test2Instance.Solve(LCDM_z0_64Test2_dir+'SkelconvOutput_LCDMz064.a.NDskl', ndim=3)

			Comparison_dir = 'lcdm_testing/Comparison_plots/'
			if Comparison == 1:
				NumConnections_list = [NConnections_64, NConnections_128]
				FilLengths_list = [FilLengths_64, FilLengths_128]
				FilPoints_list = [FilPoints_64, FilPoints_128]
				#Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, ndim=3, NumberConnections=NumConnections_list, FilamentLengths=FilLengths_list)
				ComparisonInstance = Histogram_Comparison2(savefile=0, savefigDirectory=Comparison_dir, ndim=3, model='$\mathregular{\Lambda}$CDM',\
													 Nparticles=[64, 128])
				ComparisonInstance.Solve('Connections', NumConnections_list)
				ComparisonInstance.Solve('Lengths', FilLengths_list)
				ComparisonInstance.Solve('FilamentPoints', FilPoints_list)
				plt.show()