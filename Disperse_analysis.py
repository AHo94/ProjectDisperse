### Set comment on the two below to plot. Only use if running on papsukal, nekkar etc. 
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import PolyCollection
#from matplotlib.colors import colorConverter
import os
import scipy.stats as stats
from matplotlib import colors as mcolors

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
		""" Reads the data for a given file """
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

	def Sort_arrays(self, dimensions):
		""" Sorts data to their respective arrays """
		print 'Sorting data ...'
		self.NcritPts = int(self.CriticalPoints[0][0])
		if FilamentLimit == 0:
			self.NFils = int(self.Filaments[0][0])
		else:
			self.NFils = FilamentLimit

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

	def Check_Boundaries(self):
		""" Checks the boundaries. Assuming periodic bondary. Will split segment if change is larger than 10000 """
		print 'Checking boundaries...'
		self.FilXYPositionsBoundary = []
		self.FilZPositionsBoundary = []
		XY_positions = []
		for i in range(self.NFils-1):
			check = 0
			xIndex = np.where(np.abs(np.array(self.xdimPos[i][:-1]) - np.array(self.xdimPos[i][1:])) > 10000.0)[0]
			yIndex = np.where(np.abs(np.array(self.ydimPos[i][:-1]) - np.array(self.ydimPos[i][1:])) > 10000.0)[0]
			zIndex = np.where(np.abs(np.array(self.zdimPos[i][:-1]) - np.array(self.zdimPos[i][1:])) > 10000.0)[0]

			if len(xIndex) == 2 and len(yIndex) == 0 and len(zIndex) == 0:
				x_upper = self.xdimPos[i][0:xIndex[0]+1]
				x_mid = self.xdimPos[i][xIndex[0]+1:xIndex[1]+1]
				x_lower = self.xdimPos[i][xIndex[1]+1:]
				
				y_upper = self.ydimPos[i][0:xIndex[0]+1]
				y_mid = self.ydimPos[i][xIndex[0]+1:xIndex[1]+1]
				y_lower = self.ydimPos[i][xIndex[1]+1:]
				
				z_upper = self.zdimPos[i][0:xIndex[0]+1]
				z_mid = self.zdimPos[i][xIndex[0]+1:xIndex[1]+1]
				z_lower = self.zdimPos[i][xIndex[1]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(yIndex) == 2 and len(xIndex) == 0 and len(zIndex) == 0:
				x_upper = self.xdimPos[i][0:yIndex[0]+1]
				x_mid = self.xdimPos[i][yIndex[0]+1:yIndex[1]+1]
				x_lower = self.xdimPos[i][yIndex[1]+1:]
				
				y_upper = self.ydimPos[i][0:yIndex[0]+1]
				y_mid = self.ydimPos[i][yIndex[0]+1:yIndex[1]+1]
				y_lower = self.ydimPos[i][yIndex[1]+1:]
				
				z_upper = self.zdimPos[i][0:yIndex[0]+1]
				z_mid = self.zdimPos[i][yIndex[0]+1:yIndex[1]+1]
				z_lower = self.zdimPos[i][yIndex[1]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)
			
			elif len(zIndex) == 2 and len(yIndex) == 0 and len(xIndex) == 0:
				x_upper = self.xdimPos[i][0:zIndex[0]+1]
				x_mid = self.xdimPos[i][zIndex[0]+1:zIndex[1]+1]
				x_lower = self.xdimPos[i][zIndex[1]+1:]
				
				y_upper = self.ydimPos[i][0:zIndex[0]+1]
				y_mid = self.ydimPos[i][zIndex[0]+1:zIndex[1]+1]
				y_lower = self.ydimPos[i][zIndex[1]+1:]
				
				z_upper = self.zdimPos[i][0:zIndex[0]+1]
				z_mid = self.zdimPos[i][zIndex[0]+1:zIndex[1]+1]
				z_lower = self.zdimPos[i][zIndex[1]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(xIndex) == 1 and len(yIndex) == 0 and len(zIndex) == 0:
				x_upper = self.xdimPos[i][0:xIndex[0]+1]
				x_lower = self.xdimPos[i][xIndex[0]+1:]
				
				y_upper = self.ydimPos[i][0:xIndex[0]+1]
				y_lower = self.ydimPos[i][xIndex[0]+1:]
				
				z_upper = self.zdimPos[i][0:xIndex[0]+1]
				z_lower = self.zdimPos[i][xIndex[0]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(yIndex) == 1 and len(xIndex) == 0 and len(zIndex) == 0:
				x_upper = self.xdimPos[i][0:yIndex[0]+1]
				x_lower = self.xdimPos[i][yIndex[0]+1:]
				
				y_upper = self.ydimPos[i][0:yIndex[0]+1]
				y_lower = self.ydimPos[i][yIndex[0]+1:]

				z_upper = self.zdimPos[i][0:yIndex[0]+1]
				z_lower = self.zdimPos[i][yIndex[0]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(zIndex) == 1 and len(yIndex) == 0 and len(xIndex) == 0:
				x_upper = self.xdimPos[i][0:zIndex[0]+1]
				x_lower = self.xdimPos[i][zIndex[0]+1:]
				
				y_upper = self.ydimPos[i][0:zIndex[0]+1]
				y_lower = self.ydimPos[i][zIndex[0]+1:]
				
				z_upper = self.zdimPos[i][0:zIndex[0]+1]
				z_lower = self.zdimPos[i][zIndex[0]+1:]
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(xIndex) > 0 and len(yIndex) > 0:
				if xIndex[0] < yIndex[0]:
					x_upper = self.xdimPos[i][0:xIndex[0]+1]
					x_mid = self.xdimPos[i][xIndex[0]+1:yIndex[0]+1]
					x_lower = self.xdimPos[i][yIndex[0]+1:]

					y_upper = self.ydimPos[i][0:xIndex[0]+1]
					y_mid = self.ydimPos[i][xIndex[0]+1:yIndex[0]+1]
					y_lower = self.ydimPos[i][yIndex[0]+1:]

					z_upper = self.zdimPos[i][0:xIndex[0]+1]
					z_mid = self.zdimPos[i][xIndex[0]+1:yIndex[0]+1]
					z_lower = self.zdimPos[i][yIndex[0]+1:]
				else:
					x_upper = self.xdimPos[i][0:yIndex[0]+1]
					x_mid = self.xdimPos[i][yIndex[0]+1:xIndex[0]+1]
					x_lower = self.xdimPos[i][xIndex[0]+1:]

					y_upper = self.ydimPos[i][0:yIndex[0]+1]
					y_mid = self.ydimPos[i][yIndex[0]+1:xIndex[0]+1]
					y_lower = self.ydimPos[i][xIndex[0]+1:]

					z_upper = self.zdimPos[i][0:yIndex[0]+1]
					z_mid = self.zdimPos[i][yIndex[0]+1:xIndex[0]+1]
					z_lower = self.zdimPos[i][xIndex[0]+1:]
					
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(xIndex) > 0 and len(zIndex) > 0:				
				if xIndex[0] < zIndex[0]:
					x_upper = self.xdimPos[i][0:xIndex[0]+1]
					x_mid = self.xdimPos[i][xIndex[0]+1:zIndex[0]+1]
					x_lower = self.xdimPos[i][zIndex[0]+1:]

					y_upper = self.ydimPos[i][0:xIndex[0]+1]
					y_mid = self.ydimPos[i][xIndex[0]+1:zIndex[0]+1]
					y_lower = self.ydimPos[i][zIndex[0]+1:]

					z_upper = self.zdimPos[i][0:xIndex[0]+1]
					z_mid = self.zdimPos[i][xIndex[0]+1:zIndex[0]+1]
					z_lower = self.zdimPos[i][zIndex[0]+1:]
				else:
					x_upper = self.xdimPos[i][0:zIndex[0]+1]
					x_mid = self.xdimPos[i][zIndex[0]+1:xIndex[0]+1]
					x_lower = self.xdimPos[i][xIndex[0]+1:]

					y_upper = self.ydimPos[i][0:zIndex[0]+1]
					y_mid = self.ydimPos[i][zIndex[0]+1:xIndex[0]+1]
					y_lower = self.ydimPos[i][xIndex[0]+1:]

					z_upper = self.zdimPos[i][0:zIndex[0]+1]
					z_mid = self.zdimPos[i][zIndex[0]+1:xIndex[0]+1]
					z_lower = self.zdimPos[i][xIndex[0]+1:]
					
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)

			elif len(yIndex) > 0 and len(zIndex) > 0:
				if zIndex[0] < yIndex[0]:
					x_upper = self.xdimPos[i][0:zIndex[0]+1]
					x_mid = self.xdimPos[i][zIndex[0]+1:yIndex[0]+1]
					x_lower = self.xdimPos[i][yIndex[0]+1:]

					y_upper = self.ydimPos[i][0:zIndex[0]+1]
					y_mid = self.ydimPos[i][zIndex[0]+1:yIndex[0]+1]
					y_lower = self.ydimPos[i][yIndex[0]+1:]

					z_upper = self.zdimPos[i][0:zIndex[0]+1]
					z_mid = self.zdimPos[i][zIndex[0]+1:yIndex[0]+1]
					z_lower = self.zdimPos[i][yIndex[0]+1:]
				else:
					x_upper = self.xdimPos[i][0:yIndex[0]+1]
					x_mid = self.xdimPos[i][yIndex[0]+1:zIndex[0]+1]
					x_lower = self.xdimPos[i][zIndex[0]+1:]

					y_upper = self.ydimPos[i][0:yIndex[0]+1]
					y_mid = self.ydimPos[i][yIndex[0]+1:zIndex[0]+1]
					y_lower = self.ydimPos[i][zIndex[0]+1:]

					z_upper = self.zdimPos[i][0:yIndex[0]+1]
					z_mid = self.zdimPos[i][yIndex[0]+1:zIndex[0]+1]
					z_lower = self.zdimPos[i][zIndex[0]+1:]
					
				XY_positions.append([x_upper, y_upper])
				XY_positions.append([x_mid, y_mid])
				XY_positions.append([x_lower, y_lower])
				self.FilZPositionsBoundary.append(z_upper)
				self.FilZPositionsBoundary.append(z_mid)
				self.FilZPositionsBoundary.append(z_lower)

			else:
				XY_positions.append([self.xdimPos[i][:], self.ydimPos[i][:]])
				self.FilZPositionsBoundary.append(self.zdimPos[i][:])
			
		for i in range(len(XY_positions)):
			TempPos = []
			for j in range(len(XY_positions[i][0])):
				TempPos.append([XY_positions[i][0][j], XY_positions[i][1][j]])
			self.FilXYPositionsBoundary.append(TempPos)
		# Remove single points in the arrays
		Indices_to_be_removed = []
		for j in range(len(self.FilXYPositionsBoundary)-1):
			if len(self.FilXYPositionsBoundary[j]) <= 1:
				Indices_to_be_removed.append(j)

		for k in range(len(Indices_to_be_removed)-1):
			del self.FilXYPositionsBoundary[Indices_to_be_removed[k]-k]
			del self.FilZPositionsBoundary[Indices_to_be_removed[k]-k]


		for i in range(len(self.FilXYPositionsBoundary)-1):
			if len(self.FilXYPositionsBoundary[i]) < 1:
				del self.FilXYPositionsBoundary[i]
				del self.FilZPositionsBoundary[i]

	def Filament_Length(self, dimensions):
		if dimensions == 3:
			self.FilLengths = []
			for i in range(self.NFils-1):
				TempLength = 0
				for j in range(len(self.xdimPos[i])-1):
					TempLength += np.sqrt((self.xdimPos[i][j+1] - self.xdimPos[i][j])**2.0 + (self.ydimPos[i][j+1] - self.ydimPos[i][j])**2.0\
								 + (self.zdimPos[i][j+1] - self.zdimPos[i][j])**2.0)
				self.FilLengths.append(TempLength)

	def Plot_Figures(self, filename, ndim=2):
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
			ColorArray = np.linspace(0,10,110)
		else:
			ColorArray = np.linspace(0,10,1)

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

		# Histogram of filament lengths
		LenMin = min(self.FilLengths)
		LenMax = max(self.FilLengths)
		hist = np.histogram(self.FilLengths)
		BinList_FilLengths = np.linspace(LenMin, LenMax, len(hist[0])/0.0002 + 1)
		#BinSize_FilLengths = (LenMax - LenMin)/(0.5) + 1
		#BinList_FilLengths = np.linspace(LenMin, LenMax, BinSize_FilLengths)
		#self.SortedLengths = sorted(self.FilLengths)
		self.FilLengths.sort()
		fit = stats.norm.pdf(self.FilLengths, np.mean(self.FilLengths), np.std(self.FilLengths))
		FilamentLengthsHist = plt.figure()
		#plt.hist(self.FilLengths, align='left', rwidth=1, bins=BinList_FilLengths, normed=NormedArgument)#, bins=BinList_FilLengths)
		plt.hold("on")
		plt.plot(self.FilLengths, fit)
		plt.xlabel('Length of filaments - [Mpc]')
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

		if ndim == 2:
			if PlotFilaments == 1:
				FilPositions = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection(line_segments)
				plt.xlabel('$\mathregular{x}$ - Mpc')
				plt.ylabel('$\mathregular{y}$ - Mpc')
				plt.title('Positions of the filaments with '+ self.nPart_text+ '$^3$ particles')
			if PlotFilamentsWCritPts == 1:
				FilPositions_WCritPts = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection(line_segments)
				plt.hold("on")
				plt.plot(self.CritPointXpos, self.CritPointYpos, 'ro', alpha=0.7, markersize=3)
				#plt.plot(self.CritPointXposNOTCON, self.CritPointYposNOTCON, 'go', alpha=0.7, markersize=3)
				plt.xlabel('$\mathregular{x}$ - Mpc')
				plt.ylabel('$\mathregular{y}$ - Mpc')
				plt.title('Position of the filaments with critical points shown')
		if ndim == 3:
			if PlotFilaments == 1:
				FilPositions = plt.figure()
				ax = FilPositions.gca(projection='3d')
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				#line_segments = LineCollection(self.FilXYPositionsBoundary, linestyle='solid')
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				#ax.add_collection3d(line_segments, self.FilZPositionsBoundary, zdir='z')
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				ax.set_xlabel('$\mathregular{x}$ - Mpc')
				ax.set_ylabel('$\mathregular{y}$ - Mpc')
				ax.set_zlabel('$\mathregular{z}$ - Mpc')
				plt.title('3D Position of the filaments with '+ self.nPart_text+ '$^3$ particles.')
			if PlotFilamentsWCritPts == 1:
				FilPositions_WCritPts = plt.figure()
				ax = FilPositions_WCritPts.gca(projection='3d')
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				ax.set_zlim(self.zmin, self.zmax)
				line_segments = LineCollection(self.FilamentPos, linestyle='solid', array=ColorArray, cmap=plt.cm.gist_ncar)
				ax.add_collection3d(line_segments, self.zdimPos, zdir='z')
				plt.hold("on")
				ax.plot(self.CritPointXpos, self.CritPointYpos, self.CritPointZpos, 'ro', alpha=0.7, markersize=3)
				ax.set_xlabel('$\mathregular{x}$ - Mpc')
				ax.set_ylabel('$\mathregular{y}$ - Mpc')
				ax.set_zlabel('$\mathregular{z}$ - Mpc')
				plt.title('3D Position of the filaments with critical points')
			if Projection2D == 1:
				FilPositions_2DProjection = plt.figure()
				ax = plt.axes()
				ax.set_xlim(self.xmin, self.ymax)
				ax.set_ylim(self.ymin, self.ymax)
				line_segments = LineCollection(self.FilamentPos, array=ColorArray, cmap=plt.cm.gist_ncar)
				#plt.gca().add_collection(line_segments)
				ax.autoscale()
				ax.add_collection(line_segments)
				ax.set_xlabel('$\mathregular{x}$ - Mpc')
				ax.set_ylabel('$\mathregular{y}$ - Mpc')
				plt.title('2D Position projection of the filaments')				

		if self.savefile == 1:
			ConnectedHist.savefig(self.results_dir + 'NumberFilamentConnectedHistogram' + Filenamedimension)
			FilamentLengthsHist.savefig(self.results_dir + 'FilamentLengthsHistogram' + Filenamedimension)
			FilamentPtsHis.savefig(self.results_dir + 'FilamentPointsHistogram' + Filenamedimension)
			if PlotFilaments == 1:
				FilPositions.savefig(self.results_dir + 'FilamentPositions' + Filenamedimension)
			if PlotFilamentsWCritPts == 1:
				FilPositions_WCritPts.savefig(self.results_dir + 'FilamentPositionsWithCritPts' + Filenamedimension)
			if Projection2D == 1:
				FilPositions_2DProjection.savefig(self.results_dir + '2DProjectedFilamentPosition' + Filenamedimension)
		else:
			plt.show()

	def Solve(self, filename, ndim=2):
		self.ReadFile(filename, ndim)
		self.Sort_arrays(ndim)
		self.Filament_Length(ndim)
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
			plt.hist(self.FilamentLengths[i], align='mid', rwidth=1, bins=Bin_list, normed=True, alpha=alphas[i])
		plt.xlabel('Filament lengths -[Mpc]')
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
			fit = stats.norm.pdf(histograms[i], np.mean(histograms[i]), np.std(histograms[i]))
			#plt.hist(histograms[i], align='mid', rwidth=1, bins=300, normed=True, alpha=self.alphas[i])
			plt.plot(histograms[i], fit)
		plt.xlabel('Filament lengths -[Mpc]')
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
	PlotFilaments = 1			# Set 1 to plot actual filaments
	PlotFilamentsWCritPts = 0	# Set to 1 to plot filaments with critical points
	Projection2D = 1 			# Set to 1 to plot a 2D projection of the 3D case
	FilamentColors = 1 			# Set to 1 to get different colors for different filaments
	Comparison = 1				# Set 1 if you want to compare different number of particles. Usual plots will not be plotted!
	FilamentLimit = 0			# Limits the number of filament read from file. Reads all if 0

	if HOMEPC == 0:
		file_directory = 'C:/Users/Alex/Documents/Masters_project/Disperse'
		savefile_directory = file_directory
		#solveInstance1 = Disperse_Plotter(savefile=1, savefigDirectory='Plot_Disperse_Example/', nPart=64)
		#solveInstance1.Plot("simu_2D.ND.NDnet_s3.up.NDskl.a.NDskl", ndim=2)
		#solveInstance1.Plot("simu_32_id.gad.NDnet_s3.5.up.NDskl.a.NDskl", ndim=3)
		
		LCDM_64Periodic_dir = 'lcdm_z0_testing/LCDM64_Periodic/'
		LCDM_z0_64Peri = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_64Periodic_dir + 'Plots/', nPart=64)
		NConnections_64Peri, FilLengths_64Peri, NPoints_64Peri = LCDM_z0_64Peri.Solve(LCDM_64Periodic_dir+'SkelconvOutput_LCDM64Periodic.a.NDskl', ndim=3)
		
		#LCDM_128Periodic_dir = 'lcdm_z0_testing/LCDM128_Periodic/'
		#LCDM_z0_128Peri = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_128Periodic_dir + 'Plots/', nPart=128)
		#NConnections_128Peri, FilLengths_128Peri, NPoints_128Peri = LCDM_z0_128Peri.Solve(LCDM_128Periodic_dir+'SkelconvOutput_LCDM128Periodic.a.NDskl', ndim=3)
		
		# Lots of memory usage for 256^3.
		#LCDM_256Periodic_dir = 'lcdm_z0_testing/LCDM256_Periodic/'
		#LCDM_z0_256Peri = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_256Periodic_dir + 'Plots/', nPart=256)
		#LCDM_z0_256Peri.Solve(LCDM_256Periodic_dir+'SkelconvOutput_LCDM256Periodic.a.NDskl', ndim=3)
		
		Comparison_dir = 'lcdm_z0_testing/Comparison_plots/'
		if Comparison == 1:
			NumConnections_list = [NConnections_64Peri, NConnections_128Peri]
			FilLengths_list = [FilLengths_64Peri, FilLengths_128Peri]
			FilPoints_list = [NPoints_64Peri, NPoints_128Peri]
			Histogram_Comparison(savefile=1, savefigDirectory=Comparison_dir, ndim=3, NumberConnections=NumConnections_list, FilamentLengths=FilLengths_list)
	
	if HOMEPC == 1:
		file_directory = '/mn/stornext/d5/aleh'
		savefile_directory = '/uio/hume/student-u70/aleh/Masters_project/disperse_results'
		
		LCDM_z0_64_dir = 'lcdm_testing/LCDM_z0_64Periodic/'
		LCDM_z0_64Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_64_dir+'Plots/', nPart=64)
		NConnections_64, FilLengths_64, FilPoints_64 = LCDM_z0_64Instance.Solve(LCDM_z0_64_dir+'SkelconvOutput_LCDM64Periodic.a.NDskl', ndim=3)
		
		LCDM_z0_128_dir = 'lcdm_testing/LCDM_z0_128Periodic/'
		LCDM_z0_128Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_128_dir+'Plots/', nPart=128)
		NConnections_128, FilLengths_128, FilPoints_128 = LCDM_z0_128Instance.Solve(LCDM_z0_128_dir+'SkelconvOutput_LCDM128Periodic.a.NDskl', ndim=3)
		
		"""
		LCDM_z0_256_dir = 'lcdm_testing/LCDM_z0_256Periodic/'
		LCDM_z0_256Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_256_dir+'Plots/', nPart=256)
		NConnections_256, FilLengths_256, FilPoints_256 = LCDM_z0_256Instance.Solve(LCDM_z0_256_dir+'SkelconvOutput_LCDM256Periodic.a.NDskl', ndim=3)
		
		LCDM_z0_512_dir = 'lcdm_testing/LCDM_z0_512Periodic/'
		LCDM_z0_512Instance = Disperse_Plotter(savefile=1, savefigDirectory=LCDM_z0_512_dir+'Plots/', nPart=512)
		NConnections_512, FilLengths_512, FilPoints_512 = LCDM_z0_512Instance.Solve(LCDM_z0_512_dir+'SkelconvOutput_LCDM512Periodic.a.NDskl', ndim=3)
		"""
		"""
		LCDM_z0_64_rockstar_dir = 'lcdm_rockstar/LCDM_z0_64/'
		LCDM_z0_64rockstarInstance = Disperse_Plotter(savefile=0, savefigDirectory=LCDM_z0_64_rockstar_dir+'Plots/', nPart=64)
		NN, FF, FP = LCDM_z0_64rockstarInstance.Solve(LCDM_z0_64_rockstar_dir+'SkelconvOutput_LCDMz0_64.a.NDskl', ndim=3)
		"""
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