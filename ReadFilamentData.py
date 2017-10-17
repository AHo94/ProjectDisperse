import numpy as np
import os

class read_disperse_output():
	""" 
	This class reads the skeleton output file from DisPerSE.
	All data sets are assumed to be 3D
	"""
	def __init__(self, directory, UnitConverter, FilamentLimit = False):
		self.directory = directory
		self.UnitConverter = UnitConverter
		self.dimensions = dimensions
		self.FilamentLimit = FilamentLimit

	def ReadFile(self, filename):
		""" Reads the data from the skeleton file from Disperse """
		time_start = time.clock()
		print 'Reading data for the file: ', filename, '...' 
		datafiles = open(os.path.join(self.directory, filename), 'r')
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
				self.xmin = float(MinValues[0])*self.UnitConverter
				self.xmax = float(MaxValues[0])*self.UnitConverter
				self.ymin = float(MinValues[1])*self.UnitConverter
				self.ymax = float(MaxValues[1])*self.UnitConverter
				self.zmin = float(MinValues[2])*self.UnitConverter
				self.zmax = float(MaxValues[2])*self.UnitConverter

			if '[CRITICAL POINTS]' in line:
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
					if self.FilamentLimit is not False:
						if SETLIMIT == self.FilamentLimit+1:
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

	def sort_cp_coordinates(self):
		""" Sort critical point coordinate data """
		print 'Sorting critical point coordinate data'
		CritPointXpos = []
		CritPointYpos = []
		CritPointZpos = []
		NcritPts = int(self.CriticalPoints[0][0])
		Critpts_filamentID = []
		CP_persistent_pair = []
		CP_id_of_connecting_filament = []
		Number_filaments_connecting_to_CP = []
		CP_type = []
		i = 1
		counter = 0
		while i < len(CriticalPoints):
			Info = CriticalPoints[i]
			"""
			# This saves all points
			CritPointXpos.append(float(Info[1])*self.UnitConverter)
			CritPointYpos.append(float(Info[2])*self.UnitConverter)
			CritPointZpos.append(float(Info[3])*self.UnitConverter)
			CP_persistent_pair.append(int(Info[5]))
			"""
			Filament_connections = int(CriticalPoints[i+1][0])
			Number_filaments_connecting_to_CP.append(Filament_connections)
			CP_type.append(int(Info[0]))
			Temp_filID = []
			Temp_CPID = []
			if Filament_connections != 0:
				""" Saves points which has a filament connection """
				CritPointXpos.append(float(Info[1])*self.UnitConverter)
				CritPointYpos.append(float(Info[2])*self.UnitConverter)
				CritPointZpos.append(float(Info[3])*self.UnitConverter)
				CP_persistent_pair.append(int(Info[5]))
				for j in range(1,Filament_connections+1):
					Temp_CPID.append(int(CriticalPoints[i+j+1][0]))
					Temp_filID.append(int(CriticalPoints[i+j+1][1]))
				Critpts_filamentID.append(np.array(Temp_filID))
				CP_id_of_connecting_filament.append(np.array(Temp_CPID))
			#else:
			#	break
			counter += 1
			i += 2 + Filament_connections

		Critpts_filamentID = np.array(Critpts_filamentID)
		CritPointXpos = np.asarray(CritPointXpos)
		CritPointYpos = np.asarray(CritPointYpos)
		CritPointZpos = np.asarray(CritPointZpos)
		return CritPointYpos, CritPointYpos, CritPointZpos, CP_type, CP_persistent_pair, Critpts_filamentID, CP_id_of_connecting_filament, Number_filaments_connecting_to_CP

	def sort_filament_coordinates(self):
		""" Sort filament coordinate data """
		print 'Sorting filament coordinate data'
		if self.FilamentLimit is not False:
			NFils = int(self.Filaments[0][0])
		else:
			NFils = FilamentLimit

		FilamentPos = []
		FilID = []
		PairIDS = []
		xdimPos = []
		ydimPos = []
		zdimPos = []
		NFilamentPoints = []
		k = 1
		NewID = 0
		for i in range(1, self.NFils):
			Filstuff = self.Filaments[k]
			TempPositions = []
			xtemp = []
			ytemp = []
			ztemp = []
			FilID.append(NewID)
			PairIDS.append(np.array([int(Filstuff[0]), int(Filstuff[1])]))
			for j in range(1, int(Filstuff[-1])+1):
				if k+j >= len(self.Filaments):
					break
				xPos = float(self.Filaments[k+j][0])*self.UnitConverter
				yPos = float(self.Filaments[k+j][1])*self.UnitConverter
				zPos = float(self.Filaments[k+j][2])*self.UnitConverter
				TempPositions.append(np.array([xPos, yPos]))
				xtemp.append(xPos)
				ytemp.append(yPos)
				ztemp.append(zPos)		
			NFilamentPoints.append(int(Filstuff[-1]))
			FilamentPos.append(np.array(TempPositions))
			xdimPos.append(np.array(xtemp))
			ydimPos.append(np.array(ytemp))
			zdimPos.append(np.array(ztemp))
			k += int(Filstuff[-1])+1
			NewID += 1
			if self.FilamentLimit is not False:
				if k >= self.FilamentLimit:
					NFils = len(xdimPos)
					break

		FilamentPos = np.array(FilamentPos)
		NFilamentPoints = np.array(NFilamentPoints)
		xdimPos = np.array(xdimPos)
		ydimPos = np.array(ydimPos)
		zdimPos = np.array(zdimPos)
		return FilamentPos, xdimPos, ydimPos, zdimPos, NFilamentPoints
		
	def get_data(self, filename):
		""" Recieves all info based on the filament input """
		self.ReadFile(filename)
		CP_coordinate_data = self.sort_cp_coordinates()
		Filament_coordinate_data = self.sort_filament_coordinates()

		#### In the main program, read these as single variables and have a function which unpacks them to self variables
		#### Save first to pickle and read it before unpacking!

		return CP_coordinate_data, Filament_coordinate_data