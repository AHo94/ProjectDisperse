import numpy as np
import os
import time

class read_disperse_output():
	""" 
	This class reads the skeleton output file from DisPerSE.
	All data sets are assumed to be 3D.
	"""
	def __init__(self, directory, UnitConverter, FilamentLimit=False):
		self.directory = directory
		self.UnitConverter = np.float32(UnitConverter)
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
				self.xmin = np.float32(MinValues[0])*self.UnitConverter
				self.xmax = np.float32(MaxValues[0])*self.UnitConverter
				self.ymin = np.float32(MinValues[1])*self.UnitConverter
				self.ymax = np.float32(MaxValues[1])*self.UnitConverter
				self.zmin = np.float32(MinValues[2])*self.UnitConverter
				self.zmax = np.float32(MaxValues[2])*self.UnitConverter

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
					if self.FilamentLimit:
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
		NcritPts = np.int16(self.CriticalPoints[0][0])
		Critpts_filamentID = []
		CP_persistent_pair = []
		CP_id_of_connecting_filament = []
		Number_filaments_connecting_to_CP = []
		CP_type = []
		i = 1
		counter = 0
		while i < len(self.CriticalPoints):
			Info = self.CriticalPoints[i]
			"""
			# This saves all points
			CritPointXpos.append(float(Info[1])*self.UnitConverter)
			CritPointYpos.append(float(Info[2])*self.UnitConverter)
			CritPointZpos.append(float(Info[3])*self.UnitConverter)
			CP_persistent_pair.append(int(Info[5]))
			"""
			Filament_connections = np.int16(self.CriticalPoints[i+1][0])
			Number_filaments_connecting_to_CP.append(Filament_connections)
			CP_type.append(np.int16(Info[0]))
			Temp_filID = []
			Temp_CPID = []
			if Filament_connections != 0:
				""" Saves points which has a filament connection """
				CritPointXpos.append(np.float16(Info[1])*self.UnitConverter)
				CritPointYpos.append(np.float16(Info[2])*self.UnitConverter)
				CritPointZpos.append(np.float16(Info[3])*self.UnitConverter)
				CP_persistent_pair.append(int(Info[5]))
				for j in range(1,Filament_connections+1):
					Temp_CPID.append(np.int16(self.CriticalPoints[i+j+1][0]))
					Temp_filID.append(np.int16(self.CriticalPoints[i+j+1][1]))
				Critpts_filamentID.append(np.array(Temp_filID, dtype=np.int16))
				CP_id_of_connecting_filament.append(np.array(Temp_CPID, dtype=np.int16))
			#else:
			#	break
			counter += 1
			i += 2 + Filament_connections
		#Critpts_filamentID = np.asarray(Critpts_filamentID, dtype=np.int16)
		CritPointXpos = np.asarray(CritPointXpos, dtype=np.float16)
		CritPointYpos = np.asarray(CritPointYpos, dtype=np.float16)
		CritPointZpos = np.asarray(CritPointZpos, dtype=np.float16)
		return_list = [
				CritPointXpos, CritPointYpos, CritPointZpos, 
				CP_type, CP_persistent_pair, Critpts_filamentID, 
				CP_id_of_connecting_filament, Number_filaments_connecting_to_CP
				]

		self.Critpts_filamentID_array = Critpts_filamentID
		self.CP_type = CP_type
		return return_list
		
	def sort_filament_coordinates(self):
		""" Sort filament coordinate data """
		print 'Sorting filament coordinate data'
		if not self.FilamentLimit:
			self.NFils = np.int32(self.Filaments[0][0])
		else:
			self.NFils = self.FilamentLimit

		FilamentPos = []
		FilID = []
		PairIDS = []
		xdimPos = []
		ydimPos = []
		zdimPos = []
		NFilamentPoints = []
		k = 1
		NewID = 0
		for i in range(1, self.NFils+1):
			Filstuff = self.Filaments[k]
			TempPositions = []
			xtemp = []
			ytemp = []
			ztemp = []
			FilID.append(NewID)
			PairIDS.append(np.array([np.int16(Filstuff[0]), np.int16(Filstuff[1])], dtype=np.int16))
			for j in range(1, int(Filstuff[-1])+1):
				if k+j >= len(self.Filaments):
					break
				xPos = np.float32(self.Filaments[k+j][0])*self.UnitConverter
				yPos = np.float32(self.Filaments[k+j][1])*self.UnitConverter
				zPos = np.float32(self.Filaments[k+j][2])*self.UnitConverter
				TempPositions.append(np.array([xPos, yPos],dtype=np.float32))
				xtemp.append(xPos)
				ytemp.append(yPos)
				ztemp.append(zPos)		
			NFilamentPoints.append(np.int16(Filstuff[-1]))
			FilamentPos.append(np.array(TempPositions, dtype=np.float32))
			xdimPos.append(np.array(xtemp, dtype=np.float32))
			ydimPos.append(np.array(ytemp, dtype=np.float32))
			zdimPos.append(np.array(ztemp, dtype=np.float32))
			k += int(Filstuff[-1])+1
			NewID += 1
			if self.FilamentLimit:
				if k >= self.FilamentLimit:
					self.NFils = len(xdimPos)
					break

		FilamentPos = np.array(FilamentPos)
		NFilamentPoints = np.array(NFilamentPoints, dtype=np.int32)
		xdimPos = np.array(xdimPos)
		ydimPos = np.array(ydimPos)
		zdimPos = np.array(zdimPos)
		self.PairIDS = PairIDS
		return_list = [self.NFils, FilamentPos, xdimPos, ydimPos, zdimPos, NFilamentPoints, FilID, PairIDS]
		return return_list

	def sort_cp_data(self):
		""" Sorts data under the section CRITICAL POINT DATA in the skeleton file """
		print 'Sorting data of filaments and critical points'
		header_count_critpoint = int(self.CriticalPointsData[0][0])
		Persistence_ratio = []
		Persistence_nsigmas = []
		Persistence = []
		Persistence_pair = []
		Parent_index = []
		Parent_log_index = []
		log_CritPoint_field_value = []
		CritPoint_field_value = []
		CritPoint_cell = []

		for i in range(1, len(self.CriticalPointsData)):
			if i <= header_count_critpoint:
				continue
			else:
				if i <= len(self.Critpts_filamentID_array)+header_count_critpoint:
					# Do not include critical points that are not connected to filaments
					Persistence_ratio.append(np.float32(self.CriticalPointsData[i][0]))
					Persistence_nsigmas.append(np.float32(self.CriticalPointsData[i][1]))
					Persistence.append(np.float32(self.CriticalPointsData[i][2]))
					Persistence_pair.append(np.int32(self.CriticalPointsData[i][3]))
					Parent_index.append(np.int32(self.CriticalPointsData[i][4]))
					Parent_log_index.append(np.int32(self.CriticalPointsData[i][5]))
					log_CritPoint_field_value.append(np.float32(self.CriticalPointsData[i][6]))
					CritPoint_field_value.append(np.float32(self.CriticalPointsData[i][7]))
					CritPoint_cell.append(np.float32(self.CriticalPointsData[i][8]))
		
		return_list = [
				Persistence_ratio, Persistence_nsigmas, Persistence, 
				Persistence_pair, Parent_index, Parent_log_index, 
				log_CritPoint_field_value, CritPoint_field_value, CritPoint_cell
				]
		self.Persistence_nsigmas_array = np.asarray(Persistence_nsigmas, dtype=np.float32)
		return return_list

	def sort_filament_data(self):
		""" Sorts data under the section FILAMENT DATA of the skeleton file """
		header_count_filament = int(self.FilamentsData[0][0])
		Filament_field_value = []
		orientation = []
		Filament_cell = []
		log_Filament_field_value = []
		Filament_type = []

		for j in range(1, len(self.FilamentsData)):
			if j <= header_count_filament:
				continue
			else:
				Filament_field_value.append(np.float32(self.FilamentsData[j][0]))
				orientation.append(np.int32(self.FilamentsData[j][1]))
				Filament_cell.append(np.float32(self.FilamentsData[j][2]))
				log_Filament_field_value.append(np.float32(self.FilamentsData[j][3]))
				Filament_type.append(np.int32(self.FilamentsData[j][4]))

		# Give filaments persistence nsigma value based on the nsigma value of the connected maxima CP.
		maximas = np.where(np.array(self.CP_type) == 3)[0]
		maxima_nsigmas = self.Persistence_nsigmas_array[maximas]
		Filament_sigma = np.zeros(self.NFils-1)
		for i in range(self.NFils-1):
			CP_maxima_index = self.PairIDS[i][1]
			Filament_sigma[i] = self.Persistence_nsigmas_array[CP_maxima_index]

		return_list = [
				Filament_field_value, orientation, Filament_cell,
				log_Filament_field_value, Filament_type, Filament_sigma
				]
		return return_list

	def get_data(self, filename):
		""" Recieves all info based on the filament input """
		self.ReadFile(filename)
		CP_coordinate_data = self.sort_cp_coordinates()
		Filament_coordinate_data = self.sort_filament_coordinates()
		CP_data = self.sort_cp_data()
		Filament_data = self.sort_filament_data()
		Box_boundaries = [self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax]
		#### In the main program, read these as single variables and have a function which unpacks them to self variables
		#### Save first to pickle and read it before unpacking!
		return Box_boundaries, CP_coordinate_data, Filament_coordinate_data, CP_data, Filament_data

if __name__ == '__main__':
	# Testing program
	Instance = read_disperse_output('C:/Users/Alex/Documents/Masters_project/Disperse', 1)
	a,b,c,d, e = Instance.get_data('lcdm_z0_testing/LCDM_z0_64PeriodicTesting/SkelconvOutput_LCDMz064.a.NDskl')
	# At UiO 
	#Instance = read_disperse_output('/mn/stornext/d5/aleh', 1)
	#a,b,c,d,e = Instance.get_data('lcdm_testing/LCDM_z0_64PeriodicTesting/SkelconvOutput_LCDMz064.a.NDskl')