import numpy as np
import os
import time
file_directory = '/mn/stornext/d5/aleh'

class FilamentsPerSigma():
	def __init__(self, filename):
		self.ReadFile(filename)
		self.Sort_filament_coordinates()
		self.Sort_filament_data()
		# Testing
		#self.Filaments_per_sigma(np.linspace(3,10,40))
		#self.Filaments_per_sigma2(np.linspace(4,10,5))

	def ReadFile(self, filename, dimensions=3):
		""" Reads the data from the skeleton file from Disperse """
		datafiles = open(os.path.join(file_directory, filename), 'r')
		#datafiles = open(filename, 'r')
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

	def Sort_filament_coordinates(self, dimensions=3):
		""" 
		Sorts the coordinates of the filaments and critical points to their respective arrays 
		Data to be sorted: Critical points coordinate, ID of filament and filament coordinates
		"""
		self.NFils = int(self.Filaments[0][0])

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
			self.CritPointXpos.append(float(Info[1]))
			self.CritPointYpos.append(float(Info[2]))
			self.CritPointZpos.append(float(Info[3]))
			#self.Neighbours_CP.append(int(Info[5]))	# Includes all CPs even without filaments
			Filament_connections = int(self.CriticalPoints[i+1][0])
			self.Number_filaments_connecting_to_CP.append(Filament_connections)
			Temp_filID = []
			Temp_CPID = []
			if Filament_connections != 0:
				self.Neighbours_CP.append(int(Info[5]))
				for j in range(1,Filament_connections+1):
					Temp_CPID.append(int(self.CriticalPoints[i+j+1][0]))
					Temp_filID.append(int(self.CriticalPoints[i+j+1][1]))
				self.Critpts_filamentID.append(np.array(Temp_filID))
				self.CP_id_of_connecting_filament.append(np.array(Temp_CPID))
			counter += 1
			i += 2 + Filament_connections

		self.CP_id_of_connecting_filament = np.asarray(self.CP_id_of_connecting_filament)
		self.Neighbours_CP = np.asarray(self.Neighbours_CP)

	def Sort_filament_data(self):
		""" 
		Sorts the data of the filaments and critical points.
		Data to be sorted includes: Persistence, persistence pairs, persistence ratios, field values, etc.
		"""
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

	def Filaments_per_sigma(self, sigma_array):
		"""
		Checks number of existing filaments based on sigma value
		This algorithm assumes that a critical point only has one filament connection 
		"""
		print 'Computing number of filaments as a function of sigma'
		fil_per_sig = []
		Temporary_sigmas = self.Persistence_nsigmas
		for sigmas in sigma_array:
			Filaments = []
			CPs_included = np.where(Temporary_sigmas >= sigmas)[0]
			time_start = time.clock()
			for i in CPs_included:
				for j in range(len(self.CP_id_of_connecting_filament[i])):
					if self.Neighbours_CP[i] == self.CP_id_of_connecting_filament[i][j]:
						Filaments.append(self.Critpts_filamentID[i][j])
			time_start = time.clock()
			Unique_filaments = np.unique(np.array(Filaments))
			fil_per_sig.append(len(Unique_filaments))
			Temporary_sigmas = Temporary_sigmas[CPs_included]
			self.Neighbours_CP = self.Neighbours_CP[CPs_included]

		return fil_per_sig

	def Filaments_per_sigma2(self, sigma_array):
		""" 
		Alternative filament per sigma estimator.
		This algorithm assumes that all filaments connecting to a critical point, which is to be removed, are removed as well
		"""
		print 'Computing number of filaments as a function of sigma'
		fil_per_sig = []
		Temporary_sigmas = self.Persistence_nsigmas
		for sigmas in sigma_array:
			CPs_included = np.where(Temporary_sigmas >= sigmas)[0]
			if len(CPs_included) != 0:
				fil_per_sig.append(np.sum(np.array(self.Number_filaments_connecting_to_CP)[CPs_included]))
			else:
				fil_per_sig.append(0)
		Temporary_sigmas = Temporary_sigmas[CPs_included]
		return fil_per_sig
		

if __name__ == '__main__':
	# Testing program
	FilamentsPerSigma('lcdm_testing/LCDM_z0_64PeriodicTesting/'+'SkelconvOutput_LCDMz064.a.NDskl')
	