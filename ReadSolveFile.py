
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
