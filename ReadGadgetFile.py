import numpy as np

class Read_Gadget_file():
	def __init__(self, mask_check):
		self.read_file()
		self.Create_Mask()
		self.Create_KDTree()
		self.Mask_check_list = mask_check

		return self.PartPosX, self.PartPosY, self.PartPosZ, self.DM_tree

	def read_file(self):
		""" 
		Reads the gadget files from RAMSES 
		Code by Bridget Falck
		"""
		print 'Reading gadget files ...'
		numfiles = 512
		nparticles = 512
		datadir = '/mn/stornext/d5/astcosim/N512_B256_ISISCodePaper/lcdm/z0.000/data/gadgetfiles/'
		filename = datadir+'gadget.'

		self.PartPos = np.empty((nparticles**3,3),np.float32)
		self.PartVel = np.empty((nparticles**3,3),np.float32)
		self.PartIDs = np.empty((nparticles**3),np.int32)
		istart = 0
		for i in np.arange(0, numfiles):
		    file=filename+str(i)
		    f=open(file, 'rb')
		    
		    header_size = np.fromfile(f,np.int32,1)[0] # = 256: error catch here?
		    numpart = np.fromfile(f,np.int32,6)
		    npart = numpart[1] # number of particles in this file
		    mass = np.fromfile(f,np.float64,6)
		    pmass = mass[1] # in units of 10^10 solar masses?
		    scalefact,redshift = np.fromfile(f,np.float64,2)
		    flag_sfr,flag_feedback = np.fromfile(f,np.int32,2)
		    numpart_tot = np.fromfile(f,np.int32,6)
		    ntotal = numpart_tot[1]
		    flag_cooling,num_files = np.fromfile(f,np.int32,2)
		    boxsize,omega0,omegal,hubble = np.fromfile(f,np.float64,4)
		    flag_stellarage,flag_metals,hashtabsize = np.fromfile(f,np.int32,3)
		    # read rest of header_size + 2 dummy integers:
		    dummy = np.fromfile(f,np.int32,23)
		    
		    thispos = np.fromfile(f,np.float32,3*npart)
		    thispos = np.reshape(thispos, [npart, 3])
		    
		    # read velocities
		    dummy = np.fromfile(f,np.int32,2)
		    thisvel = np.fromfile(f,np.float32,3*npart)
		    thisvel = np.reshape(thisvel, [npart, 3])

		    # read IDs
		    dummy = np.fromfile(f,np.int32,2)
		    thisID = np.fromfile(f,np.int32,npart)
		    f.close()
		    
		    self.PartPos[istart:(istart+npart),:] = thispos
		    self.PartVel[istart:(istart+npart),:] = thisvel
		    self.PartIDs[istart:(istart+npart)] = thisID
		    istart = istart + npart

		print 'finished reading particles, '
		self.PartPos = self.PartPos/1000.0

		
	def Create_Mask(self):
		""" Creates a mask for the dark matter particles """
		MaskXdir = self.Mask_check_list[0]
		MaskYdir = self.Mask_check_list[1]
		MaskZdir = self.Mask_check_list[2]
		if not MaskXdir and not MaskYdir and MaskZdir:
			self.mask = np.logical_and(self.PartPos[:,2] < UpperBoundaryZDir, self.PartPos[:,2] > LowerBoundaryZDir)
		elif not MaskXdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(self.PartPos[:,0] < UpperBoundaryXDir, self.PartPos[:,0] > LowerBoundaryXDir)
		elif not MaskYdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(self.PartPos[:,1] < UpperBoundaryYDir, self.PartPos[:,1] > LowerBoundaryYDir)
		elif MaskXdir and not MaskYdir and MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.PartPos[:,2] < UpperBoundaryZDir, self.PartPos[:,2] > LowerBoundaryZDir),\
									   np.logical_and(self.PartPos[:,0] < UpperBoundaryXDir, self.PartPos[:,0] > LowerBoundaryXDir))
		elif MaskXdir and MaskYdir and not MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.PartPos[:,1] < UpperBoundaryYDir, self.PartPos[:,1] > LowerBoundaryYDir),\
									   np.logical_and(self.PartPos[:,0] < UpperBoundaryXDir, self.PartPos[:,0] > LowerBoundaryXDir))
		elif not MaskXdir and MaskYdir and MaskZdir:
			self.mask = np.logical_and(np.logical_and(self.PartPos[:,2] < UpperBoundaryZDir, self.PartPos[:,2] > LowerBoundaryZDir),\
									   np.logical_and(self.PartPos[:,1] < UpperBoundaryYDir, self.PartPos[:,1] > LowerBoundaryYDir))
		else:
			self.mask = False

		if self.mask is not False:
			self.PartPosX = self.PartPos[self.mask,0]
			self.PartPosY = self.PartPos[self.mask,1]
			self.PartPosZ = self.PartPos[self.mask,2]
		else:
			self.PartPosX = self.PartPos[:,0]
			self.PartPosY = self.PartPos[:,1]
			self.PartPosZ = self.PartPos[:,2]

	def Create_KDTree(self):
		""" Creates a KDTree of all the dark matter particles """
		DM_points = np.dstack((self.PartPos[:,0].ravel(), self.PartPos[:,1].ravel(), self.PartPos[:,2].ravel()))
		self.DM_tree = spatial.KDTree(DM_points[0])