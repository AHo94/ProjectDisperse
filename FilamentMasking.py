import numpy as np
import time

class FilamentMasking():
	"""
	This class masks the filament given a masking boundary.
	Filaments within the given boundary are included, while the rest are ignored.
	Also includes filaments which gets 'cut off' outside the boundary.
	"""
	def __init__(self, FilamentPos, xpoints, ypoints, zpoints, lengths, NumFilaments, mask_dir, boundaries):
		self.FilamentPos = FilamentPos
		self.xdimPos = xpoints
		self.ydimPos = ypoints
		self.zdimPos = zpoints
		self.NFils = NumFilaments
		self.Boundary_list = boundaries
		self.Masking_directions = mask_dir
		self.LengthSplitFilament = lengths

	def Mask_slices(self):
		"""
		Creates a mask to plot a slice of the filament box. Boundary of slice chosen arbitrarily.
		The masking includes filaments that are within the given boundary.
		Also includes masking where the filament are cut off outside the boundary. 
		Filament length input must be an array that is of equal size as the position arrays.
		If filaments are split, then the filament length array must contain identical lengths for the split filament.
		e.g filament I is split into two parts with length 2, then LengSplitFilament = [..., 2, 2, ...] 
		"""
		UpperBoundaryXDir, UpperBoundaryYDir, UpperBoundaryZDir, LowerBoundaryXDir, LowerBoundaryYDir, LowerBoundaryZDir = self.Boundary_list
		MaskXdir, MaskYdir, MaskZdir = self.Masking_directions 

		time_start = time.clock()
		print 'Computing masks'
		self.zdimMasked = []
		self.MaskedFilamentSegments = []
		self.MaskedLengths = []

		self.CutOffFilamentSegments = []
		self.CutOffzDim = []
		self.CutOffLengths = []

		def Add_filament(index):
			self.zdimMasked.append(self.zdimPos[index])
			self.MaskedFilamentSegments.append(self.FilamentPos[index])
			self.MaskedLengths.append(self.LengthSplitFilament[index])

		def Cutoff_filaments(Indices_list, index):
			FilSegmentTemp = []
			zDimTemp = []
			for idx in Indices_list:
				zDimTemp.append(self.zdimPos[index][idx])
				FilSegmentTemp.append(np.array([self.xdimPos[index][idx], self.ydimPos[index][idx]]))
			self.CutOffFilamentSegments.append(FilSegmentTemp)
			self.CutOffzDim.append(zDimTemp)
			TempLen = 0
			for j in range(len(zDimTemp)-1):
				TempLen += np.sqrt((FilSegmentTemp[j+1][0] - FilSegmentTemp[j][0])**2 + (FilSegmentTemp[j+1][1] - FilSegmentTemp[j][1])**2\
						 + (zDimTemp[j+1] - zDimTemp[j])**2)
			self.CutOffLengths.append(TempLen)
		
		for i in range(len(self.xdimPos)):
			if MaskXdir and not MaskYdir and not MaskZdir:
				Indices = np.where(np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)))[0]
				Segment_check = (np.array(self.xdimPos[i]) > LowerBoundaryXDir).any() and (np.array(self.xdimPos[i]) < UpperBoundaryXDir).any()
			elif not MaskXdir and MaskYdir and not MaskZdir:
				Indices = np.where(np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir)))[0]
				Segment_check = (np.array(self.ydimPos[i]) > LowerBoundaryYDir).any() and (np.array(self.ydimPos[i]) < UpperBoundaryYDir).any()
			elif not MaskXdir and not MaskYdir and MaskZdir:
				Indices = np.where(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)))[0]
				Segment_check = (np.array(self.zdimPos[i]) > LowerBoundaryZDir).any() and (np.array(self.zdimPos[i]) < UpperBoundaryZDir).any()
			elif MaskXdir and MaskYdir and not MaskZdir:
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)),\
						np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir))))[0]
				Segment_check = (np.array(self.xdimPos[i]) > LowerBoundaryXDir).any() and (np.array(self.xdimPos[i]) < UpperBoundaryXDir).any() \
								and (np.array(self.ydimPos[i]) > LowerBoundaryYDir).any() and (np.array(self.ydimPos[i]) < UpperBoundaryYDir).any()
			elif MaskXdir and not MaskYdir and MaskZdir:
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)),\
						np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir))))[0]
				Segment_check = (np.array(self.xdimPos[i]) > LowerBoundaryXDir).any() and (np.array(self.xdimPos[i]) < UpperBoundaryXDir).any() \
								and (np.array(self.zdimPos[i]) > LowerBoundaryZDir).any() and (np.array(self.zdimPos[i]) < UpperBoundaryZDir).any()
			elif not MaskXdir and MaskYdir and MaskZdir:
				Indices = np.where(np.logical_and(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)),\
						np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir))))[0]
				Segment_check = (np.array(self.zdimPos[i]) > LowerBoundaryZDir).any() and (np.array(self.zdimPos[i]) < UpperBoundaryZDir).any() \
								and (np.array(self.ydimPos[i]) > LowerBoundaryYDir).any() and (np.array(self.ydimPos[i]) < UpperBoundaryYDir).any()
			elif MaskXdir and MaskYdir and MaskZdir:
				Indices = np.where(
					np.logical_and(np.logical_and(np.greater(self.ydimPos[i], LowerBoundaryYDir), np.less(self.ydimPos[i], UpperBoundaryYDir)),
					np.logical_and(np.logical_and(np.greater(self.zdimPos[i], LowerBoundaryZDir), np.less(self.zdimPos[i], UpperBoundaryZDir)),\
					np.logical_and(np.greater(self.xdimPos[i], LowerBoundaryXDir), np.less(self.xdimPos[i], UpperBoundaryXDir)))))[0]
				Segment_check = (np.array(self.xdimPos[i]) > LowerBoundaryXDir).any() and (np.array(self.xdimPos[i]) < UpperBoundaryXDir).any() \
								and (np.array(self.ydimPos[i]) > LowerBoundaryYDir).any() and (np.array(self.ydimPos[i]) < UpperBoundaryYDir).any() \
								and (np.array(self.zdimPos[i]) > LowerBoundaryZDir).any() and (np.array(self.zdimPos[i]) < UpperBoundaryZDir).any()
			if Segment_check:
				Add_filament(i)
			if not len(Indices) == 0:
				Cutoff_filaments(Indices, i)

		MaskedFilamentSegments = np.array(self.MaskedFilamentSegments)
		MaskedLengths = np.array(self.MaskedLengths)
		zdimMasked = np.array(self.zdimMasked)
		CutOffFilamentSegments = np.array(self.CutOffFilamentSegments)
		CutOffLengths = np.array(self.CutOffLengths)
		CutOffzDim = np.array(self.CutOffzDim)
		print 'Masking time: ', time.clock() - time_start, 's'
		return MaskedFilamentSegments, MaskedLengths, zdimMasked, CutOffFilamentSegments, CutOffLengths, CutOffzDim
