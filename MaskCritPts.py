import numpy as np

def Mask_CPs(xdimPos, ydimPos, zdimPos, Masklist, MaskCheck):
	""" Function that masks the critical points """
	print 'Masking critical points'
	UpperBoundaryXDir, UpperBoundaryYDir, UpperBoundaryZDir, LowerBoundaryXDir, LowerBoundaryYDir, LowerBoundaryZDir = Masklist
	MaskXdir, MaskYdir, MaskZdir = MaskCheck

	if MaskXdir and not MaskYdir and not MaskZdir:
		Indices = np.where(np.logical_and(np.greater(xdimPos, LowerBoundaryXDir), np.less(xdimPos, UpperBoundaryXDir)))[0]
	elif not MaskXdir and MaskYdir and not MaskZdir:
		Indices = np.where(np.logical_and(np.greater(ydimPos, LowerBoundaryYDir), np.less(ydimPos, UpperBoundaryYDir)))[0]
	elif not MaskXdir and not MaskYdir and MaskZdir:
		Indices = np.where(np.logical_and(np.greater(zdimPos, LowerBoundaryZDir), np.less(zdimPos, UpperBoundaryZDir)))[0]
	elif MaskXdir and MaskYdir and not MaskZdir:
		Indices = np.where(np.logical_and(np.logical_and(np.greater(xdimPos, LowerBoundaryXDir), np.less(xdimPos, UpperBoundaryXDir)),\
				np.logical_and(np.greater(ydimPos, LowerBoundaryYDir), np.less(ydimPos, UpperBoundaryYDir))))[0]
	elif MaskXdir and not MaskYdir and MaskZdir:
		Indices = np.where(np.logical_and(np.logical_and(np.greater(xdimPos, LowerBoundaryXDir), np.less(xdimPos, UpperBoundaryXDir)),\
				np.logical_and(np.greater(zdimPos, LowerBoundaryZDir), np.less(zdimPos, UpperBoundaryZDir))))[0]
	elif not MaskXdir and MaskYdir and MaskZdir:
		Indices = np.where(np.logical_and(np.logical_and(np.greater(zdimPos, LowerBoundaryZDir), np.less(zdimPos, UpperBoundaryZDir)),\
				np.logical_and(np.greater(ydimPos, LowerBoundaryYDir), np.less(ydimPos, UpperBoundaryYDir))))[0]
	elif MaskXdir and MaskYdir and MaskZdir:
		Indices = np.where(
			np.logical_and(np.logical_and(np.greater(ydimPos, LowerBoundaryYDir), np.less(ydimPos, UpperBoundaryYDir)),
			np.logical_and(np.logical_and(np.greater(zdimPos, LowerBoundaryZDir), np.less(zdimPos, UpperBoundaryZDir)),\
			np.logical_and(np.greater(xdimPos, LowerBoundaryXDir), np.less(xdimPos, UpperBoundaryXDir)))))[0]

	return xdimPos[Indices], ydimPos[Indices], zdimPos[Indices]