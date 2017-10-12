import numpy as np
import os
import matplotlib.pyplot as plt 

file_directory = 'C:/Users/Alex/Documents/Masters_project/Disperse/lcdm_z0_testing/LCDM_z0_64PeriodicTesting/'
filename_ascii = 'DelaunayOutput_LCDMz064.NDnet.a.NDnet'
filename_binary = 'DelaunayOutput_LCDMz064.NDnet'
UnitConverter = 1

def Read_Delaunay_ascii(file_dir, file_name):
	""" Reading the ascii file, still work in progress """
	datafiles = open(os.path.join(file_dir, file_name), 'r')
	Skip_lines = 0
	NVertex_check = 0
	Vertex_x = []
	Vertex_y = []
	Vertex_z = []
	for line in datafiles:
		if Skip_lines < 4:
			Skip_lines += 1
		else:
			if NVertex_check == 0:
				NumVertices = int(line)
				NVertex_check = 1
			else:
				Vertex = 0

	datafiles.close()

def Read_Delaunay_binary(file_dir, file_name):
	""" Reads binary NDnet file """
	f = open(os.path.join(file_dir, file_name), 'rb')
	# Info stuff, see documentation for more info
	dummy = np.fromfile(f, np.int32, 1)
	tag = np.fromfile(f, 'c', 16)
	dummy = np.fromfile(f, np.int32, 2)
	ndims = np.fromfile(f, np.int32, 1)[0]
	ndims_net = np.fromfile(f, np.int32, 1)[0]
	dummy = np.fromfile(f, np.int32, 2)
	comment = np.fromfile(f, 'c', 80)
	periodicity = np.fromfile(f, np.int32, 1)[0]
	isSimpComplex = np.fromfile(f, np.int32, 1)[0]
	x0 = np.fromfile(f, np.float64, ndims)
	delta = np.fromfile(f, np.float64, ndims)
	index_size = np.fromfile(f, np.int32, 1)[0]
	cumindex_size = np.fromfile(f, np.int32, 1)[0]
	dummy_ext = np.fromfile(f, 'c', 152)
	if index_size == 8:
		NDNET_UINT = np.int64
	elif index_size == 4:
		NDNET_UINT = np.int32
	elif NDNET_UINT == 2:
		NDNET_UINT = np.int16

	if cumindex_size == 8:
		NDNET_IDCUMT = np.int64
	elif cumindex_size == 4:
		NDNET_IDCUMT = np.int32
	elif cumNDNET_UINT == 2:
		NDNET_IDCUMT = np.int16
	nvertex = np.fromfile(f, NDNET_UINT, 1)[0]
	dummy = np.fromfile(f, np.int32, 2)
	v_coords = np.fromfile(f, np.float32, ndims*nvertex)
	Vertex_x = v_coords[0::3]
	Vertex_y = v_coords[1::3]
	Vertex_z = v_coords[2::3]
	dummy = np.fromfile(f, np.int32, 2)
	nfaces = np.fromfile(f, NDNET_UINT, ndims+1)
	dummy = np.fromfile(f, np.int32, 2)
	haveVertexFromFace = np.fromfile(f, np.int32, ndims+1)
	dummy = np.fromfile(f, np.int32, 1)
	f_vertexIndex = []
	for i in range(len(haveVertexFromFace)):
		if haveVertexFromFace[i] == 1:
			dummy = np.fromfile(f, np.int32, 1)
			f_vertexIndex.append(np.fromfile(f, NDNET_UINT, (i+1)*nfaces[i]))
			dummy = np.fromfile(f, np.int32, 1)
	dummy = np.fromfile(f, np.int32, 1)
	haveFaceFromVertex = np.fromfile(f, np.int32, ndims+1)
	dummy = np.fromfile(f, np.int32, 1)
	for i in range(len(haveFaceFromVertex)):
		if haveFaceFromVertex[i] == 1:
			dummy = np.fromfile(f, np.int32, 1)
			numFaceIndexCum = np.fromfile(f, NDNET_IDCUMT, nvertex+1)
			dummy = np.fromfile(f, np.int32, 2)
			v_faceIndex = np.fromfile(f, NDNET_UINT, numFaceIndexCum)
			dummy = np.fromfile(f, np.int32, 1)
	dummy = np.fromfile(f, np.int32, 1)
	haveFaceFromFace = np.fromfile(f, np.int32, (ndims+1)**2)
	dummy = np.fromfile(f, np.int32, 1)
	f.close()

	#plt.hist2d(Vertex_x, Vertex_y, bins=150)
	#plt.show()
	
if __name__ == '__main__':
	#Read_Delaunay_ascii(file_directory, filename_ascii)
	Read_Delaunay_binary(file_directory, filename_binary)