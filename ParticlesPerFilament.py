import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def PartPerFilament(filaments, filamentID, particlepos):
	fil = filaments[0]
	xmin = np.min(fil[:,0])
	xmax = np.max(fil[:,0])
	ymin = np.min(fil[:,1])
	ymax = np.max(fil[:,1])
	zmin = np.min(fil[:,2])
	zmax = np.max(fil[:,2])

	FilamentPos = np.column_stack((fil[:,0], fil[:,1]))
	#x = fil[:,0]
	#y = fil[:,1]
	#FilamentPos = []
	#for i in range(len(x)):
	#	FilamentPos.append(np.array([x[i], y[i]]))
	#print FilamentPos
	fig, ax = plt.subplots()
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	line = LineCollection([FilamentPos], linestyle='solid')
	ax.add_collection(line)
	plt.show()