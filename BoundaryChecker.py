import numpy as np
import time

class BoundaryChecker():
	"""
	A class that is used to ensure that the filaments passes the boundary with periodic boundary conditions.
	This class and algorithm assumes that the box is a square box, where the upper and lower boundaries for all directions are the same.
	"""
	def __init__(self, lower_boundary, upper_boundary, xpoints, ypoints, zpoints, FilID, NumFilaments):
		self.UpperBound = upper_boundary
		self.LowerBound = lower_boundary
		self.xdimPos = xpoints
		self.ydimPos = ypoints
		self.zdimPos = zpoints
		self.FilID = FilID
		self.NFils = NumFilaments

	def Get_periodic_boundary(self):
		"""
		This function checks whether a filament crosses the boundary or not. 
		If a filament crosses the boundary, the filament will be split into two/three filaments, i.e different lists/arrays.
		More details of the algorithm, see paper.
		Function also computes the length of the filament.
		"""
		time_start = time.clock()
		print 'Checking boundaries'
		BoxSize = self.UpperBound - self.LowerBound
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
		self.LengthSplitFilament = []
		self.FilamentIDs = []

		def New_points(P1, P2, boundary_dir, boundary):
			""" Computes the points of the two non-boundary crossing coordinates at the boundary. """
			if boundary == 0:
				DifferenceAdd = -BoxSize
				BoxBoundary = self.LowerBound
			elif boundary == 1:
				DifferenceAdd = BoxSize
				BoxBoundary = self.UpperBound

			if boundary_dir == 'x':
				t_variable = (BoxBoundary - P1[0])/(P2[0] - P1[0] - BoxSize)
				y_coord = P1[1] + (P2[1] - P1[1])*t_variable
				z_coord = P1[2] + (P2[2] - P1[2])*t_variable
				return y_coord, z_coord
			elif boundary_dir == 'y':
				t_variable = (BoxBoundary - P1[1])/(P2[1] - P1[1] - BoxSize)
				x_coord = P1[0] + (P2[0] - P1[0])*t_variable
				z_coord = P1[2] + (P2[2] - P1[2])*t_variable
				return x_coord, z_coord
			elif boundary_dir == 'z':
				t_variable = (BoxBoundary - P1[2])/(P2[2] - P1[2] - BoxSize)
				x_coord = P1[0] + (P2[0] - P1[0])*t_variable
				y_coord = P1[1] + (P2[1] - P1[1])*t_variable
				return x_coord, y_coord
			else:
				raise ValueError('boundary_dir argument not set properly')

		def Append_points(x1, x2, y1, y2, z1, z2):
			if SplitFilament == 1:
				xyTemp.append(np.array([x1, y1]))
				xTemp.append(x1)
				yTemp.append(y1)
				zTemp.append(z1)
				xyNewTemp.append(np.array([x2,y2]))
				xNewTemp.append(x2)
				yNewTemp.append(y2)
				zNewTemp.append(z2)
			if SplitFilament == 2:
				xyNewTemp.append(np.array([x1,y1]))
				xNewTemp.append(x1)
				yNewTemp.append(y1)
				zNewTemp.append(z1)
				xyNewTemp2.append(np.array([x2,y2]))
				xNewTemp2.append(x2)
				yNewTemp2.append(y2)
				zNewTemp2.append(z2)

		for i in range(self.NFils-1):
			SplitFilament = 0
			xBoundary = 0
			yBoundary = 0
			zBoundary = 0
			xyTemp = []
			xTemp = []
			yTemp = []
			zTemp = []
			xyNewTemp = []
			xNewTemp = []
			yNewTemp = []
			zNewTemp = []
			if SplitFilament == 2:
				xyNewTemp2 = []
				xNewTemp2 = []
				yNewTemp2 = []
				zNewTemp2 = []
			for j in range(len(self.xdimPos[i])-1):
				xDiff = np.abs(self.xdimPos[i][j+1] - self.xdimPos[i][j])
				yDiff = np.abs(self.ydimPos[i][j+1] - self.ydimPos[i][j])
				zDiff = np.abs(self.zdimPos[i][j+1] - self.zdimPos[i][j])

				if SplitFilament == 0:
					xyTemp.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xTemp.append(self.xdimPos[i][j])
					yTemp.append(self.ydimPos[i][j])
					zTemp.append(self.zdimPos[i][j])

				elif SplitFilament == 1:
					if xBoundary == 1 or yBoundary == 1 or zBoundary == 1:
						Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
						Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
					if xBoundary == 1:
						Boundary_check = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', 0)
						Xpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						xBoundary = 2
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
					elif yBoundary == 1:
						Boundary_check = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', 0)
						Ypoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_check else self.LowerBound
						yBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint, Zpoint)
					elif zBoundary == 1:
						Boundary_check = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', 0)
						Zpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						zBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint, Ypoint, Zpoint1, Zpoint2)
					xyNewTemp.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp.append(self.xdimPos[i][j])
					yNewTemp.append(self.ydimPos[i][j])
					zNewTemp.append(self.zdimPos[i][j])

				elif SplitFilament == 2:
					if xBoundary == 1 or yBoundary == 1 or zBoundary == 1:
						Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
						Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
					if xBoundary == 1:
						Boundary_check = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', 0)
						Xpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						xBoundary = 2
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
					elif yBoundary == 1:
						Boundary_check = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', 0)
						Ypoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_check else self.LowerBound
						yBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint, Zpoint)
					elif zBoundary == 1:
						Boundary_check = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', 0)
						Zpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						zBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint, Ypoint, Zpoint1, Zpoint2)
					xyNewTemp2.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp2.append(self.xdimPos[i][j])
					yNewTemp2.append(self.ydimPos[i][j])
					zNewTemp2.append(self.zdimPos[i][j])

				# Check if boundary has been crossed
				if xDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					xBoundary = 1
				if yDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					yBoundary = 1
				if zDiff > 0.5*BoxSize:
					if SplitFilament == 1:
						SplitFilament = 2
					SplitFilament += 1
					zBoundary = 1

			# Adds final points to the positions
			if SplitFilament == 0:
				xyTemp.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xTemp.append(self.xdimPos[i][-1])
				yTemp.append(self.ydimPos[i][-1])
				zTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 1:
				xyNewTemp.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xNewTemp.append(self.xdimPos[i][-1])
				yNewTemp.append(self.ydimPos[i][-1])
				zNewTemp.append(self.zdimPos[i][-1])
			elif SplitFilament == 2:
				xyNewTemp2.append(np.asarray([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				xNewTemp2.append(self.xdimPos[i][-1])
				yNewTemp2.append(self.ydimPos[i][-1])
				zNewTemp2.append(self.zdimPos[i][-1])

			# Adds positions to temporary arrays. Also adds ID to array
			if len(zTemp) > 1:
				self.FilamentIDs.append(self.FilID[i])
				FilPosTemp.append(xyTemp)
				xPosTemp.append(xTemp)
				yPosTemp.append(yTemp)
				zPosTemp.append(zTemp)
			if SplitFilament == 1:
				if len(zNewTemp) > 1:
					self.FilamentIDs.append(self.FilID[i])
					FilPosTemp.append(xyNewTemp)
					xPosTemp.append(xNewTemp)
					yPosTemp.append(yNewTemp)
					zPosTemp.append(zNewTemp)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					self.FilamentIDs.append(self.FilID[i])
					FilPosTemp.append(xyNewTemp2)
					xPosTemp.append(xNewTemp2)
					yPosTemp.append(yNewTemp2)
					zPosTemp.append(zNewTemp2)

			# Compute the length of the whole filament
			TempLength = 0
			if len(zTemp) > 1:
				for k in range(len(zTemp)-1):
					TempLength += np.sqrt((xyTemp[k+1][0]-xyTemp[k][0])**2 + (xyTemp[k+1][1] - xyTemp[k][1])**2 + (zTemp[k+1]-zTemp[k])**2)
			if len(zNewTemp) > 1:
				for k in range(len(zNewTemp)-1):
					TempLength += np.sqrt((xyNewTemp[k+1][0]-xyNewTemp[k][0])**2 + (xyNewTemp[k+1][1] - xyNewTemp[k][1])**2 \
										+ (zNewTemp[k+1]-zNewTemp[k])**2)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					for k in range(len(zNewTemp2)-1):
						TempLength += np.sqrt((xyNewTemp2[k+1][0]-xyNewTemp2[k][0])**2 + (xyNewTemp2[k+1][1] - xyNewTemp2[k][1])**2 \
											+ (zNewTemp2[k+1]-zNewTemp2[k])**2)
			self.FilLengths.append(TempLength)
			# Seperating length to arrays to give filaments that is split the same total length
			if len(zTemp) > 1:
				self.LengthSplitFilament.append(TempLength)
			if len(zNewTemp) > 1: 
				self.LengthSplitFilament.append(TempLength)
			if SplitFilament == 2:
				if len(zNewTemp2) > 1:
					self.LengthSplitFilament.append(TempLength)
				
		FilamentIDs = np.array(self.FilamentIDs)
		FilamentPos = np.array(FilPosTemp)
		xdimPosNew = np.array(xPosTemp)
		ydimPosNew = np.array(yPosTemp)
		zdimPosNew = np.array(zPosTemp)
		LengthSplitFilament = np.array(self.LengthSplitFilament)
		FilLengths = np.asarray(self.FilLengths)
		print 'Boundary check time:', time.clock() - time_start, 's'
		return FilamentIDs, FilamentPos, xdimPosNew, ydimPosNew, zdimPosNew, LengthSplitFilament, FilLengths
	