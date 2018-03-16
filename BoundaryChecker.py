import numpy as np
import time

class BoundaryChecker():
	"""
	A class that is used to ensure that the filaments passes the boundary with periodic boundary conditions.
	This class and algorithm assumes that the box is a square box, where the upper and lower boundaries for all directions are the same.
	"""
	def __init__(self, xpoints, ypoints, zpoints, FilID, NumFilaments):
		self.BoxSize = 256.0
		self.xdimPos = xpoints
		self.ydimPos = ypoints
		self.zdimPos = zpoints
		self.FilID = FilID
		self.NFils = NumFilaments

	def Get_periodic_boundary_old(self):
		"""
		This function checks whether a filament crosses the boundary or not. 
		If a filament crosses the boundary, the filament will be split into two/three filaments, i.e different lists/arrays.
		More details of the algorithm, see paper.
		Function also computes the length of the filament.
		"""
		time_start = time.clock()
		print 'Checking boundaries'
		BoxSize = self.BoxSize
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
		self.LengthSplitFilament = []
		self.FilamentIDs = []

		def New_points(P1, P2, boundary_dir, boundary):
			""" Computes the points of the two non-boundary crossing coordinates at the boundary. """
			if boundary:
				DifferenceAdd = -BoxSize
				BoxBoundary = self.LowerBound
			else:
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
			elif SplitFilament == 2:
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
					if xBoundary == 1 and yBoundary != 1 and zBoundary != 1:
						#print 'only x'
						Boundary_check = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_check)
						Xpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						xBoundary = 2
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
					elif yBoundary == 1 and xBoundary != 1 and zBoundary != 1:
						#print 'only y'
						Boundary_check = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', Boundary_check)
						Ypoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_check else self.LowerBound
						yBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint, Zpoint)
					elif zBoundary == 1 and xBoundary != 1 and yBoundary != 1:
						#print 'only z'
						Boundary_check = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_check)
						Zpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						zBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint, Ypoint, Zpoint1, Zpoint2)
					elif xBoundary == 1 and yBoundary == 1 and zBoundary != 1:
						print 'x and y', i
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_checkX)
						Xpoint, Zpoint123 = New_points(Point1, Point2, 'y', Boundary_checkY)
						print Zpoint, Zpoint123
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint, Zpoint)
						xBoundary = 2
						yBoundary = 2
					elif xBoundary != 1 and yBoundary == 1 and zBoundary == 1:
						print 'z and y', i
						Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_checkZ)
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', Boundary_checkY)
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint1, Zpoint2)
						zBoundary = 2
						yBoundary = 2
					elif xBoundary == 1 and yBoundary != 1 and zBoundary == 1:
						print 'x and z', i
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_checkX)
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_checkZ)
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
						xBoundary = 2
						zBoundary = 2
					elif xBoundary == 1 and yBoundary == 1 and zBoundary == 1:
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Bounadry_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint, Zpoint)
						xBoundary = 2
						yBoundary = 2
						zBoundary = 2
					xyNewTemp.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp.append(self.xdimPos[i][j])
					yNewTemp.append(self.ydimPos[i][j])
					zNewTemp.append(self.zdimPos[i][j])

				elif SplitFilament == 2:
					if xBoundary == 1 or yBoundary == 1 or zBoundary == 1:
						Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
						Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
					if xBoundary == 1 and yBoundary != 1 and zBoundary != 1:
						Boundary_check = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_check)
						Xpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						xBoundary = 2
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
					elif yBoundary == 1 and xBoundary != 1 and zBoundary != 1:
						Boundary_check = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', Boundary_check)
						Ypoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_check else self.LowerBound
						yBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint, Zpoint)
					elif zBoundary == 1 and xBoundary != 1 and yBoundary != 1:
						Boundary_check = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_check)
						Zpoint1 = self.LowerBound if Boundary_check else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_check else self.LowerBound
						zBoundary = 2
						Append_points(Xpoint, Xpoint, Ypoint, Ypoint, Zpoint1, Zpoint2)
					elif xBoundary == 1 and yBoundary == 1 and zBoundary != 1:
						print 'x and y', i
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_checkX)
						Xpoint, Zpoint123 = New_points(Point1, Point2, 'y', Boundary_checkY)
						print Zpoint, Zpoint123
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint, Zpoint)
						xBoundary = 2
						yBoundary = 2
					elif xBoundary != 1 and yBoundary == 1 and zBoundary == 1:
						print 'z and y', i
						Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_checkZ)
						Xpoint, Zpoint = New_points(Point1, Point2, 'y', Boundary_checkY)
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Append_points(Xpoint, Xpoint, Ypoint1, Ypoint2, Zpoint1, Zpoint2)
						zBoundary = 2
						yBoundary = 2
					elif xBoundary == 1 and yBoundary != 1 and zBoundary == 1:
						print 'x and z', i
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Ypoint, Zpoint = New_points(Point1, Point2, 'x', Boundary_checkX)
						Xpoint, Ypoint = New_points(Point1, Point2, 'z', Boundary_checkZ)
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint, Ypoint, Zpoint, Zpoint)
						xBoundary = 2
						zBoundary = 2
					elif xBoundary == 1 and yBoundary == 1 and zBoundary == 1:
						Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Bounadry_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*BoxSize
						Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
						Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
						Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
						Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
						Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
						Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
						Append_points(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint, Zpoint)
						xBoundary = 2
						yBoundary = 2
						zBoundary = 2
					xyNewTemp2.append(np.asarray([self.xdimPos[i][j], self.ydimPos[i][j]]))
					xNewTemp2.append(self.xdimPos[i][j])
					yNewTemp2.append(self.ydimPos[i][j])
					zNewTemp2.append(self.zdimPos[i][j])

				# Check if boundary has been crossed
				if xDiff > 0.5*BoxSize:
					if SplitFilament == 0 and not (yBoundary or zBoundary):
						SplitFilament = 1
					elif SplitFilament == 1 and not (yBoundary or zBoundary):
						SplitFilament = 2
					xBoundary = 1
				if yDiff > 0.5*BoxSize:
					if SplitFilament == 0 and not (xBoundary or zBoundary):
						SplitFilament = 1
					elif SplitFilament == 1 and not (xBoundary or zBoundary):
						SplitFilament = 2
					yBoundary = 1
				if zDiff > 0.5*BoxSize:
					if SplitFilament == 0 and not (xBoundary or yBoundary):
						SplitFilament = 1
					elif SplitFilament == 1 and not (xBoundary or yBoundary):
						SplitFilament = 2
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
			if i == 424:
				print TempLength
				print xTemp, xNewTemp, xNewTemp2
				print yTemp, yNewTemp, yNewTemp2
				print zTemp, zNewTemp, zNewTemp2
				print 
				print self.xdimPos[i]
				print self.ydimPos[i]
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
	
	def New_points(self, P1, P2, boundary_dir, boundary):
		""" 
		Computes the points of the two non-boundary crossing coordinates at the boundary. 
		See notes for more info.
		If filament crosses two boundaries at once, it is assumed that the t variable gives the same value for the third point.
		"""
		if boundary:
			DifferenceAdd = -self.BoxSize
			BoxBoundary = self.LowerBound
		else:
			DifferenceAdd = self.BoxSize
			BoxBoundary = self.UpperBound
		
		if boundary_dir == 'x':
			t_variable = (BoxBoundary - P1[0])/(P2[0] - P1[0] - self.BoxSize)
			y_coord = P1[1] + (P2[1] - P1[1])*t_variable
			z_coord = P1[2] + (P2[2] - P1[2])*t_variable
			return y_coord, z_coord
		elif boundary_dir == 'y':
			t_variable = (BoxBoundary - P1[1])/(P2[1] - P1[1] - self.BoxSize)
			x_coord = P1[0] + (P2[0] - P1[0])*t_variable
			z_coord = P1[2] + (P2[2] - P1[2])*t_variable
			return x_coord, z_coord
		elif boundary_dir == 'z':
			t_variable = (BoxBoundary - P1[2])/(P2[2] - P1[2] - self.BoxSize)
			x_coord = P1[0] + (P2[0] - P1[0])*t_variable
			y_coord = P1[1] + (P2[1] - P1[1])*t_variable
			return x_coord, y_coord
		else:
			raise ValueError('boundary_dir argument not set properly')

	def check_bc(self, diff):
		Check1 = np.where(diff <= -self.BoxSize/2.0)[0]
		Check2 = np.where(diff >= self.BoxSize/2.0)[0]
		if len(Check1):
			diff[Check1] += self.BoxSize
		if len(Check2):
			diff[Check2] -= self.BoxSize
		return diff

	def compute_lengths(self, xpos, ypos, zpos):
		""" Computes the length of the whole filament. This takes into account for periodic boundaries. """
		length = 0
		diffx = xpos[1:] - xpos[:-1]
		diffy = ypos[1:] - ypos[:-1]
		diffz = zpos[1:] - zpos[:-1]
		
		#diffx = self.check_bc(diffx)
		#diffy = self.check_bc(diffy)
		#diffz = self.scheck_bc(diffz)
		diffx[diffx <= -self.BoxSize/2.0] += self.BoxSize
		diffx[diffx >= self.BoxSize/2.0] -= self.BoxSize
		diffy[diffy <= -self.BoxSize/2.0] += self.BoxSize
		diffy[diffy >= self.BoxSize/2.0] -= self.BoxSize
		diffz[diffz <= -self.BoxSize/2.0] += self.BoxSize
		diffz[diffz >= self.BoxSize/2.0] -= self.BoxSize
		for j in range(len(diffx)):
			length += np.sqrt(diffx[j]**2 + diffy[j]**2 + diffz[j]**2)

		return length

	def Get_periodic_boundary(self):
		"""
		This function checks whether a filament crosses the boundary or not. 
		If a filament crosses the boundary, the filament will be split into two/three filaments, i.e different lists/arrays.
		More details of the algorithm, see paper.
		Function also computes the length of the filament.
		"""
		print 'Checking boundaries'
		time_start = time.clock()
		FilPosTemp = []
		xPosTemp = []
		yPosTemp = []
		zPosTemp = []
		self.FilLengths = []
		self.LengthSplitFilament = []
		self.FilamentIDs = []

		def get_boundary_point(i, j, xbound, ybound, zbound):
			""" Interpolates boundary points for the non-boundary crossing axes. """
			Point1 = np.array([self.xdimPos[i][j-1], self.ydimPos[i][j-1], self.zdimPos[i][j-1]])
			Point2 = np.array([self.xdimPos[i][j], self.ydimPos[i][j], self.zdimPos[i][j]])
			if xbound and not (ybound and zbound):
				Boundary_check = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Ypoint1, Zpoint1 = self.New_points(Point1, Point2, 'x', Boundary_check)
				Ypoint2 = Ypoint1
				Zpoint2 = Zpoint1
				Xpoint1 = self.LowerBound if Boundary_check else self.UpperBound
				Xpoint2 = self.UpperBound if Boundary_check else self.LowerBound
			elif ybound and not (xbound and zbound):
				Boundary_check = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Xpoint1, Zpoint1 = self.New_points(Point1, Point2, 'y', Boundary_check)
				Xpoint2 = Xpoint1
				Zpoint2 = Zpoint1
				Ypoint1 = self.LowerBound if Boundary_check else self.UpperBound
				Ypoint2 = self.UpperBound if Boundary_check else self.LowerBound
			elif zbound and not (xbound and ybound):
				Boundary_check = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Xpoint1, Ypoint1 = self.New_points(Point1, Point2, 'z', Boundary_check)
				Xpoint2 = Xpoint1
				Ypoint2 = Ypoint1
				Zpoint1 = self.LowerBound if Boundary_check else self.UpperBound
				Zpoint2 = self.UpperBound if Boundary_check else self.LowerBound
			elif (xbound and ybound) and not zbound:
				Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Ypoint, Zpoint1 = self.New_points(Point1, Point2, 'x', Boundary_checkX)
				Zpoint2 = Zpoint1
				Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
				Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
				Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
				Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
			elif (xbound and zbound) and not ybound:
				Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Ypoint1, Zpoint = self.New_points(Point1, Point2, 'x', Boundary_checkX)
				Ypoint2 = Ypoint1
				Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
				Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
				Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
				Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
			elif (ybound and zbound) and not xbound:
				Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Xpoint1, Zpoint = self.New_points(Point1, Point2, 'y', Boundary_checkY)
				Xpoint2 = Xpoint1
				Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
				Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
				Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
				Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
			elif xbound and ybound and zbound:
				Boundary_checkX = np.abs(self.xdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Boundary_checkY = np.abs(self.ydimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Boundary_checkZ = np.abs(self.zdimPos[i][j-1] - self.LowerBound) < 0.5*self.BoxSize
				Xpoint1 = self.LowerBound if Boundary_checkX else self.UpperBound
				Xpoint2 = self.UpperBound if Boundary_checkX else self.LowerBound
				Ypoint1 = self.LowerBound if Boundary_checkY else self.UpperBound
				Ypoint2 = self.UpperBound if Boundary_checkY else self.LowerBound
				Zpoint1 = self.LowerBound if Boundary_checkZ else self.UpperBound
				Zpoint2 = self.UpperBound if Boundary_checkZ else self.LowerBound
			else:
				raise ValueError('Something went wrong with setting xbound, ybound or zbound!')
			add_boundary_point(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint1, Zpoint2, i, j)

		def add_boundary_point(Xpoint1, Xpoint2, Ypoint1, Ypoint2, Zpoint1, Zpoint2, i, j):
			""" 
			Adds the points to respective array for the first boundary point.
			Clears the temporary lists containing the coordinate points and adds the second boundary points.
			"""
			# Add the new interpolated points to the old arrays
			self.xyTemp.append(np.array([Xpoint1, Ypoint1]))
			self.xTemp.append(Xpoint1)
			self.yTemp.append(Ypoint1)
			self.zTemp.append(Zpoint1)

			FilPosTemp.append(self.xyTemp)
			xPosTemp.append(self.xTemp)
			yPosTemp.append(self.yTemp)
			zPosTemp.append(self.zTemp)
			self.FilamentIDs.append(self.FilID[i])
			# Create new temporary arrays
			self.xyTemp = []
			self.xTemp = []
			self.yTemp = []
			self.zTemp = []
			
			# Add new extrapolated point to the new arrays
			self.xyTemp.append(np.array([Xpoint2, Ypoint2]))
			self.xTemp.append(Xpoint2)
			self.yTemp.append(Ypoint2)
			self.zTemp.append(Zpoint2)
			# Add current point to the new arrays
			self.xyTemp.append(np.array([self.xdimPos[i][j], self.ydimPos[i][j]]))
			self.xTemp.append(self.xdimPos[i][j])
			self.yTemp.append(self.ydimPos[i][j])
			self.zTemp.append(self.zdimPos[i][j])

		for i in range(self.NFils):
			xBoundary = 0
			yBoundary = 0
			zBoundary = 0
			diffx = self.xdimPos[i][1:] - self.xdimPos[i][:-1]
			diffy = self.ydimPos[i][1:] - self.ydimPos[i][:-1]
			diffz = self.zdimPos[i][1:] - self.zdimPos[i][:-1]
			self.xyTemp = []
			self.xTemp = []
			self.yTemp = []
			self.zTemp = []
			SplitFilament = 0
			SplitCount = 1
			Filament_length = 0
			for j in range(len(diffx)):
				if SplitFilament:
					get_boundary_point(i, j, xBoundary, yBoundary, zBoundary)
					SplitCount += 1
					SplitFilament = 0
					xBoundary = 0
					yBoundary = 0
					zBoundary = 0
				else:
					self.xyTemp.append(np.array([self.xdimPos[i][j], self.ydimPos[i][j]]))
					self.xTemp.append(self.xdimPos[i][j])
					self.yTemp.append(self.ydimPos[i][j])
					self.zTemp.append(self.zdimPos[i][j])

				if np.abs(diffx[j]) > 0.5*self.BoxSize:
					SplitFilament = 1
					xBoundary = 1
				if np.abs(diffy[j]) > 0.5*self.BoxSize:
					SplitFilament = 1
					yBoundary = 1
				if np.abs(diffz[j]) > 0.5*self.BoxSize:
					SplitFilament = 1
					zBoundary = 1
					
			# Add final points
			if xBoundary or yBoundary or zBoundary:
				get_boundary_point(i, j+1, xBoundary, yBoundary, zBoundary)
				SplitCount += 1   # May not be correct here. Make sure to check the lengths of arrays!
			else:
				self.xyTemp.append(np.array([self.xdimPos[i][-1], self.ydimPos[i][-1]]))
				self.xTemp.append(self.xdimPos[i][-1])
				self.yTemp.append(self.ydimPos[i][-1])
				self.zTemp.append(self.zdimPos[i][-1])
			FilPosTemp.append(np.array(self.xyTemp))
			xPosTemp.append(self.xTemp)
			yPosTemp.append(self.yTemp)
			zPosTemp.append(self.zTemp)
			self.FilamentIDs.append(self.FilID[i])
			# Compute length of the whole filament
			Filament_length = self.compute_lengths(self.xdimPos[i], self.ydimPos[i], self.zdimPos[i])
			self.FilLengths.append(Filament_length)

			# Adds filament length to a different array.
			# Corresponds to the colorbar of the filament. Ensures that split filaments retains their total length in the plots.
			for k in range(SplitCount):
				self.LengthSplitFilament.append(Filament_length)

		FilamentIDs = np.array(self.FilamentIDs)
		FilamentPos = np.array(FilPosTemp)
		xdimPosNew = np.array(xPosTemp)
		ydimPosNew = np.array(yPosTemp)
		zdimPosNew = np.array(zPosTemp)
		LengthSplitFilament = np.array(self.LengthSplitFilament)
		FilLengths = np.array(self.FilLengths)
		print 'Boundary check time:', time.clock() - time_start, 's'
		return FilamentIDs, FilamentPos, xdimPosNew, ydimPosNew, zdimPosNew, LengthSplitFilament, FilLengths