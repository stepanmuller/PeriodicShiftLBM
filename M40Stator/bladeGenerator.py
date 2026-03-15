import sys
import numpy as np
from matplotlib import pyplot as plt
import copy
from stl import mesh
import pyvista as pv

def getCamberlinePointAndNormal(t, controlPoints, w):
	A = controlPoints[0]
	B = controlPoints[1]
	C = controlPoints[2]
	b0 = (1 - t)**2
	b1 = 2 * (1 - t) * t * w
	b2 = t**2
	b0Derived = -2 * (1 - t)
	b1Derived = 2 * w - 4 * t * w
	b2Derived = 2 * t
	nominator = b0 * A + b1 * B + b2 * C
	denominator = b0 + b1 + b2
	point = nominator / denominator
	nominatorDerived = b0Derived * A + b1Derived * B + b2Derived * C
	denominatorDerived = b0Derived + b1Derived + b2Derived
	tangent = ( nominatorDerived * denominator - nominator * denominatorDerived ) / (denominator**2)
	normal = np.array([-tangent[1], tangent[0]]) / np.linalg.norm(tangent)
	return point, normal
	
def getThicknessMultiplier(s):
	s = s * 0.8 # so that s goes from 0 to 0.8, resulting in non zero thickness at the trailing edge
	thickness = 0.594689181*(0.298222773*s**0.5 - 0.127125232*s - 0.357907906*s**2 + 0.291984971*s**3 - 0.105174606*s**4)
	thickness = thickness / 0.06 # normalizing thickness to 1
	return thickness

def getProfile( parameters, thickness ):
	
	gamma, c, beta, b, w = parameters
	Az = -(c/2) * np.cos(np.radians(gamma))
	As = -(c/2) * np.sin(np.radians(gamma))
	Cz = (c/2) * np.cos(np.radians(gamma))
	Cs = (c/2) * np.sin(np.radians(gamma))
	Bz = Az + b * np.cos(np.radians(gamma + beta))
	Bs = As + b * np.sin(np.radians(gamma + beta))
	controlPoints = np.array([ [Az, As], [Bz, Bs], [Cz, Cs] ])
	
	resolution = 300
	tArray = np.linspace(0, 1, resolution)
	camberlinePoints = np.zeros((resolution, 2))
	camberlineNormal = np.zeros((resolution, 2))
	
	profilePointsTop = np.zeros((resolution, 2))
	profilePointsBottom = np.zeros((resolution, 2))
	
	for iteration in range(10):
		camberlineCurvilinearDistance = np.zeros(resolution)
		for i, t in enumerate(tArray):
			camberlinePoints[i], camberlineNormal[i] = getCamberlinePointAndNormal(t, controlPoints, w)
			if i > 0:
				ds = np.linalg.norm(camberlinePoints[i] - camberlinePoints[i-1])
				camberlineCurvilinearDistance[i] = camberlineCurvilinearDistance[i-1] + ds
		camberlineCurvilinearDistance = camberlineCurvilinearDistance / camberlineCurvilinearDistance[-1]
		
		thicknessArray = np.zeros_like(tArray)
		for i, t in enumerate(tArray):
			thicknessArray[i] = thickness * getThicknessMultiplier( camberlineCurvilinearDistance[i] )
		
		profilePointsTop = camberlinePoints + 0.5 * camberlineNormal * thicknessArray[:, np.newaxis]
		profilePointsBottom = camberlinePoints - 0.5 * camberlineNormal * thicknessArray[:, np.newaxis]
		
		profileCurvilinearDistance = np.zeros(resolution)
		for i, t in enumerate(tArray):
			if i > 0:
				dsTop = np.linalg.norm(profilePointsTop[i] - profilePointsTop[i-1])
				dsBottom = np.linalg.norm(profilePointsBottom[i] - profilePointsBottom[i-1])
				profileCurvilinearDistance[i] = profileCurvilinearDistance[i-1] + (dsTop + dsBottom) / 2
		profileCurvilinearDistance = profileCurvilinearDistance / profileCurvilinearDistance[-1]

		# Now, because weight of the bezier middle point might be high, and thickness is mapped from zero at the start, points are not evenly spaced -> rearrange t
		tArray = np.interp(np.linspace(0, 1, resolution), profileCurvilinearDistance, tArray)
	
	profile = np.concatenate( [profilePointsBottom[::-1, :], profilePointsTop[1:, :]], axis=0 ) #this reverses bottom points and removes duplicate point from the top points
	
	return profile
	
def visualizeProfile( profile ):
	plt.scatter(profile[:, 0], profile[:, 1], color="black")
	ax = plt.gca()
	ax.set_aspect("equal")
	plt.show()

def transformProfileTo3D( profile, R ):
	Uz = profile[:, 0]
	Us = profile[:, 1]
	Ux = R * np.sin( Us / R )
	Uy = R * np.cos( Us / R )
	profile3D = np.array((Ux, Uy, Uz)).T
	return profile3D
	
def visualizeSTL( filename ):
	mesh = pv.read(filename)
	mesh.plot(color="white", show_edges=True)

def generateBlade(parameterMatrix):
	RIn = 7
	RMid = 11
	ROut = 15
	thicknessIn = 2
	thicknessMid = 2
	thicknessOut = 2
	
	profile0 = getProfile(parameterMatrix[0], thicknessIn)
	profile1 = getProfile(parameterMatrix[1], thicknessMid)
	profile2 = getProfile(parameterMatrix[2], thicknessOut)
	
	profile0 = transformProfileTo3D(profile0, RIn)
	profile1 = transformProfileTo3D(profile1, RMid)
	profile2 = transformProfileTo3D(profile2, ROut)
	
	resolution = 100
	qArray = np.linspace(-0.05, 1.05, resolution)
	bladePoints = []
	for q in qArray:
		bladeQRow = (2*q**2-3*q+1) * profile0 + (-4*q**2+4*q) * profile1 + (2*q**2-q) * profile2
		bladePoints.append(bladeQRow)
	bladePoints = np.array(bladePoints)

	triangleList = []
	
	# Main surface
	for q in range(len(bladePoints)-1):
		for t in range(len(bladePoints[0])-1):
			pt0 = bladePoints[q, t]
			pt1 = bladePoints[q, t+1]
			pt2 = bladePoints[q+1, t+1]
			pt3 = bladePoints[q+1, t]
			triangleList.append([pt0, pt1, pt2])
			triangleList.append([pt0, pt2, pt3])
			
	# Closing the trailing edge
	for q in range(len(bladePoints)-1):
		pt0 = bladePoints[q, 0]
		pt1 = bladePoints[q+1, 0]
		pt2 = bladePoints[q+1, -1]
		pt3 = bladePoints[q, -1]
		triangleList.append([pt0, pt1, pt2])
		triangleList.append([pt0, pt2, pt3])
		
	# Closing the hub
	tMid = int(len(bladePoints[0])/2)
	for t in range(tMid-1):
		pt0 = bladePoints[0, t]
		pt1 = bladePoints[0, -1-t]
		pt2 = bladePoints[0, -2-t]
		pt3 = bladePoints[0, t+1]
		triangleList.append([pt0, pt1, pt2])
		triangleList.append([pt0, pt2, pt3])
	pt0 = bladePoints[0, tMid-1]
	pt1 = bladePoints[0, tMid]
	pt2 = bladePoints[0, tMid+1]
	triangleList.append([pt1, pt0, pt2])
	
	# Closing the tip
	tMid = int(len(bladePoints[0])/2)
	for t in range(tMid-1):
		pt0 = bladePoints[-1, t]
		pt1 = bladePoints[-1, -1-t]
		pt2 = bladePoints[-1, -2-t]
		pt3 = bladePoints[-1, t+1]
		triangleList.append([pt1, pt0, pt2])
		triangleList.append([pt2, pt0, pt3])
	pt0 = bladePoints[-1, tMid-1]
	pt1 = bladePoints[-1, tMid]
	pt2 = bladePoints[-1, tMid+1]
	triangleList.append([pt0, pt1, pt2])
				
	triangles = np.array(triangleList)

	surface_mesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
	surface_mesh.vectors = triangles
	surface_mesh.save("blade.STL")

if __name__ == "__main__":
	gamma = [45, 67, 74]
	c = [20, 20, 20]
	beta = [16, 3, 1]
	b = [10, 10, 10]
	w = [1, 1, 1]
	parameterMatrixList = [gamma, c, beta, b, w]
	parameterMatrix = np.array(parameterMatrixList).T
	generateBlade(parameterMatrix)
	visualizeSTL( "blade.STL" )
	
	gamma = [-3.5, -19, -15]
	c = [20, 26.875, 38.75]
	beta = [-3.5, -41.5, -37.5]
	b = [10, 10, 10]
	w = [1, 1, 1]
	parameterMatrixList = [gamma, c, beta, b, w]
	parameterMatrix = np.array(parameterMatrixList).T
	generateBlade(parameterMatrix)
	visualizeSTL( "blade.STL" )
