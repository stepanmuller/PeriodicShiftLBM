import sys
import numpy as np
from matplotlib import pyplot as plt
import copy
from stl import mesh
	
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
	n = 5 # higher = flatter, goes to full thickness more quickly
	return (1 - np.abs(2 * s - 1)**n)**(1/n)

def getProfile( parameters ):
	
	h, a, b, l, w, s, thickness = parameters
	controlPoints = np.array([ [0, 0], [a, b], [l, h] ])
	
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
		
	return profilePointsTop, profilePointsBottom
	
def visualizeProfile(curve_top, curve_bottom):
	plt.scatter(curve_top[:, 0], curve_top[:, 1], color="black")
	plt.scatter(curve_bottom[:, 0], curve_bottom[:, 1], color="gray")
	ax = plt.gca()
	ax.set_aspect("equal")
	plt.show()

if __name__ == "__main__":
	inputCount = len(sys.argv)
	if inputCount < 7:
		h = 10
		a = 10
		b = 3
		l = 15
		w = 2
		s = 15
		t = 1
	else:
		nothing, h, a, b, l, w, s, t = sys.argv
	parameters = [h, a, b, l, w, s, t]
	profilePointsTop, profilePointsBottom = getProfile( parameters )
	if inputCount < 7:
		visualizeProfile( profilePointsTop, profilePointsBottom )
		
	width = 1 * h
	triangleList = []
	# Main surface top
	for i in range(len(profilePointsTop)-1):
		zy0 = profilePointsTop[i]
		zy1 = profilePointsTop[i+1]
		pt0 = [0, zy0[1], zy0[0]]
		pt1 = [0, zy1[1], zy1[0]]
		pt2 = [width, zy0[1], zy0[0]]
		pt3 = [width, zy1[1], zy1[0]]
		triangleList.append([pt0, pt1, pt2])
		triangleList.append([pt1, pt2, pt3])
			
	for j in range(len(surfaceMatrix)-1):
		i = 0
		pt1 = surfaceMatrix[j, i]
		pt2 = surfaceMatrix[j+1, i]
		pt3 = np.array([0, -27-b/2, -80])
		triangleList.append([pt1, pt2, pt3])

	for j in range(len(surfaceMatrix)-1):
		i = len(surfaceMatrix[j])-1
		pt1 = surfaceMatrix[j, i]
		pt2 = surfaceMatrix[j+1, i]
		pt3 = pt1.copy()
		pt3[-1] = 100
		pt4 = pt2.copy()
		pt4[-1] = 100
		pt5 = np.array([0, 0, 100])
		triangleList.append([pt1, pt2, pt3])
		triangleList.append([pt2, pt3, pt4])
		triangleList.append([pt3, pt4, pt5])
				
	triangles = np.array(triangleList)

	trianglesMirror = triangles.copy()
	trianglesMirror[:, :, 0] = - trianglesMirror[:, :, 0]

	triangles = np.concatenate([triangles, trianglesMirror], axis=0)

	surface_mesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
	surface_mesh.vectors = triangles
	filename = "intake.STL"
	surface_mesh.save(filename)
	print(f"Surface saved as {filename}")
	

