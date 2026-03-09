import sys
import numpy as np
from matplotlib import pyplot as plt
import copy
from stl import mesh
from STLVisualizer import show_stl
	
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
	
	h, a, b, l, w, thickness, s = parameters
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
	
	# Correction for exact start and end match
	profilePointsTop[0] = profilePointsBottom[0]
	profilePointsTop[-1] = profilePointsBottom[-1]
	
	return profilePointsTop, profilePointsBottom
	
def visualizeProfile(curve_top, curve_bottom):
	plt.scatter(curve_top[:, 0], curve_top[:, 1], color="black")
	plt.scatter(curve_bottom[:, 0], curve_bottom[:, 1], color="gray")
	ax = plt.gca()
	ax.set_aspect("equal")
	plt.show()

def generateVane(parameters):
	profilePointsTop, profilePointsBottom = getProfile( parameters )
	h, a, b, l, w, thickness, s = parameters
	width = h
	triangleList = []
	# Main surface top
	for i in range(len(profilePointsTop)-1):
		xy0 = profilePointsTop[i]
		xy1 = profilePointsTop[i+1]
		pt0 = [xy0[0], xy0[1], 0]
		pt1 = [xy1[0], xy1[1], 0]
		pt2 = [xy0[0], xy0[1], width]
		pt3 = [xy1[0], xy1[1], width]
		triangleList.append([pt0, pt1, pt2])
		triangleList.append([pt1, pt2, pt3])	
	# Main surface bottom
	for i in range(len(profilePointsBottom)-1):
		xy0 = profilePointsBottom[i]
		xy1 = profilePointsBottom[i+1]
		pt0 = [xy0[0], xy0[1], 0]
		pt1 = [xy1[0], xy1[1], 0]
		pt2 = [xy0[0], xy0[1], width]
		pt3 = [xy1[0], xy1[1], width]
		triangleList.append([pt0, pt1, pt2])
		triangleList.append([pt1, pt2, pt3])
	# Closing the side
	for i in range(len(profilePointsTop)-2):
		if i > 0:
			xy0 = profilePointsTop[i]
			xy1 = profilePointsTop[i+1]
			xy2 = profilePointsBottom[i]
			xy3 = profilePointsBottom[i+1]
			pt0 = [xy0[0], xy0[1], width]
			pt1 = [xy1[0], xy1[1], width]
			pt2 = [xy2[0], xy2[1], width]
			pt3 = [xy3[0], xy3[1], width]
			triangleList.append([pt0, pt1, pt2])
			triangleList.append([pt1, pt2, pt3])
	# Final front and rear side triangles
	xy0 = profilePointsTop[0]
	xy1 = profilePointsTop[1]
	xy2 = profilePointsBottom[1]
	pt0 = [xy0[0], xy0[1], width]
	pt1 = [xy1[0], xy1[1], width]
	pt2 = [xy2[0], xy2[1], width]
	triangleList.append([pt0, pt1, pt2])
	xy0 = profilePointsTop[-1]
	xy1 = profilePointsTop[-2]
	xy2 = profilePointsBottom[-2]
	pt0 = [xy0[0], xy0[1], width]
	pt1 = [xy1[0], xy1[1], width]
	pt2 = [xy2[0], xy2[1], width]
	triangleList.append([pt0, pt1, pt2])
				
	triangles = np.array(triangleList)

	trianglesMirror = triangles.copy()
	trianglesMirror[:, :, -1] = - trianglesMirror[:, :, -1]

	triangles = np.concatenate([triangles, trianglesMirror], axis=0)

	surface_mesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
	surface_mesh.vectors = triangles
	filename = "vane.STL"
	surface_mesh.save(filename)
	print(f"Vane saved as {filename}")

if __name__ == "__main__":
	inputCount = len(sys.argv)
	if inputCount < 7:
		h = 10
		a = 10
		b = 3
		l = 15
		w = 2
		t = 1
		s = 15
	else:
		h, a, b, l, w, t, s = map(float, sys.argv[1:8])
	parameters = [h, a, b, l, w, t, s]
	generateVane(parameters)
	if inputCount < 7:
		show_stl( "vane.STL", [0] )
