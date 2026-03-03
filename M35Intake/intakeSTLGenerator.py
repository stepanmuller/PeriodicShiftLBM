import numpy as np
from stl import mesh

def generateIntakeSTL( parameters ):
	c, b, w, By, Ny = parameters
	tResolution = 20
	sResolution = 100
	
	def getBernstein3(t):
		return np.array([
			(1 - t)**3,
			3 * t * (1 - t)**2,
			3 * t**2 * (1 - t),
			t**3 ])

	def getSurfacePoint(controlMatrix, weightMatrix, t, s):
		Bt = getBernstein3(t)
		Bs = getBernstein3(s)
		bernsteinMatrix = np.outer(Bt, Bs)
		weights = weightMatrix * bernsteinMatrix
		weightedPoints = weights[:, :, np.newaxis] * controlMatrix
		point = np.sum(weightedPoints, axis=(0, 1))
		denominator = np.sum(weights)
		point = point / denominator
		return point

	A = np.array([0, -27, -80])
	B = np.array([0, By, 30-c])
	C = np.array([0, 15, 30-c])
	D = np.array([0, 15, 30])
	E = np.array([-w, -27, -80])
	F = np.array([-w, By, 30-c])
	G = np.array([-29.28, 15, 30-c])
	H = np.array([-29.28, 15, 30])
	I = np.array([-w, -27-b, -80])
	J = np.array([-w, Ny, 30-c])
	K = np.array([-29.28, -15, 30-c])
	L = np.array([-29.28, -15, 30])
	M = np.array([0, -27-b, -80])
	N = np.array([0, Ny, 30-c])
	O = np.array([0, -15, 30-c])
	P = np.array([0, -15, 30])

	controlPointMatrix = np.array([
	[A, B, C, D],
	[E, F, G, H],
	[I, J, K, L],
	[M, N, O, P]])

	weightMatrix = np.array([
	[1, 1, 1, 1],
	[3, 3, 0.35, 0.35],
	[3, 3, 0.35, 0.35],
	[1, 1, 1, 1]
	])

	tRange = np.linspace(0, 1, tResolution)
	sRange = np.linspace(0, 1, sResolution)

	surfaceMatrix = np.zeros((tResolution, sResolution, 3))

	for j in range(len(surfaceMatrix)):
		for i in range(len(surfaceMatrix[j])):
			surfaceMatrix[j, i] = getSurfacePoint(controlPointMatrix, weightMatrix, tRange[j], sRange[i])
			if np.abs(surfaceMatrix[j, i, 0]) < 0.0001:
				surfaceMatrix[j, i, 0] = 0.0
			
	triangleList = []
	for j in range(len(surfaceMatrix)-1):
		for i in range(len(surfaceMatrix[j])-1):
			pt1 = surfaceMatrix[j, i]
			pt2 = surfaceMatrix[j+1, i]
			pt3 = surfaceMatrix[j+1, i+1]
			pt4 = surfaceMatrix[j, i+1]
			triangleList.append([pt1, pt2, pt3])
			triangleList.append([pt1, pt3, pt4])
			
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
		pt3[-1] = 60
		pt4 = pt2.copy()
		pt4[-1] = 60
		pt5 = np.array([0, 0, 60])
		triangleList.append([pt1, pt2, pt3])
		triangleList.append([pt1, pt3, pt4])
		triangleList.append([pt3, pt4, pt5])
				
	triangles = np.array(triangleList)

	trianglesMirror = triangles.copy()
	trianglesMirror[:, :, 0] = - trianglesMirror[:, :, 0]

	triangles = np.concatenate([triangles, trianglesMirror], axis=0)

	surface_mesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
	surface_mesh.vectors = triangles
	filename = "intakeChannel.STL"
	surface_mesh.save(filename)
	print(f"Surface saved as {filename}")

c = 40
b = 15
w = 20
By = 0
Ny = -20
parameters = [c, b, w, By, Ny]
generateIntakeSTL(parameters)
