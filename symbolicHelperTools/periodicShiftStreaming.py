import copy

nx = 4
ny = 3

velocities = [
	[0, 0],
	[1, 0],
	[0, 1],
	[-1, 0],
	[0, -1],
	[1, 1],
	[-1, 1],
	[-1, -1],
	[1, -1]	
	]

alphabet = "abcdefghijklmnopqrstuvwxyz"

cell_names = []
cell_values = []
print("cell names")
for i in range(ny):
	cell_names.append([])
	cell_values.append([])
	for j in range(nx):
		index = i * nx + j
		cell_names[i].append(alphabet[index])
		cell_values[i].append(index)
	print(cell_names[-1])

print("memory layout")
flat_cell_names = [x for sub in cell_names for x in sub]
print(flat_cell_names)
print("values before streaming for every direction")
flat_cell_values = [x for sub in cell_values for x in sub]
print(flat_cell_values)
print()

for velocity in velocities:
	print("result of streaming with velocity", velocity)
	streamed_cell_values = copy.deepcopy(cell_values)
	for i in range(ny):
		for j in range(nx):
			source_i = i + velocity[1]
			source_j = j - velocity[0]
			if source_i < 0 or source_i >= ny:
				streamed_cell_values[i][j] = "?"
				continue
			if source_j < 0 or source_j >= nx:
				streamed_cell_values[i][j] = "?"
				continue
			streamed_cell_values[i][j] = cell_values[source_i][source_j]
	flat_streamed_cell_values = [x for sub in streamed_cell_values for x in sub]
	print(flat_streamed_cell_values)
