// Reads known distributions
// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation
// Example: cell has outer normal [1, -1, 0]
// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f22 = f22ArrayView[shiftedIndex[22]];
	float f23 = f23ArrayView[shiftedIndex[23]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f22 = f22ArrayView[shiftedIndex[22]];
	float f24 = f24ArrayView[shiftedIndex[24]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f21 = f21ArrayView[shiftedIndex[21]];
	float f24 = f24ArrayView[shiftedIndex[24]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f21 = f21ArrayView[shiftedIndex[21]];
	float f24 = f24ArrayView[shiftedIndex[24]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f21 = f21ArrayView[shiftedIndex[21]];
	float f23 = f23ArrayView[shiftedIndex[23]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f22 = f22ArrayView[shiftedIndex[22]];
	float f23 = f23ArrayView[shiftedIndex[23]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f24 = f24ArrayView[shiftedIndex[24]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f22 = f22ArrayView[shiftedIndex[22]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f23 = f23ArrayView[shiftedIndex[23]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f21 = f21ArrayView[shiftedIndex[21]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f22 = f22ArrayView[shiftedIndex[22]];
}
elif (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f22 = f22ArrayView[shiftedIndex[22]];
	float f23 = f23ArrayView[shiftedIndex[23]];
}
elif (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f23 = f23ArrayView[shiftedIndex[23]];
}
elif (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f20 = f20ArrayView[shiftedIndex[20]];
	float f21 = f21ArrayView[shiftedIndex[21]];
}
elif (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f21 = f21ArrayView[shiftedIndex[21]];
	float f24 = f24ArrayView[shiftedIndex[24]];
}
elif (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f19 = f19ArrayView[shiftedIndex[19]];
	float f24 = f24ArrayView[shiftedIndex[24]];
}
elif (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f26 = f26ArrayView[shiftedIndex[26]];
}
elif (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f25 = f25ArrayView[shiftedIndex[25]];
}
elif (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f17 = f17ArrayView[shiftedIndex[17]];
	float f24 = f24ArrayView[shiftedIndex[24]];
}
elif (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f9 = f9ArrayView[shiftedIndex[9]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f20 = f20ArrayView[shiftedIndex[20]];
}
elif (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f12 = f12ArrayView[shiftedIndex[12]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f22 = f22ArrayView[shiftedIndex[22]];
}
elif (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f1 = f1ArrayView[shiftedIndex[1]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f7 = f7ArrayView[shiftedIndex[7]];
	float f16 = f16ArrayView[shiftedIndex[16]];
	float f18 = f18ArrayView[shiftedIndex[18]];
	float f23 = f23ArrayView[shiftedIndex[23]];
}
elif (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f3 = f3ArrayView[shiftedIndex[3]];
	float f6 = f6ArrayView[shiftedIndex[6]];
	float f10 = f10ArrayView[shiftedIndex[10]];
	float f13 = f13ArrayView[shiftedIndex[13]];
	float f15 = f15ArrayView[shiftedIndex[15]];
	float f19 = f19ArrayView[shiftedIndex[19]];
}
elif (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	float f0 = f0ArrayView[shiftedIndex[0]];
	float f2 = f2ArrayView[shiftedIndex[2]];
	float f4 = f4ArrayView[shiftedIndex[4]];
	float f5 = f5ArrayView[shiftedIndex[5]];
	float f8 = f8ArrayView[shiftedIndex[8]];
	float f11 = f11ArrayView[shiftedIndex[11]];
	float f14 = f14ArrayView[shiftedIndex[14]];
	float f21 = f21ArrayView[shiftedIndex[21]];
}