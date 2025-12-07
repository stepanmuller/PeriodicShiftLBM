// Reads known distributions
// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation
// Example: cell has outer normal [1, -1, 0]
// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f22 = f22ArrayView[shiftedIndex[22]];
}
else if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f23 = f23ArrayView[shiftedIndex[23]];
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f23 = f23ArrayView[shiftedIndex[23]];
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f21 = f21ArrayView[shiftedIndex[21]];
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f24 = f24ArrayView[shiftedIndex[24]];
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f24 = f24ArrayView[shiftedIndex[24]];
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f26 = f26ArrayView[shiftedIndex[26]];
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f25 = f25ArrayView[shiftedIndex[25]];
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f24 = f24ArrayView[shiftedIndex[24]];
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f20 = f20ArrayView[shiftedIndex[20]];
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f22 = f22ArrayView[shiftedIndex[22]];
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f23 = f23ArrayView[shiftedIndex[23]];
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f19 = f19ArrayView[shiftedIndex[19]];
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Reading known distributions
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f21 = f21ArrayView[shiftedIndex[21]];
}
