// Periodic BC, can also work as zero gradient depending how source cell is chosen
// Reads unknown distributions from source cell
// Writes them into local cell
// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation
// Example: cell has outer normal [1, -1, 0]
// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction
// In this case it is only used to identify unknown distributions and might not match the geometrical meaning

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading unknown distributions from source cell
	const float f2 = f2ArrayView[shiftedSourceIndex[2]];
	const float f8 = f8ArrayView[shiftedSourceIndex[8]];
	const float f10 = f10ArrayView[shiftedSourceIndex[10]];
	const float f11 = f11ArrayView[shiftedSourceIndex[11]];
	const float f15 = f15ArrayView[shiftedSourceIndex[15]];
	const float f19 = f19ArrayView[shiftedSourceIndex[19]];
	const float f21 = f21ArrayView[shiftedSourceIndex[21]];
	const float f24 = f24ArrayView[shiftedSourceIndex[24]];
	const float f25 = f25ArrayView[shiftedSourceIndex[25]];
	// Writing them into local cell
	f2ArrayView[shiftedIndex[2]] = f2;
	f8ArrayView[shiftedIndex[8]] = f8;
	f10ArrayView[shiftedIndex[10]] = f10;
	f11ArrayView[shiftedIndex[11]] = f11;
	f15ArrayView[shiftedIndex[15]] = f15;
	f19ArrayView[shiftedIndex[19]] = f19;
	f21ArrayView[shiftedIndex[21]] = f21;
	f24ArrayView[shiftedIndex[24]] = f24;
	f25ArrayView[shiftedIndex[25]] = f25;
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Reading unknown distributions from source cell
	const float f5 = f5ArrayView[shiftedSourceIndex[5]];
	const float f11 = f11ArrayView[shiftedSourceIndex[11]];
	const float f14 = f14ArrayView[shiftedSourceIndex[14]];
	const float f16 = f16ArrayView[shiftedSourceIndex[16]];
	const float f18 = f18ArrayView[shiftedSourceIndex[18]];
	const float f20 = f20ArrayView[shiftedSourceIndex[20]];
	const float f21 = f21ArrayView[shiftedSourceIndex[21]];
	const float f23 = f23ArrayView[shiftedSourceIndex[23]];
	const float f25 = f25ArrayView[shiftedSourceIndex[25]];
	// Writing them into local cell
	f5ArrayView[shiftedIndex[5]] = f5;
	f11ArrayView[shiftedIndex[11]] = f11;
	f14ArrayView[shiftedIndex[14]] = f14;
	f16ArrayView[shiftedIndex[16]] = f16;
	f18ArrayView[shiftedIndex[18]] = f18;
	f20ArrayView[shiftedIndex[20]] = f20;
	f21ArrayView[shiftedIndex[21]] = f21;
	f23ArrayView[shiftedIndex[23]] = f23;
	f25ArrayView[shiftedIndex[25]] = f25;
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Reading unknown distributions from source cell
	const float f3 = f3ArrayView[shiftedSourceIndex[3]];
	const float f7 = f7ArrayView[shiftedSourceIndex[7]];
	const float f10 = f10ArrayView[shiftedSourceIndex[10]];
	const float f13 = f13ArrayView[shiftedSourceIndex[13]];
	const float f18 = f18ArrayView[shiftedSourceIndex[18]];
	const float f19 = f19ArrayView[shiftedSourceIndex[19]];
	const float f22 = f22ArrayView[shiftedSourceIndex[22]];
	const float f23 = f23ArrayView[shiftedSourceIndex[23]];
	const float f25 = f25ArrayView[shiftedSourceIndex[25]];
	// Writing them into local cell
	f3ArrayView[shiftedIndex[3]] = f3;
	f7ArrayView[shiftedIndex[7]] = f7;
	f10ArrayView[shiftedIndex[10]] = f10;
	f13ArrayView[shiftedIndex[13]] = f13;
	f18ArrayView[shiftedIndex[18]] = f18;
	f19ArrayView[shiftedIndex[19]] = f19;
	f22ArrayView[shiftedIndex[22]] = f22;
	f23ArrayView[shiftedIndex[23]] = f23;
	f25ArrayView[shiftedIndex[25]] = f25;
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
	// Reading unknown distributions from source cell
	const float f1 = f1ArrayView[shiftedSourceIndex[1]];
	const float f7 = f7ArrayView[shiftedSourceIndex[7]];
	const float f9 = f9ArrayView[shiftedSourceIndex[9]];
	const float f12 = f12ArrayView[shiftedSourceIndex[12]];
	const float f16 = f16ArrayView[shiftedSourceIndex[16]];
	const float f20 = f20ArrayView[shiftedSourceIndex[20]];
	const float f22 = f22ArrayView[shiftedSourceIndex[22]];
	const float f23 = f23ArrayView[shiftedSourceIndex[23]];
	const float f26 = f26ArrayView[shiftedSourceIndex[26]];
	// Writing them into local cell
	f1ArrayView[shiftedIndex[1]] = f1;
	f7ArrayView[shiftedIndex[7]] = f7;
	f9ArrayView[shiftedIndex[9]] = f9;
	f12ArrayView[shiftedIndex[12]] = f12;
	f16ArrayView[shiftedIndex[16]] = f16;
	f20ArrayView[shiftedIndex[20]] = f20;
	f22ArrayView[shiftedIndex[22]] = f22;
	f23ArrayView[shiftedIndex[23]] = f23;
	f26ArrayView[shiftedIndex[26]] = f26;
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Reading unknown distributions from source cell
	const float f6 = f6ArrayView[shiftedSourceIndex[6]];
	const float f12 = f12ArrayView[shiftedSourceIndex[12]];
	const float f13 = f13ArrayView[shiftedSourceIndex[13]];
	const float f15 = f15ArrayView[shiftedSourceIndex[15]];
	const float f17 = f17ArrayView[shiftedSourceIndex[17]];
	const float f19 = f19ArrayView[shiftedSourceIndex[19]];
	const float f22 = f22ArrayView[shiftedSourceIndex[22]];
	const float f24 = f24ArrayView[shiftedSourceIndex[24]];
	const float f26 = f26ArrayView[shiftedSourceIndex[26]];
	// Writing them into local cell
	f6ArrayView[shiftedIndex[6]] = f6;
	f12ArrayView[shiftedIndex[12]] = f12;
	f13ArrayView[shiftedIndex[13]] = f13;
	f15ArrayView[shiftedIndex[15]] = f15;
	f17ArrayView[shiftedIndex[17]] = f17;
	f19ArrayView[shiftedIndex[19]] = f19;
	f22ArrayView[shiftedIndex[22]] = f22;
	f24ArrayView[shiftedIndex[24]] = f24;
	f26ArrayView[shiftedIndex[26]] = f26;
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Reading unknown distributions from source cell
	const float f4 = f4ArrayView[shiftedSourceIndex[4]];
	const float f8 = f8ArrayView[shiftedSourceIndex[8]];
	const float f9 = f9ArrayView[shiftedSourceIndex[9]];
	const float f14 = f14ArrayView[shiftedSourceIndex[14]];
	const float f17 = f17ArrayView[shiftedSourceIndex[17]];
	const float f20 = f20ArrayView[shiftedSourceIndex[20]];
	const float f21 = f21ArrayView[shiftedSourceIndex[21]];
	const float f24 = f24ArrayView[shiftedSourceIndex[24]];
	const float f26 = f26ArrayView[shiftedSourceIndex[26]];
	// Writing them into local cell
	f4ArrayView[shiftedIndex[4]] = f4;
	f8ArrayView[shiftedIndex[8]] = f8;
	f9ArrayView[shiftedIndex[9]] = f9;
	f14ArrayView[shiftedIndex[14]] = f14;
	f17ArrayView[shiftedIndex[17]] = f17;
	f20ArrayView[shiftedIndex[20]] = f20;
	f21ArrayView[shiftedIndex[21]] = f21;
	f24ArrayView[shiftedIndex[24]] = f24;
	f26ArrayView[shiftedIndex[26]] = f26;
}