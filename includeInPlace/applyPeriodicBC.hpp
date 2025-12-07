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
	ArrayView[shiftedIndex[2]] = f2;
	ArrayView[shiftedIndex[8]] = f8;
	ArrayView[shiftedIndex[10]] = f10;
	ArrayView[shiftedIndex[11]] = f11;
	ArrayView[shiftedIndex[15]] = f15;
	ArrayView[shiftedIndex[19]] = f19;
	ArrayView[shiftedIndex[21]] = f21;
	ArrayView[shiftedIndex[24]] = f24;
	ArrayView[shiftedIndex[25]] = f25;
}
elif (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
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
	ArrayView[shiftedIndex[5]] = f5;
	ArrayView[shiftedIndex[11]] = f11;
	ArrayView[shiftedIndex[14]] = f14;
	ArrayView[shiftedIndex[16]] = f16;
	ArrayView[shiftedIndex[18]] = f18;
	ArrayView[shiftedIndex[20]] = f20;
	ArrayView[shiftedIndex[21]] = f21;
	ArrayView[shiftedIndex[23]] = f23;
	ArrayView[shiftedIndex[25]] = f25;
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
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
	ArrayView[shiftedIndex[3]] = f3;
	ArrayView[shiftedIndex[7]] = f7;
	ArrayView[shiftedIndex[10]] = f10;
	ArrayView[shiftedIndex[13]] = f13;
	ArrayView[shiftedIndex[18]] = f18;
	ArrayView[shiftedIndex[19]] = f19;
	ArrayView[shiftedIndex[22]] = f22;
	ArrayView[shiftedIndex[23]] = f23;
	ArrayView[shiftedIndex[25]] = f25;
}
elif (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
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
	ArrayView[shiftedIndex[1]] = f1;
	ArrayView[shiftedIndex[7]] = f7;
	ArrayView[shiftedIndex[9]] = f9;
	ArrayView[shiftedIndex[12]] = f12;
	ArrayView[shiftedIndex[16]] = f16;
	ArrayView[shiftedIndex[20]] = f20;
	ArrayView[shiftedIndex[22]] = f22;
	ArrayView[shiftedIndex[23]] = f23;
	ArrayView[shiftedIndex[26]] = f26;
}
elif (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
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
	ArrayView[shiftedIndex[6]] = f6;
	ArrayView[shiftedIndex[12]] = f12;
	ArrayView[shiftedIndex[13]] = f13;
	ArrayView[shiftedIndex[15]] = f15;
	ArrayView[shiftedIndex[17]] = f17;
	ArrayView[shiftedIndex[19]] = f19;
	ArrayView[shiftedIndex[22]] = f22;
	ArrayView[shiftedIndex[24]] = f24;
	ArrayView[shiftedIndex[26]] = f26;
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
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
	ArrayView[shiftedIndex[4]] = f4;
	ArrayView[shiftedIndex[8]] = f8;
	ArrayView[shiftedIndex[9]] = f9;
	ArrayView[shiftedIndex[14]] = f14;
	ArrayView[shiftedIndex[17]] = f17;
	ArrayView[shiftedIndex[20]] = f20;
	ArrayView[shiftedIndex[21]] = f21;
	ArrayView[shiftedIndex[24]] = f24;
	ArrayView[shiftedIndex[26]] = f26;
}