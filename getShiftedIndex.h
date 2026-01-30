// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

__device__ __constant__ int cxArray[27] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
__device__ __constant__ int cyArray[27] = { 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
__device__ __constant__ int czArray[27] = { 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

__cuda_callable__ void getShiftedIndex( 	const int &iCell, const int &jCell, const int &kCell, 
											int (&shiftedIndex)[27], InfoStruct &Info )
// Performs periodic shift streaming by shifting the pointer index
{
	const int cell = kCell * Info.cellCountX * Info.cellCountY + jCell * Info.cellCountX + iCell;	
	int modX = Info.iterationsFinished % Info.cellCountX;
	int modY = Info.iterationsFinished % Info.cellCountY;
	int modZ = Info.iterationsFinished % Info.cellCountZ;
	for ( int direction = 0; direction < 27; direction++ )
	{	
		int shiftX = modX * ( - cxArray[direction] ) + Info.cellCountX;
		int shiftY = modY * ( - cyArray[direction] ) + Info.cellCountY;
		int shiftZ = modZ * ( - czArray[direction] ) + Info.cellCountZ;
		const int shift = shiftZ * Info.cellCountX * Info.cellCountY + shiftY * Info.cellCountX + shiftX;
		shiftedIndex[direction] = (cell + shift) % (Info.cellCountX * Info.cellCountY * Info.cellCountZ);
	}
}
