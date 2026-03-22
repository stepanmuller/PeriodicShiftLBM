// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

IntArrayType cxArray{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
IntArrayType cyArray{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
IntArrayType czArray{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

// Perform periodic shift streaming by adjusting F.shifter
void applyStreaming( GridStruct& Grid )
{
	auto shifterView = Grid.shifter.getView();
	InfoStruct Info = Grid.Info;
	auto cxArrayView = cxArray.getConstView();
	auto cyArrayView = cyArray.getConstView();
	auto czArrayView = czArray.getConstView();
	auto directionLambda = [=] __cuda_callable__ ( const int direction ) mutable
	{
		int shift = shifterView[direction];
		shift -= cxArrayView[direction];
		shift -= Info.cellCountX * cyArrayView[direction];
		shift -= Info.cellCountY * Info.cellCountX * czArrayView[direction];
		if ( shift < 0 ) shift += Info.cellCount;
		if ( shift < 0 ) shift += Info.cellCount; // one more time to allow cellCountZ = 1
		if ( shift >= Info.cellCount ) shift -= Info.cellCount;
		if ( shift >= Info.cellCount ) shift -= Info.cellCount; // one more time to allow cellCountZ = 1
		shifterView[direction] = shift;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, 27, directionLambda );
}

__cuda_callable__ void getShiftedIndex( const int& cell, int (&shiftedIndex)[27], IntArrayConstViewType shifterView, const InfoStruct &Info )
{
    for ( int direction = 0; direction < 27; direction++ ) 
		{
			const int shift = shifterView[direction];
			shiftedIndex[direction] = cell + shift;
			if (shiftedIndex[direction] >= Info.cellCount) shiftedIndex[direction] -= Info.cellCount;
		}	
}

// Esotwist version: Just flip the EsotwistFlipper to alter between odd / even iterations
void applyStreaming( DIADGridStruct& Grid )
{
	Grid.esotwistFlipper = !Grid.esotwistFlipper;
}

/*
__cuda_callable__ void getEsotwistReadIndex( const int& cell, int (&cellIndex)[27], int (&fIndex)[27], DIADEsotwistNbrStruct &Nbr, const bool esotwistFlipper, const InfoStruct &Info )
{
	cellIndex = {	cell, cell, Nbr.i, Nbr.k, cell, Nbr.j, cell, 
					Nbr.k, Nbr.i, cell, Nbr.ik, Nbr.ij, cell, Nbr.k, Nbr.j, Nbr.i, Nbr.j, cell, Nbr.jk, 
					Nbr.ik, Nbr.j, Nbr.ij, Nbr.k, Nbr.jk, Nbr.i, Nbr.ijk, cell };
	if ( esotwistFlipper )
	{
		fIndex = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	}
	else
	{
		fIndex = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
	}
}

__cuda_callable__ void getEsotwistWriteIndex( const int& cell, int (&cellIndex)[27], int (&fIndex)[27], DIADEsotwistNbrStruct &Nbr, const bool esotwistFlipper, const InfoStruct &Info )
{
	cellIndex = {	cell, Nbr.i, cell, cell, Nbr.k, cell, Nbr.j, 
					Nbr.i, Nbr.k, Nbr.ik, cell, cell, Nbr.ij, Nbr.j, Nbr.k, Nbr.j, Nbr.i, Nbr.jk, cell, 
					Nbr.j, Nbr.ik, Nbr.k, Nbr.ij, Nbr.i, Nbr.jk, cell, Nbr.ijk };
	if ( esotwistFlipper )
	{
		fIndex = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
		
	}
	else
	{
		fIndex = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	}
}
*/

__cuda_callable__ void getEsotwistReadIndex( const int& cell, int (&cellIndex)[27], int (&fIndex)[27], DIADEsotwistNbrStruct &Nbr, const bool esotwistFlipper, const InfoStruct &Info )
{
    cellIndex[0]  = cell;
    cellIndex[1]  = cell;
    cellIndex[2]  = Nbr.i;
    cellIndex[3]  = Nbr.k;
    cellIndex[4]  = cell;
    cellIndex[5]  = Nbr.j;
    cellIndex[6]  = cell;
    cellIndex[7]  = Nbr.k;
    cellIndex[8]  = Nbr.i;
    cellIndex[9]  = cell;
    cellIndex[10] = Nbr.ik;
    cellIndex[11] = Nbr.ij;
    cellIndex[12] = cell;
    cellIndex[13] = Nbr.k;
    cellIndex[14] = Nbr.j;
    cellIndex[15] = Nbr.i;
    cellIndex[16] = Nbr.j;
    cellIndex[17] = cell;
    cellIndex[18] = Nbr.jk;
    cellIndex[19] = Nbr.ik;
    cellIndex[20] = Nbr.j;
    cellIndex[21] = Nbr.ij;
    cellIndex[22] = Nbr.k;
    cellIndex[23] = Nbr.jk;
    cellIndex[24] = Nbr.i;
    cellIndex[25] = Nbr.ijk;
    cellIndex[26] = cell;

    if ( esotwistFlipper )
    {
        fIndex[0]  = 0;   fIndex[1]  = 1;   fIndex[2]  = 2;
        fIndex[3]  = 3;   fIndex[4]  = 4;   fIndex[5]  = 5;
        fIndex[6]  = 6;   fIndex[7]  = 7;   fIndex[8]  = 8;
        fIndex[9]  = 9;   fIndex[10] = 10;  fIndex[11] = 11;
        fIndex[12] = 12;  fIndex[13] = 13;  fIndex[14] = 14;
        fIndex[15] = 15;  fIndex[16] = 16;  fIndex[17] = 17;
        fIndex[18] = 18;  fIndex[19] = 19;  fIndex[20] = 20;
        fIndex[21] = 21;  fIndex[22] = 22;  fIndex[23] = 23;
        fIndex[24] = 24;  fIndex[25] = 25;  fIndex[26] = 26;
    }
    else
    {
        fIndex[0]  = 0;   fIndex[1]  = 2;   fIndex[2]  = 1;
        fIndex[3]  = 4;   fIndex[4]  = 3;   fIndex[5]  = 6;
        fIndex[6]  = 5;   fIndex[7]  = 8;   fIndex[8]  = 7;
        fIndex[9]  = 10;  fIndex[10] = 9;   fIndex[11] = 12;
        fIndex[12] = 11;  fIndex[13] = 14;  fIndex[14] = 13;
        fIndex[15] = 16;  fIndex[16] = 15;  fIndex[17] = 18;
        fIndex[18] = 17;  fIndex[19] = 20;  fIndex[20] = 19;
        fIndex[21] = 22;  fIndex[22] = 21;  fIndex[23] = 24;
        fIndex[24] = 23;  fIndex[25] = 26;  fIndex[26] = 25;
    }
}

__cuda_callable__ void getEsotwistWriteIndex( const int& cell, int (&cellIndex)[27], int (&fIndex)[27], DIADEsotwistNbrStruct &Nbr, const bool esotwistFlipper, const InfoStruct &Info )
{
    cellIndex[0]  = cell;
    cellIndex[1]  = Nbr.i;
    cellIndex[2]  = cell;
    cellIndex[3]  = cell;
    cellIndex[4]  = Nbr.k;
    cellIndex[5]  = cell;
    cellIndex[6]  = Nbr.j;
    cellIndex[7]  = Nbr.i;
    cellIndex[8]  = Nbr.k;
    cellIndex[9]  = Nbr.ik;
    cellIndex[10] = cell;
    cellIndex[11] = cell;
    cellIndex[12] = Nbr.ij;
    cellIndex[13] = Nbr.j;
    cellIndex[14] = Nbr.k;
    cellIndex[15] = Nbr.j;
    cellIndex[16] = Nbr.i;
    cellIndex[17] = Nbr.jk;
    cellIndex[18] = cell;
    cellIndex[19] = Nbr.j;
    cellIndex[20] = Nbr.ik;
    cellIndex[21] = Nbr.k;
    cellIndex[22] = Nbr.ij;
    cellIndex[23] = Nbr.i;
    cellIndex[24] = Nbr.jk;
    cellIndex[25] = cell;
    cellIndex[26] = Nbr.ijk;

    if ( esotwistFlipper )
    {
        fIndex[0]  = 0;   fIndex[1]  = 2;   fIndex[2]  = 1;
        fIndex[3]  = 4;   fIndex[4]  = 3;   fIndex[5]  = 6;
        fIndex[6]  = 5;   fIndex[7]  = 8;   fIndex[8]  = 7;
        fIndex[9]  = 10;  fIndex[10] = 9;   fIndex[11] = 12;
        fIndex[12] = 11;  fIndex[13] = 14;  fIndex[14] = 13;
        fIndex[15] = 16;  fIndex[16] = 15;  fIndex[17] = 18;
        fIndex[18] = 17;  fIndex[19] = 20;  fIndex[20] = 19;
        fIndex[21] = 22;  fIndex[22] = 21;  fIndex[23] = 24;
        fIndex[24] = 23;  fIndex[25] = 26;  fIndex[26] = 25;
    }
    else
    {
        fIndex[0]  = 0;   fIndex[1]  = 1;   fIndex[2]  = 2;
        fIndex[3]  = 3;   fIndex[4]  = 4;   fIndex[5]  = 5;
        fIndex[6]  = 6;   fIndex[7]  = 7;   fIndex[8]  = 8;
        fIndex[9]  = 9;   fIndex[10] = 10;  fIndex[11] = 11;
        fIndex[12] = 12;  fIndex[13] = 13;  fIndex[14] = 14;
        fIndex[15] = 15;  fIndex[16] = 16;  fIndex[17] = 17;
        fIndex[18] = 18;  fIndex[19] = 19;  fIndex[20] = 20;
        fIndex[21] = 21;  fIndex[22] = 22;  fIndex[23] = 23;
        fIndex[24] = 24;  fIndex[25] = 25;  fIndex[26] = 26;
    }
}
