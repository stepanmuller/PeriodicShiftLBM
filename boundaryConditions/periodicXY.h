// Periodic boundary condition for X and Y. Note that in Z, the periodic boundary condition is automatically applied
// by the periodic shift streaming, to enable it just use the Marker.periodicZ

// To use periodic BC in X and Y, call the relevant function just after streaming and enable Marker.periodicX or Y

// Velocity set
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

void applyBCPeriodicX( GridStruct &Grid )
{
	const InfoStruct Info = Grid.Info;
	auto fArrayView = Grid.fArray.getView();
	auto shifterView = Grid.shifter.getConstView();
	
	const int positiveCxDistributions[9] = { 1, 7, 9, 12, 16, 20, 22, 23, 26 };
	const int negativeCxDistributions[9] = { 2, 8, 10, 11, 15, 19, 21, 24, 25 };
	
	IntPairType start = IntPairType{ 0, 0 };
	IntPairType end = IntPairType{ Info.cellCountY, Info.cellCountZ };
	
	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int jCell = doubleIndex[0];
		const int kCell = doubleIndex[1];
		const int iCellLeft = 0;
		const int iCellRight = Info.cellCountX-1;
		int cellLeft, cellRight;
		getCellIndex( cellLeft, iCellLeft, jCell, kCell, Info );
		getCellIndex( cellRight, iCellRight, jCell, kCell, Info );
		
		int shiftedIndexLeft[27];
		getShiftedIndex( cellLeft, shiftedIndexLeft, shifterView, Info );
		int shiftedIndexRight[27];
		getShiftedIndex( cellRight, shiftedIndexRight, shifterView, Info );
		
		// READ DISTRIBUTIONS WITH POSITIVE CX FROM THE RIGHT CELL, WRITE TO THE LEFT CELL
		for (int i = 0; i < 9; i++) 
		{
			fArrayView( positiveCxDistributions[i], shiftedIndexLeft[positiveCxDistributions[i]] )
			= 
			fArrayView( positiveCxDistributions[i], shiftedIndexRight[positiveCxDistributions[i]] );
		}
		
		// READ DISTRIBUTIONS WITH NEGATIVE CX FROM THE LEFT CELL, WRITE TO THE RIGHT CELL
		for (int i = 0; i < 9; i++) 
		{
			fArrayView( negativeCxDistributions[i], shiftedIndexRight[negativeCxDistributions[i]] )
			= 
			fArrayView( negativeCxDistributions[i], shiftedIndexLeft[negativeCxDistributions[i]] );
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

void applyBCPeriodicY( GridStruct &Grid )
{
	const InfoStruct Info = Grid.Info;
	auto fArrayView = Grid.fArray.getView();
	auto shifterView = Grid.shifter.getConstView();
	
	const int positiveCyDistributions[9] = { 6, 12, 13, 15, 17, 19, 22, 24, 26 };
	const int negativeCyDistributions[9] = { 5, 11, 14, 16, 18, 20, 21, 23, 25 };
	
	IntPairType start = IntPairType{ 0, 0 };
	IntPairType end = IntPairType{ Info.cellCountX, Info.cellCountZ };
	
	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int iCell = doubleIndex[0];
		const int kCell = doubleIndex[1];
		const int jCellLeft = 0;
		const int jCellRight = Info.cellCountY-1;
		int cellLeft, cellRight;
		getCellIndex( cellLeft, iCell, jCellLeft, kCell, Info );
		getCellIndex( cellRight, iCell, jCellRight, kCell, Info );
		
		int shiftedIndexLeft[27];
		getShiftedIndex( cellLeft, shiftedIndexLeft, shifterView, Info );
		int shiftedIndexRight[27];
		getShiftedIndex( cellRight, shiftedIndexRight, shifterView, Info );
		
		// READ DISTRIBUTIONS WITH POSITIVE CY FROM THE RIGHT CELL, WRITE TO THE LEFT CELL
		for (int i = 0; i < 9; i++)
		{
			fArrayView( positiveCyDistributions[i], shiftedIndexLeft[positiveCyDistributions[i]] )
			= 
			fArrayView( positiveCyDistributions[i], shiftedIndexRight[positiveCyDistributions[i]] );
		}
		
		// READ DISTRIBUTIONS WITH NEGATIVE CY FROM THE LEFT CELL, WRITE TO THE RIGHT CELL
		for (int i = 0; i < 9; i++)
		{
			fArrayView( negativeCyDistributions[i], shiftedIndexRight[negativeCyDistributions[i]] )
			= 
			fArrayView( negativeCyDistributions[i], shiftedIndexLeft[negativeCyDistributions[i]] );
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
