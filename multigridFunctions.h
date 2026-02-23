// Velocity set. On grid interfaces, read and write only the distributions which will get streamed into the opposing grid.
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

void applyCoarseFineGridCommunication( GridStruct &GridCoarse, GridStruct &GridFine )
{
	auto fArrayViewCoarse = GridCoarse.fArray.getView();
	auto shifterViewCoarse = GridCoarse.shifter.getConstView();
	InfoStruct InfoCoarse = GridCoarse.Info;
	auto fArrayViewFine  = GridFine.fArray.getView();
	auto shifterViewFine  = GridFine.shifter.getConstView();
	InfoStruct InfoFine = GridFine.Info;
	
	const int positiveCxDistributions[9] = { 1, 7, 9, 12, 16, 20, 22, 23, 26 };
	const int negativeCxDistributions[9] = { 2, 8, 10, 11, 15, 19, 21, 24, 25 };
	const int positiveCyDistributions[9] = { 6, 12, 13, 15, 17, 19, 22, 24, 26 };
	const int negativeCyDistributions[9] = { 5, 11, 14, 16, 18, 20, 21, 23, 25 };
	const int positiveCzDistributions[9] = { 4, 8, 9, 14, 17, 20, 21, 24, 26 };
	const int negativeCzDistributions[9] = { 3, 7, 10, 13, 18, 19, 22, 23, 25 };

	
	int outerNormalX = 0;
	int outerNormalY = 0;
	int outerNormalZ = 0;
	IntTripleType start; 
	IntTripleType end;
	
	// FIRST, READ FROM COARSE GRID, WRITE TO FINE GRID
	
	auto writeToFineCellLambda = [=] __cuda_callable__ ( const IntTripleType& tripleIndex ) mutable
	{
		const int iFine = tripleIndex[0];
		const int jFine = tripleIndex[1];
		const int kFine = tripleIndex[2];
		int cellFine;
		getCellIndex( cellFine, iFine, jFine, kFine, InfoFine );
		int shiftedIndexFine[27];
		getShiftedIndex( cellFine, shiftedIndexFine, shifterViewFine, InfoFine );
		
		MarkerStruct Marker;
		getMarkers( iFine, jFine, kFine, Marker, InfoFine );
		
		if ( Marker.ghost )
		{
			int direction[9];
			for ( int i = 0; i < 9; i++ )
			{
				if ( outerNormalX == -1 ) direction[i] = positiveCxDistributions[i];
				else if ( outerNormalX == 1 ) direction[i] = negativeCxDistributions[i];
				else if ( outerNormalY == -1 ) direction[i] = positiveCyDistributions[i];
				else if ( outerNormalY == 1 ) direction[i] = negativeCyDistributions[i];
				else if ( outerNormalZ == -1 ) direction[i] = positiveCzDistributions[i];
				else if ( outerNormalZ == 1 ) direction[i] = negativeCzDistributions[i];
			}
			
			const int iCoarse = iFine / 2 + InfoCoarse.iSubgridStart;
			const int jCoarse = jFine / 2 + InfoCoarse.jSubgridStart;
			const int kCoarse = kFine / 2 + InfoCoarse.kSubgridStart;
			int cellCoarse;
			getCellIndex( cellCoarse, iCoarse, jCoarse, kCoarse, InfoCoarse );
			int shiftedIndexCoarse[27];
			getShiftedIndex( cellCoarse, shiftedIndexCoarse, shifterViewCoarse, InfoCoarse );
			
			float f[27];
			for (int i = 0; i < 9; i++) f[direction[i]] = fArrayViewCoarse( direction[i], shiftedIndexCoarse[direction[i]] );	
			for (int i = 0; i < 9; i++) fArrayViewFine( direction[i], shiftedIndexFine[direction[i]] ) = f[direction[i]];	
		}
	};
	
	// Start X
	outerNormalX = -1;
	start = IntTripleType{ 0, 0, 0 };
	end = IntTripleType{ 2, InfoFine.cellCountY, InfoFine.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End X
	outerNormalX = 1;
	start = IntTripleType{ InfoFine.cellCountX-2, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX, InfoFine.cellCountY, InfoFine.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// Start Y
	outerNormalX = 0;
	outerNormalY = -1;
	start = IntTripleType{ 0, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX, 2, InfoFine.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End Y
	outerNormalY = 1;
	start = IntTripleType{ 0, InfoFine.cellCountY-2, 0 };
	end = IntTripleType{ InfoFine.cellCountX, InfoFine.cellCountY, InfoFine.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// Start Z
	outerNormalX = 0;
	outerNormalZ = -1;
	start = IntTripleType{ 0, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX, InfoFine.cellCountY, 2 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End Z
	outerNormalX = 1;
	start = IntTripleType{ 0, 0, InfoFine.cellCountZ-2 };
	end = IntTripleType{ InfoFine.cellCountX, InfoFine.cellCountY, InfoFine.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	outerNormalZ = 0;
	// SECOND, READ FROM FINE GRID, TAKE AVERAGE, WRITE TO COARSE GRID

	auto writeToCoarseCellLambda = [=] __cuda_callable__ ( const IntTripleType& tripleIndex ) mutable
	{
		const int iCoarse = tripleIndex[0];
		const int jCoarse = tripleIndex[1];
		const int kCoarse = tripleIndex[2];
		int cellCoarse;
		getCellIndex( cellCoarse, iCoarse, jCoarse, kCoarse, InfoCoarse );
		int shiftedIndexCoarse[27];
		getShiftedIndex( cellCoarse, shiftedIndexCoarse, shifterViewCoarse, InfoCoarse );
		
		const int iFineFirst = (iCoarse - InfoCoarse.iSubgridStart) * 2;
		const int jFineFirst = (jCoarse - InfoCoarse.jSubgridStart) * 2;
		const int kFineFirst = (kCoarse - InfoCoarse.kSubgridStart) * 2;
		
		float f[27] = {0};
		
		for ( int kAdd = 0; kAdd <= 1; kAdd++ )
		{
			for ( int jAdd = 0; jAdd <= 1; jAdd++ )
			{
				for ( int iAdd = 0; iAdd <= 1; iAdd++ )
				{
					const int iFine = iFineFirst + iAdd;
					const int jFine = jFineFirst + jAdd;
					const int kFine = kFineFirst + kAdd;
					int cellFine;
					getCellIndex( cellFine, iFine, jFine, kFine, InfoFine );
					int shiftedIndexFine[27];
					getShiftedIndex( cellFine, shiftedIndexFine, shifterViewFine, InfoFine );
					for (int direction = 0; direction < 27; direction++) f[direction] += fArrayViewFine( direction, shiftedIndexFine[direction] );	
				}
			}
		}
		for (int direction = 0; direction < 27; direction++) fArrayViewCoarse( direction, shiftedIndexCoarse[direction] ) = f[direction] * 0.125f;
	};
	
	// Start X
	start = IntTripleType{ InfoCoarse.iSubgridStart+1, InfoCoarse.jSubgridStart+1, InfoCoarse.kSubgridStart+1 };
	end = IntTripleType{ InfoCoarse.iSubgridStart+2, InfoCoarse.jSubgridEnd-1, InfoCoarse.kSubgridEnd-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
	// End X
	start = IntTripleType{ InfoCoarse.iSubgridEnd-2, InfoCoarse.jSubgridStart+1, InfoCoarse.kSubgridStart+1 };
	end = IntTripleType{ InfoCoarse.iSubgridEnd-1, InfoCoarse.jSubgridEnd-1, InfoCoarse.kSubgridEnd-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
	// Start Y
	start = IntTripleType{ InfoCoarse.iSubgridStart+1, InfoCoarse.jSubgridStart+1, InfoCoarse.kSubgridStart+1 };
	end = IntTripleType{ InfoCoarse.iSubgridEnd-1, InfoCoarse.jSubgridStart+2, InfoCoarse.kSubgridEnd-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
	// End-1 Y
	start = IntTripleType{ InfoCoarse.iSubgridStart+1, InfoCoarse.jSubgridEnd-2, InfoCoarse.kSubgridStart+1 };
	end = IntTripleType{ InfoCoarse.iSubgridEnd-1, InfoCoarse.jSubgridEnd-1, InfoCoarse.kSubgridEnd-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
	// Start+1 Z
	start = IntTripleType{ InfoCoarse.iSubgridStart+1, InfoCoarse.jSubgridStart+1, InfoCoarse.kSubgridStart+1 };
	end = IntTripleType{ InfoCoarse.iSubgridEnd-1, InfoCoarse.jSubgridEnd-1, InfoCoarse.kSubgridStart+2 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
	// End-1 Z
	start = IntTripleType{ InfoCoarse.iSubgridStart+1, InfoCoarse.jSubgridStart+1, InfoCoarse.kSubgridEnd-2 };
	end = IntTripleType{ InfoCoarse.iSubgridEnd-1, InfoCoarse.jSubgridEnd-1, InfoCoarse.kSubgridEnd-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToCoarseCellLambda );
}
