void applyCoarseFineGridCommunication( GridStruct &GridCoarse, GridStruct &GridFine )
{
	auto fArrayViewCoarse = GridCoarse.fArray.getView();
	auto shifterViewCoarse = GridCoarse.shifter.getConstView();
	InfoStruct InfoCoarse = GridCoarse.Info;
	auto fArrayViewFine  = GridFine.fArray.getView();
	auto shifterViewFine  = GridFine.shifter.getConstView();
	InfoStruct InfoFine = GridFine.Info;
	
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
		
		const int iCoarse = iFine / 2 + InfoCoarse.iSubgridStart;
		const int jCoarse = jFine / 2 + InfoCoarse.jSubgridStart;
		const int kCoarse = kFine / 2 + InfoCoarse.kSubgridStart;
		int cellCoarse;
		getCellIndex( cellCoarse, iCoarse, jCoarse, kCoarse, InfoCoarse );
		int shiftedIndexCoarse[27];
		getShiftedIndex( cellCoarse, shiftedIndexCoarse, shifterViewCoarse, InfoCoarse );
		
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayViewCoarse( direction, shiftedIndexCoarse[direction] );	
		for (int direction = 0; direction < 27; direction++) fArrayViewFine( direction, shiftedIndexFine[direction] ) = f[direction];	
	};
	
	// Start X
	IntTripleType start{ 0, 0, 0 };
	IntTripleType end{ 2, InfoFine.cellCountY-1, InfoFine.cellCountZ-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End X
	start = IntTripleType{ InfoFine.cellCountX-3, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX-1, InfoFine.cellCountY-1, InfoFine.cellCountZ-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// Start Y
	start = IntTripleType{ 0, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX-1, 2, InfoFine.cellCountZ-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End Y
	start = IntTripleType{ 0, InfoFine.cellCountY-3, 0 };
	end = IntTripleType{ InfoFine.cellCountX-1, InfoFine.cellCountY-1, InfoFine.cellCountZ-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// Start Z
	start = IntTripleType{ 0, 0, 0 };
	end = IntTripleType{ InfoFine.cellCountX-1, InfoFine.cellCountY-1, 2 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	// End Z
	start = IntTripleType{ 0, 0, InfoFine.cellCountZ-3 };
	end = IntTripleType{ InfoFine.cellCountX-1, InfoFine.cellCountY-1, InfoFine.cellCountZ-1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, writeToFineCellLambda );
	
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
