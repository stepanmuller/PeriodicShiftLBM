#include <TNL/Algorithms/sort.h>

// Builds IJK from Info by filling the whole domain	
void buildIJKFromInfo( IJKArrayStruct &IJK, InfoStruct &Info )
{
	IJK.iArray.setSize( Info.cellCount );
	IJK.jArray.setSize( Info.cellCount );
	IJK.kArray.setSize( Info.cellCount );
	auto iView = IJK.iArray.getView();
	auto jView = IJK.jArray.getView();
	auto kView = IJK.kArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
		iView[ cell ] = iCell;
		jView[ cell ] = jCell;
		kView[ cell ] = kCell;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, Info.cellCount, cellLambda );
}

void shiftIJK( IJKArrayStruct &IJK, const int cx, const int cy, const int cz )
{
	auto iView = IJK.iArray.getView();
	auto jView = IJK.jArray.getView();
	auto kView = IJK.kArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		iView[ cell ] = iView[ cell ] + cx;
		jView[ cell ] = jView[ cell ] + cy;
		kView[ cell ] = kView[ cell ] + cz;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, IJK.iArray.getSize(), cellLambda );
}

void shiftIJKPeriodic( IJKArrayStruct &IJK, const int cx, const int cy, const int cz, InfoStruct &Info )
{
	auto iView = IJK.iArray.getView();
	auto jView = IJK.jArray.getView();
	auto kView = IJK.kArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		iView[ cell ] = (iView[ cell ] + cx) % Info.cellCountX;
		jView[ cell ] = (jView[ cell ] + cy) % Info.cellCountY;
		kView[ cell ] = (kView[ cell ] + cz) % Info.cellCountZ;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, IJK.iArray.getSize(), cellLambda );
}

// sort IJK arrays as ascending k, j, i so that k changes the slowest
void sortIJK( IJKArrayStruct &IJK )
{
	const int cellCount = IJK.iArray.getSize();
	auto iView = IJK.iArray.getView();
	auto jView = IJK.jArray.getView();
	auto kView = IJK.kArray.getView();

    auto comparisonLambda = [=] __cuda_callable__ ( const size_t a, const size_t b )
	{
		if ( kView[ a ] < kView[ b ] ) return true;
		else if ( kView[ a ] > kView[ b ] ) return false;
		else
		{
			if ( jView[ a ] < jView[ b ] ) return true;
			else if ( jView[ a ] > jView[ b ] ) return false;
			else
			{
				if ( iView[ a ] < iView[ b ] ) return true;
				else if ( iView[ a ] > iView[ b ] ) return false;
				else return false;
			}
		}
	};
	auto swapLambda = [=] __cuda_callable__ ( const size_t a, const size_t b ) mutable
	{
		TNL::swap( iView[ a ], iView[ b ] );
		TNL::swap( jView[ a ], jView[ b ] );
		TNL::swap( kView[ a ], kView[ b ] );
	};
	TNL::Algorithms::sort<TNL::Devices::Cuda, size_t>( 0, cellCount, comparisonLambda, swapLambda );
}

// final sort for finest grid
void sortFinestGrid( DIADGridStruct &Grid, BoolArrayType &keepCellMarkerArray )
{
	const int cellCount = Grid.Info.cellCount;
	auto iView = Grid.IJK.iArray.getView();
	auto jView = Grid.IJK.jArray.getView();
	auto kView = Grid.IJK.kArray.getView();
	auto bouncebackMarkerView = Grid.bouncebackMarkerArray.getView();
	auto keepCellMarkerView = keepCellMarkerArray.getView();

    auto comparisonLambda = [=] __cuda_callable__ ( const size_t a, const size_t b )
	{
		if ( keepCellMarkerView[ a ] > keepCellMarkerView[ b ] ) return true;
		else if ( keepCellMarkerView[ a ] < keepCellMarkerView[ b ] ) return false;
		else
		{
			if ( kView[ a ] < kView[ b ] ) return true;
			else if ( kView[ a ] > kView[ b ] ) return false;
			else
			{
				if ( jView[ a ] < jView[ b ] ) return true;
				else if ( jView[ a ] > jView[ b ] ) return false;
				else
				{
					if ( iView[ a ] < iView[ b ] ) return true;
					else if ( iView[ a ] > iView[ b ] ) return false;
					else return false;
				}
			}
		}
	};
	auto swapLambda = [=] __cuda_callable__ ( const size_t a, const size_t b ) mutable
	{
		TNL::swap( iView[ a ], iView[ b ] );
		TNL::swap( jView[ a ], jView[ b ] );
		TNL::swap( kView[ a ], kView[ b ] );
		TNL::swap( keepCellMarkerView[ a ], keepCellMarkerView[ b ] );
		TNL::swap( bouncebackMarkerView[ a ], bouncebackMarkerView[ b ] );
	};
	TNL::Algorithms::sort<TNL::Devices::Cuda, size_t>( 0, cellCount, comparisonLambda, swapLambda );
}

// sort for all coarse grid levels
void sortCoarseGrid( DIADGridStruct &Grid, BoolArrayType &keepCellMarkerArray, BoolArrayType &fineToCoarseMarkerArray, BoolArrayType &coarseToFineMarkerArray )
{
	const int cellCount = Grid.Info.cellCount;
	auto iView = Grid.IJK.iArray.getView();
	auto jView = Grid.IJK.jArray.getView();
	auto kView = Grid.IJK.kArray.getView();
	auto keepCellMarkerView = keepCellMarkerArray.getView();
	auto fineToCoarseMarkerView = fineToCoarseMarkerArray.getView();
	auto coarseToFineMarkerView = coarseToFineMarkerArray.getView();

    auto comparisonLambda = [=] __cuda_callable__ ( const size_t a, const size_t b )
	{
		if ( keepCellMarkerView[ a ] > keepCellMarkerView[ b ] ) return true;
		else if ( keepCellMarkerView[ a ] < keepCellMarkerView[ b ] ) return false;
		else
		{
			if ( kView[ a ] < kView[ b ] ) return true;
			else if ( kView[ a ] > kView[ b ] ) return false;
			else
			{
				if ( jView[ a ] < jView[ b ] ) return true;
				else if ( jView[ a ] > jView[ b ] ) return false;
				else
				{
					if ( iView[ a ] < iView[ b ] ) return true;
					else if ( iView[ a ] > iView[ b ] ) return false;
					else return false;
				}
			}
		}
	};
	auto swapLambda = [=] __cuda_callable__ ( const size_t a, const size_t b ) mutable
	{
		TNL::swap( iView[ a ], iView[ b ] );
		TNL::swap( jView[ a ], jView[ b ] );
		TNL::swap( kView[ a ], kView[ b ] );
		TNL::swap( keepCellMarkerView[ a ], keepCellMarkerView[ b ] );
		TNL::swap( fineToCoarseMarkerView[ a ], fineToCoarseMarkerView[ b ] );
		TNL::swap( coarseToFineMarkerView[ a ], coarseToFineMarkerView[ b ] );
	};
	TNL::Algorithms::sort<TNL::Devices::Cuda, size_t>( 0, cellCount, comparisonLambda, swapLambda );
}

int countMarkerCells( BoolArrayType &markerArray )
{
	const int cellCount = markerArray.getSize();
	auto markerView = markerArray.getView();
	
	auto fetch = [ = ] __cuda_callable__( const int cell )
		{
			return (int)markerView[ cell ];
		};
		auto reduction = [] __cuda_callable__( const int& a, const int& b )
		{
			return a + b;
		};	
	int result = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, cellCount, fetch, reduction, 0 );
	return result;
}

int countInvalidIndexes( IntArrayType &intArray )
{
	const int cellCount = intArray.getSize();
	auto intView = intArray.getView();
	
	auto fetch = [ = ] __cuda_callable__( const int cell )
		{
			if (intView[ cell ] < 0) return 1;
			else return 0;
		};
		auto reduction = [] __cuda_callable__( const int& a, const int& b )
		{
			return a + b;
		};	
	int result = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, cellCount, fetch, reduction, 0 );
	return result;
}

// Helper for the findMatchingIJKIndex function
void getIJKStepArray( IJKArrayStruct &Source, IntArrayType &kStepArray )
{
	const int sourceCellCount = Source.iArray.getSize();
	
	IntArrayTypeCPU kSourceArrayCPU;
	kSourceArrayCPU = Source.kArray;
	IntArrayTypeCPU kStepArrayCPU = IntArrayTypeCPU( sourceCellCount );
	int counter = 0;
	int lastIndex = -2;
	for ( int cell = 0; cell < sourceCellCount; cell++ )
	{
		int currentIndex = kSourceArrayCPU[ cell ];
		if ( currentIndex != lastIndex )
		{
			kStepArrayCPU[ counter ] = cell;
			lastIndex = currentIndex;
			counter++;
		}
	}
	const int kStepCount = counter;
	kStepArrayCPU.resize( kStepCount );
	kStepArray = kStepArrayCPU;
}

// Receives sorted IJKSource that !must! be already sorted with key k, j, i so that k index changes the slowest
// For each Wanted, find a matching cell in Source. If its found, write its index to resultArray. If there is no such cell, write -2.		
void findMatchingIJKIndex( 	IJKArrayStruct &Wanted, IJKArrayStruct &Source, IntArrayType &resultArray, IntArrayType &kStepArray )
{
	const int wantedCellCount = Wanted.iArray.getSize();
	const int sourceCellCount = Source.iArray.getSize();
	const int kStepCount = kStepArray.getSize();
	
	auto iWantedView = Wanted.iArray.getConstView();
	auto jWantedView = Wanted.jArray.getConstView();
	auto kWantedView = Wanted.kArray.getConstView();
	auto iSourceView = Source.iArray.getConstView();
	auto jSourceView = Source.jArray.getConstView();
	auto kSourceView = Source.kArray.getConstView();
	auto resultView = resultArray.getView();
	auto kStepView = kStepArray.getConstView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cellWanted ) mutable
	{
		if ( resultView[ cellWanted ] >= 0 ) return; // Do not overwrite the cell if its already valid from before
		const int iWanted = iWantedView[ cellWanted ];
		const int jWanted = jWantedView[ cellWanted ];
		const int kWanted = kWantedView[ cellWanted ];
		int start = 0;
		int end = kStepCount;
		int result = -2;
		// Now by using the kIndexChanges array, we will reduce the search interval
		while ( end > start )
		{
			int half = start + ( end - start ) / 2;
			if ( kSourceView[kStepView[half]] == kWanted )
			{
				result = half;
				break;
			}
			if ( kSourceView[kStepView[half]] > kWanted ) end = half;
			else start = half + 1;
		}
		if ( start < kStepCount && kSourceView[kStepView[start]] == kWanted )
		{
			result = start;
		}
		if ( result == -2 )
		{
			resultView[ cellWanted ] = -2;
			return;
		}
		
		start = kStepView[ result ];
		end = sourceCellCount;
		if ( result < kStepCount - 1 ) end = kStepView[ result + 1 ];
		result = -2;
		
		// Now find the exact match for j and i in the already reduced interval
		// We search for both simultaneously to avoid writing multiple binary searches
		while ( end > start )
		{
			int half = start + ( end - start ) / 2;
			int jHalf = jSourceView[ half ];
			int iHalf = iSourceView[ half ];
			
			if ( jHalf == jWanted && iHalf == iWanted )
			{
				result = half;
				break;
			}
			
			// Move the end bound if j is too big, OR if j matches but i is too big
			if ( jHalf > jWanted || ( jHalf == jWanted && iHalf > iWanted ) ) end = half;
			else start = half + 1;
		}
		
		if ( result == -2 && start < sourceCellCount )
		{
			if ( jSourceView[ start ] == jWanted && iSourceView[ start ] == iWanted )
			{
				result = start;
			}
		}
		
		resultView[ cellWanted ] = result;
	};
	
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, wantedCellCount, cellLambda );
}

// Do a repetetive fluid marking in order to find which of our coarse cells could possibly contain at least one finest (maximum refinement level) fluid cell. 
// To do this we will be repeatedly temporarily shifting the origin of our grid, to simulate being on a finer grid.
void markWhereFinestFluidIs( BoolArrayType &fluidMarkerArray, IJKArrayStruct &IJK, std::vector<STLStruct> STLs, InfoStruct &Info )
{
	const int shiftCount = std::pow(2, (gridLevelCount - Info.gridID - 1));
	const float finestRes = Info.res / ( std::pow(2, (gridLevelCount - Info.gridID - 1)) );
	const float shiftStart = - 0.5f * Info.res + 0.5f * finestRes;
	const float oxOriginal = Info.ox;
	const float oyOriginal = Info.oy;
	const float ozOriginal = Info.oz;
	fluidMarkerArray.setValue( 0 );
	BoolArrayType markerArray = BoolArrayType( Info.cellCount );
	markerArray.setValue( 0 );
	BoolArrayType markerArraySTL = BoolArrayType( Info.cellCount );
	markerArraySTL.setValue( 0 );
	for ( int shiftCountX = 0; shiftCountX < shiftCount; shiftCountX++ )
	{
		for ( int shiftCountY = 0; shiftCountY < shiftCount; shiftCountY++ )
		{
			for ( int shiftCountZ = 0; shiftCountZ < shiftCount; shiftCountZ++ )
			{
				Info.ox = oxOriginal + shiftStart + shiftCountX * finestRes;
				Info.oy = oyOriginal + shiftStart + shiftCountY * finestRes;
				Info.oz = ozOriginal + shiftStart + shiftCountZ * finestRes;				
				ApplyMarkersFromFunction( markerArray, IJK, Info );
				for ( int STLIndex = 0; STLIndex < STLs.size(); STLIndex++ )
				{
					const bool insideMarkerValue = 1;
					ApplyMarkersInsideSTL( markerArraySTL, IJK, STLs[STLIndex], insideMarkerValue, Info );
					sumBoolArrays( markerArray, markerArraySTL, markerArray );
				}
				invertBoolArray( markerArray );
				sumBoolArrays( fluidMarkerArray, markerArray, fluidMarkerArray );
			}
		}
	}
	Info.ox = oxOriginal;
	Info.oy = oyOriginal;
	Info.oz = ozOriginal;
}

void getDIADNeighbours( IJKArrayStruct &IJK, DIADNeighboursStruct &Neighbours )
{
	const int cellCount = IJK.iArray.getSize();
	Neighbours.iPlusArray.setSize(cellCount);
	Neighbours.jPlusArray.setSize(cellCount);
	Neighbours.kPlusArray.setSize(cellCount);
	Neighbours.iMinusArray.setSize(cellCount);
	Neighbours.jMinusArray.setSize(cellCount);
	Neighbours.kMinusArray.setSize(cellCount);
	Neighbours.iPlusArray.setValue( -2 );
	Neighbours.jPlusArray.setValue( -2 );
	Neighbours.kPlusArray.setValue( -2 );
	Neighbours.iMinusArray.setValue( -2 );
	Neighbours.jMinusArray.setValue( -2 );
	Neighbours.kMinusArray.setValue( -2 );
	
	IntArrayType kStepArray;
	getIJKStepArray( IJK, kStepArray );

	IJKArrayStruct IJKShifted;
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 1, 0, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.iPlusArray, kStepArray );

	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 1, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.jPlusArray, kStepArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 0, 1 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.kPlusArray, kStepArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, -1, 0, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.iMinusArray, kStepArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, -1, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.jMinusArray, kStepArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 0, -1 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.kMinusArray, kStepArray );
}

void getDIADEsotwistConnections( DIADGridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	IJKArrayStruct IJK = Grid.IJK;
	const int cellCount = Info.cellCount;
	Grid.EsotwistConnections.iNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.jNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.kNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.ijNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.ikNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.jkNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.ijkNbrArray.setSize( cellCount );
	Grid.EsotwistConnections.iNbrArray.setValue( -2 );
	Grid.EsotwistConnections.jNbrArray.setValue( -2 );
	Grid.EsotwistConnections.kNbrArray.setValue( -2 );
	Grid.EsotwistConnections.ijNbrArray.setValue( -2 );
	Grid.EsotwistConnections.ikNbrArray.setValue( -2 );
	Grid.EsotwistConnections.jkNbrArray.setValue( -2 );
	Grid.EsotwistConnections.ijkNbrArray.setValue( -2 );
	
	IntArrayType kStepArray;
	getIJKStepArray( IJK, kStepArray );
	
	IJKArrayStruct IJKShifted;
	int invalidNeighbourCount;
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 0, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.iNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.iNbrArray );
	}

	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 1, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.jNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.jNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 0, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.kNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.kNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 1, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.ijNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.ijNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 0, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.ikNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.ikNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 1, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.jkNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.jkNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 1, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistConnections.ijkNbrArray, kStepArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistConnections.ijkNbrArray );
	}
}

// Mark target cell as 1 if at least one of its neighbours in source is 1
void markDIADNeighbours( DIADNeighboursStruct &Neighbours, BoolArrayType &sourceArray, BoolArrayType &targetArray )
{
	const int cellCount = Neighbours.iPlusArray.getSize();
	
	auto iPlusView = Neighbours.iPlusArray.getConstView();
	auto jPlusView = Neighbours.jPlusArray.getConstView();
	auto kPlusView = Neighbours.kPlusArray.getConstView();
	auto iMinusView = Neighbours.iMinusArray.getConstView();
	auto jMinusView = Neighbours.jMinusArray.getConstView();
	auto kMinusView = Neighbours.kMinusArray.getConstView();
	
	auto sourceView = sourceArray.getConstView();
	auto targetView = targetArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		const int ip = iPlusView[ cell ];
		const int jp = jPlusView[ cell ];
		const int kp = kPlusView[ cell ];
		const int im = iMinusView[ cell ];
		const int jm = jMinusView[ cell ];
		const int km = kMinusView[ cell ];
		// Now I need to browse all directions including all diagonals (the iPlus, jPlus etc are indexes of neighbours only in the main directions)
		// The neighbour arrays can also contain indexes lower than zero, that is when there no neighbour exists
		
		// ---------------------------------------------------------
        // Direct neighbors (6 Faces)
        // ---------------------------------------------------------
        if (ip >= 0 && sourceView[ip]) { targetView[cell] = true; return; }
        if (im >= 0 && sourceView[im]) { targetView[cell] = true; return; }
        if (jp >= 0 && sourceView[jp]) { targetView[cell] = true; return; }
        if (jm >= 0 && sourceView[jm]) { targetView[cell] = true; return; }
        if (kp >= 0 && sourceView[kp]) { targetView[cell] = true; return; }
        if (km >= 0 && sourceView[km]) { targetView[cell] = true; return; }

        // ---------------------------------------------------------
        // Direct diagonals (12 edges) and indirect diagonals (8 corners)
        // ---------------------------------------------------------
		// 1. Traverse starting along X (+i, -i)
		// Covers paths: X->Y->Z, X->Z->Y and X->Y, X->Z
		// ---------------------------------------------------------
		if (ip >= 0) 
		{
			// Path: X -> Z -> Y
			const int ip_kp = kPlusView[ip];
			if (ip_kp >= 0) {
				if (sourceView[ip_kp]) { targetView[cell] = true; return; }
				const int ip_kp_jp = jPlusView[ip_kp];
				if (ip_kp_jp >= 0 && sourceView[ip_kp_jp]) { targetView[cell] = true; return; }
				const int ip_kp_jm = jMinusView[ip_kp];
				if (ip_kp_jm >= 0 && sourceView[ip_kp_jm]) { targetView[cell] = true; return; }
			}
			const int ip_km = kMinusView[ip];
			if (ip_km >= 0) {
				if (sourceView[ip_km]) { targetView[cell] = true; return; }
				const int ip_km_jp = jPlusView[ip_km];
				if (ip_km_jp >= 0 && sourceView[ip_km_jp]) { targetView[cell] = true; return; }
				const int ip_km_jm = jMinusView[ip_km];
				if (ip_km_jm >= 0 && sourceView[ip_km_jm]) { targetView[cell] = true; return; }
			}
			
			// Path: X -> Y -> Z
			const int ip_jp = jPlusView[ip];
			if (ip_jp >= 0) {
				if (sourceView[ip_jp]) { targetView[cell] = true; return; }
				const int ip_jp_kp = kPlusView[ip_jp];
				if (ip_jp_kp >= 0 && sourceView[ip_jp_kp]) { targetView[cell] = true; return; }
				const int ip_jp_km = kMinusView[ip_jp];
				if (ip_jp_km >= 0 && sourceView[ip_jp_km]) { targetView[cell] = true; return; }
			}
			const int ip_jm = jMinusView[ip];
			if (ip_jm >= 0) {
				if (sourceView[ip_jm]) { targetView[cell] = true; return; }
				const int ip_jm_kp = kPlusView[ip_jm];
				if (ip_jm_kp >= 0 && sourceView[ip_jm_kp]) { targetView[cell] = true; return; }
				const int ip_jm_km = kMinusView[ip_jm];
				if (ip_jm_km >= 0 && sourceView[ip_jm_km]) { targetView[cell] = true; return; }
			}
		}

		if (im >= 0) 
		{
			// Path: -X -> Z -> Y
			const int im_kp = kPlusView[im];
			if (im_kp >= 0) {
				if (sourceView[im_kp]) { targetView[cell] = true; return; }
				const int im_kp_jp = jPlusView[im_kp];
				if (im_kp_jp >= 0 && sourceView[im_kp_jp]) { targetView[cell] = true; return; }
				const int im_kp_jm = jMinusView[im_kp];
				if (im_kp_jm >= 0 && sourceView[im_kp_jm]) { targetView[cell] = true; return; }
			}
			const int im_km = kMinusView[im];
			if (im_km >= 0) {
				if (sourceView[im_km]) { targetView[cell] = true; return; }
				const int im_km_jp = jPlusView[im_km];
				if (im_km_jp >= 0 && sourceView[im_km_jp]) { targetView[cell] = true; return; }
				const int im_km_jm = jMinusView[im_km];
				if (im_km_jm >= 0 && sourceView[im_km_jm]) { targetView[cell] = true; return; }
			}
			
			// Path: -X -> Y -> Z
			const int im_jp = jPlusView[im];
			if (im_jp >= 0) {
				if (sourceView[im_jp]) { targetView[cell] = true; return; }
				const int im_jp_kp = kPlusView[im_jp];
				if (im_jp_kp >= 0 && sourceView[im_jp_kp]) { targetView[cell] = true; return; }
				const int im_jp_km = kMinusView[im_jp];
				if (im_jp_km >= 0 && sourceView[im_jp_km]) { targetView[cell] = true; return; }
			}
			const int im_jm = jMinusView[im];
			if (im_jm >= 0) {
				if (sourceView[im_jm]) { targetView[cell] = true; return; }
				const int im_jm_kp = kPlusView[im_jm];
				if (im_jm_kp >= 0 && sourceView[im_jm_kp]) { targetView[cell] = true; return; }
				const int im_jm_km = kMinusView[im_jm];
				if (im_jm_km >= 0 && sourceView[im_jm_km]) { targetView[cell] = true; return; }
			}
		}

		// ---------------------------------------------------------
		// 2. Traverse starting along Y (+j, -j)
		// Covers paths: Y->X->Z, Y->Z->X and Y->X, Y->Z
		// ---------------------------------------------------------
		if (jp >= 0) 
		{
			// Path: Y -> X -> Z
			const int jp_ip = iPlusView[jp];
			if (jp_ip >= 0) {
				if (sourceView[jp_ip]) { targetView[cell] = true; return; }
				const int jp_ip_kp = kPlusView[jp_ip];
				if (jp_ip_kp >= 0 && sourceView[jp_ip_kp]) { targetView[cell] = true; return; }
				const int jp_ip_km = kMinusView[jp_ip];
				if (jp_ip_km >= 0 && sourceView[jp_ip_km]) { targetView[cell] = true; return; }
			}
			const int jp_im = iMinusView[jp];
			if (jp_im >= 0) {
				if (sourceView[jp_im]) { targetView[cell] = true; return; }
				const int jp_im_kp = kPlusView[jp_im];
				if (jp_im_kp >= 0 && sourceView[jp_im_kp]) { targetView[cell] = true; return; }
				const int jp_im_km = kMinusView[jp_im];
				if (jp_im_km >= 0 && sourceView[jp_im_km]) { targetView[cell] = true; return; }
			}

			// Path: Y -> Z -> X
			const int jp_kp = kPlusView[jp];
			if (jp_kp >= 0) {
				if (sourceView[jp_kp]) { targetView[cell] = true; return; }
				const int jp_kp_ip = iPlusView[jp_kp];
				if (jp_kp_ip >= 0 && sourceView[jp_kp_ip]) { targetView[cell] = true; return; }
				const int jp_kp_im = iMinusView[jp_kp];
				if (jp_kp_im >= 0 && sourceView[jp_kp_im]) { targetView[cell] = true; return; }
			}
			const int jp_km = kMinusView[jp];
			if (jp_km >= 0) {
				if (sourceView[jp_km]) { targetView[cell] = true; return; }
				const int jp_km_ip = iPlusView[jp_km];
				if (jp_km_ip >= 0 && sourceView[jp_km_ip]) { targetView[cell] = true; return; }
				const int jp_km_im = iMinusView[jp_km];
				if (jp_km_im >= 0 && sourceView[jp_km_im]) { targetView[cell] = true; return; }
			}
		}

		if (jm >= 0) 
		{
			// Path: -Y -> X -> Z
			const int jm_ip = iPlusView[jm];
			if (jm_ip >= 0) {
				if (sourceView[jm_ip]) { targetView[cell] = true; return; }
				const int jm_ip_kp = kPlusView[jm_ip];
				if (jm_ip_kp >= 0 && sourceView[jm_ip_kp]) { targetView[cell] = true; return; }
				const int jm_ip_km = kMinusView[jm_ip];
				if (jm_ip_km >= 0 && sourceView[jm_ip_km]) { targetView[cell] = true; return; }
			}
			const int jm_im = iMinusView[jm];
			if (jm_im >= 0) {
				if (sourceView[jm_im]) { targetView[cell] = true; return; }
				const int jm_im_kp = kPlusView[jm_im];
				if (jm_im_kp >= 0 && sourceView[jm_im_kp]) { targetView[cell] = true; return; }
				const int jm_im_km = kMinusView[jm_im];
				if (jm_im_km >= 0 && sourceView[jm_im_km]) { targetView[cell] = true; return; }
			}

			// Path: -Y -> Z -> X
			const int jm_kp = kPlusView[jm];
			if (jm_kp >= 0) {
				if (sourceView[jm_kp]) { targetView[cell] = true; return; }
				const int jm_kp_ip = iPlusView[jm_kp];
				if (jm_kp_ip >= 0 && sourceView[jm_kp_ip]) { targetView[cell] = true; return; }
				const int jm_kp_im = iMinusView[jm_kp];
				if (jm_kp_im >= 0 && sourceView[jm_kp_im]) { targetView[cell] = true; return; }
			}
			const int jm_km = kMinusView[jm];
			if (jm_km >= 0) {
				if (sourceView[jm_km]) { targetView[cell] = true; return; }
				const int jm_km_ip = iPlusView[jm_km];
				if (jm_km_ip >= 0 && sourceView[jm_km_ip]) { targetView[cell] = true; return; }
				const int jm_km_im = iMinusView[jm_km];
				if (jm_km_im >= 0 && sourceView[jm_km_im]) { targetView[cell] = true; return; }
			}
		}

		// ---------------------------------------------------------
		// 3. Traverse starting along Z (+k, -k)
		// Covers paths: Z->X->Y, Z->Y->X and Z->X, Z->Y
		// ---------------------------------------------------------
		if (kp >= 0) 
		{
			// Path: Z -> X -> Y
			const int kp_ip = iPlusView[kp];
			if (kp_ip >= 0) {
				if (sourceView[kp_ip]) { targetView[cell] = true; return; }
				const int kp_ip_jp = jPlusView[kp_ip];
				if (kp_ip_jp >= 0 && sourceView[kp_ip_jp]) { targetView[cell] = true; return; }
				const int kp_ip_jm = jMinusView[kp_ip];
				if (kp_ip_jm >= 0 && sourceView[kp_ip_jm]) { targetView[cell] = true; return; }
			}
			const int kp_im = iMinusView[kp];
			if (kp_im >= 0) {
				if (sourceView[kp_im]) { targetView[cell] = true; return; }
				const int kp_im_jp = jPlusView[kp_im];
				if (kp_im_jp >= 0 && sourceView[kp_im_jp]) { targetView[cell] = true; return; }
				const int kp_im_jm = jMinusView[kp_im];
				if (kp_im_jm >= 0 && sourceView[kp_im_jm]) { targetView[cell] = true; return; }
			}

			// Path: Z -> Y -> X
			const int kp_jp = jPlusView[kp];
			if (kp_jp >= 0) {
				if (sourceView[kp_jp]) { targetView[cell] = true; return; }
				const int kp_jp_ip = iPlusView[kp_jp];
				if (kp_jp_ip >= 0 && sourceView[kp_jp_ip]) { targetView[cell] = true; return; }
				const int kp_jp_im = iMinusView[kp_jp];
				if (kp_jp_im >= 0 && sourceView[kp_jp_im]) { targetView[cell] = true; return; }
			}
			const int kp_jm = jMinusView[kp];
			if (kp_jm >= 0) {
				if (sourceView[kp_jm]) { targetView[cell] = true; return; }
				const int kp_jm_ip = iPlusView[kp_jm];
				if (kp_jm_ip >= 0 && sourceView[kp_jm_ip]) { targetView[cell] = true; return; }
				const int kp_jm_im = iMinusView[kp_jm];
				if (kp_jm_im >= 0 && sourceView[kp_jm_im]) { targetView[cell] = true; return; }
			}
		}

		if (km >= 0) 
		{
			// Path: -Z -> X -> Y
			const int km_ip = iPlusView[km];
			if (km_ip >= 0) {
				if (sourceView[km_ip]) { targetView[cell] = true; return; }
				const int km_ip_jp = jPlusView[km_ip];
				if (km_ip_jp >= 0 && sourceView[km_ip_jp]) { targetView[cell] = true; return; }
				const int km_ip_jm = jMinusView[km_ip];
				if (km_ip_jm >= 0 && sourceView[km_ip_jm]) { targetView[cell] = true; return; }
			}
			const int km_im = iMinusView[km];
			if (km_im >= 0) {
				if (sourceView[km_im]) { targetView[cell] = true; return; }
				const int km_im_jp = jPlusView[km_im];
				if (km_im_jp >= 0 && sourceView[km_im_jp]) { targetView[cell] = true; return; }
				const int km_im_jm = jMinusView[km_im];
				if (km_im_jm >= 0 && sourceView[km_im_jm]) { targetView[cell] = true; return; }
			}

			// Path: -Z -> Y -> X
			const int km_jp = jPlusView[km];
			if (km_jp >= 0) {
				if (sourceView[km_jp]) { targetView[cell] = true; return; }
				const int km_jp_ip = iPlusView[km_jp];
				if (km_jp_ip >= 0 && sourceView[km_jp_ip]) { targetView[cell] = true; return; }
				const int km_jp_im = iMinusView[km_jp];
				if (km_jp_im >= 0 && sourceView[km_jp_im]) { targetView[cell] = true; return; }
			}
			const int km_jm = jMinusView[km];
			if (km_jm >= 0) {
				if (sourceView[km_jm]) { targetView[cell] = true; return; }
				const int km_jm_ip = iPlusView[km_jm];
				if (km_jm_ip >= 0 && sourceView[km_jm_ip]) { targetView[cell] = true; return; }
				const int km_jm_im = iMinusView[km_jm];
				if (km_jm_im >= 0 && sourceView[km_jm_im]) { targetView[cell] = true; return; }
			}
		}
    };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, cellCount, cellLambda );
}

void buildDIADGrid( std::vector<DIADGridStruct> &grids, std::vector<STLStruct> STLs, const int level )
{
	DIADGridStruct &Grid = grids[level];
	std::cout << "Initial cell count on level " << level <<" : " << Grid.Info.cellCount << std::endl;
	
	const int kCut = Grid.Info.cellCountZ / 2;

	sortIJK( Grid.IJK );
	DIADNeighboursStruct Neighbours;
	getDIADNeighbours( Grid.IJK, Neighbours );	
		
	BoolArrayType fluidMarkerArray = BoolArrayType( Grid.Info.cellCount );
	markWhereFinestFluidIs( fluidMarkerArray, Grid.IJK, STLs, Grid.Info );
	
	// BORDER
	BoolArrayType borderMarkerArray = BoolArrayType( Grid.Info.cellCount );
	borderMarkerArray.setValue( 0 );
	markDIADNeighbours( Neighbours, fluidMarkerArray, borderMarkerArray );
	BoolArrayType fluidInverseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	fluidInverseMarkerArray = fluidMarkerArray;
	invertBoolArray( fluidInverseMarkerArray );
	multiplyBoolArrays( fluidInverseMarkerArray, borderMarkerArray, borderMarkerArray );
	
	// KEEP CELL ARRAY -> SO FAR KEEP FLUID AND BORDER
	BoolArrayType keepCellMarkerArray = BoolArrayType( Grid.Info.cellCount );
	keepCellMarkerArray = fluidMarkerArray;
	sumBoolArrays( keepCellMarkerArray, borderMarkerArray, keepCellMarkerArray );
	
	// FINEST GRID BRANCH
	if ( level == gridLevelCount - 1 )
	{
		Grid.bouncebackMarkerArray = borderMarkerArray;
		sortFinestGrid( Grid, keepCellMarkerArray );
		Grid.Info.cellCount = countMarkerCells( keepCellMarkerArray );
		Grid.IJK.iArray.resize(Grid.Info.cellCount);
		Grid.IJK.jArray.resize(Grid.Info.cellCount);
		Grid.IJK.kArray.resize(Grid.Info.cellCount);
		Grid.bouncebackMarkerArray.resize(Grid.Info.cellCount);
		getDIADEsotwistConnections( Grid );
		//PLOT
		exportSectionCutPlotXY( Grid, kCut, level );
		system("python3 ../include/plotter/plotter.py");
		//PLOT END
		std::cout << "Final cell count on level " << level <<" : " << Grid.Info.cellCount << std::endl;
		return;
	}
	
	// COARSE GRID BRANCH
	// THICK REFINEMENT REGION
	BoolArrayType refinementMarkerArray = BoolArrayType( Grid.Info.cellCount );
	refinementMarkerArray = borderMarkerArray;
	BoolArrayType newRefinementMarkerArray = BoolArrayType( Grid.Info.cellCount );
	newRefinementMarkerArray.setValue( 0 );
	const int thickness = std::max({wallRefinementSpan, 2}) + (gridLevelCount - level - 1);
	for ( int layer = 0; layer < thickness; layer++ )
	{
		markDIADNeighbours( Neighbours, refinementMarkerArray, newRefinementMarkerArray );
		multiplyBoolArrays( newRefinementMarkerArray, keepCellMarkerArray, refinementMarkerArray );
	}
	
	// FINE TO COARSE COMMUNICATION INTERFACE
	BoolArrayType fineToCoarseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	fineToCoarseMarkerArray.setValue( 0 );
	markDIADNeighbours( Neighbours, refinementMarkerArray, fineToCoarseMarkerArray );
	multiplyBoolArrays( fluidMarkerArray, fineToCoarseMarkerArray, fineToCoarseMarkerArray );
	BoolArrayType refinementInverseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	refinementInverseMarkerArray = refinementMarkerArray;
	invertBoolArray( refinementInverseMarkerArray );
	multiplyBoolArrays( fineToCoarseMarkerArray, refinementInverseMarkerArray, fineToCoarseMarkerArray );
	
	// COARSE TO FINE COMMUNICATION INTERFACE
	BoolArrayType coarseToFineMarkerArray = BoolArrayType( Grid.Info.cellCount );
	coarseToFineMarkerArray.setValue( 0 );
	markDIADNeighbours( Neighbours, fineToCoarseMarkerArray, coarseToFineMarkerArray );
	multiplyBoolArrays( fluidMarkerArray, coarseToFineMarkerArray, coarseToFineMarkerArray );
	BoolArrayType fineToCoarseInverseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	fineToCoarseInverseMarkerArray = fineToCoarseMarkerArray;
	invertBoolArray( fineToCoarseInverseMarkerArray );
	multiplyBoolArrays( coarseToFineMarkerArray, fineToCoarseInverseMarkerArray, coarseToFineMarkerArray );
	multiplyBoolArrays( coarseToFineMarkerArray, refinementInverseMarkerArray, coarseToFineMarkerArray );
	
	// NOW ALSO ADD THE INTERFACE INTO THE REFINEMENT AREA
	sumBoolArrays( refinementMarkerArray, fineToCoarseMarkerArray, refinementMarkerArray );
	sumBoolArrays( refinementMarkerArray, coarseToFineMarkerArray, refinementMarkerArray );
	refinementInverseMarkerArray = refinementMarkerArray;
	invertBoolArray( refinementInverseMarkerArray );
	
	// ASSEMBLE IJK FOR THE NEXT LEVEL FINER GRID BASED ON THE REFINEMENT AREA
	const int cellCountFine = 8 * countMarkerCells( refinementMarkerArray );
	grids[level+1].Info.cellCount = cellCountFine;
	BoolArrayTypeCPU refinementMarkerArrayCPU = BoolArrayTypeCPU( Grid.Info.cellCount );
	refinementMarkerArrayCPU = refinementMarkerArray;
	IJKArrayStructCPU IJKCPU = IJKArrayStructCPU( Grid.IJK );
	IJKArrayStructCPU IJKFineCPU;
	IJKFineCPU.iArray.setSize( cellCountFine );
	IJKFineCPU.jArray.setSize( cellCountFine );
	IJKFineCPU.kArray.setSize( cellCountFine );
	int cellFine = 0;
	for ( int cellCoarse = 0; cellCoarse < Grid.Info.cellCount; cellCoarse++ )
	{
		if ( refinementMarkerArrayCPU[ cellCoarse ] )
		{
			const int iCoarse = IJKCPU.iArray[ cellCoarse ];
			const int jCoarse = IJKCPU.jArray[ cellCoarse ];
			const int kCoarse = IJKCPU.kArray[ cellCoarse ];
			const int iFine = iCoarse * 2;
			const int jFine = jCoarse * 2;
			const int kFine = kCoarse * 2;
			IJKFineCPU.iArray[ cellFine ] = iFine;
			IJKFineCPU.jArray[ cellFine ] = jFine;
			IJKFineCPU.kArray[ cellFine ] = kFine;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine + 1;
			IJKFineCPU.jArray[ cellFine ] = jFine;
			IJKFineCPU.kArray[ cellFine ] = kFine;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine;
			IJKFineCPU.jArray[ cellFine ] = jFine + 1;
			IJKFineCPU.kArray[ cellFine ] = kFine;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine + 1;
			IJKFineCPU.jArray[ cellFine ] = jFine + 1;
			IJKFineCPU.kArray[ cellFine ] = kFine;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine;
			IJKFineCPU.jArray[ cellFine ] = jFine;
			IJKFineCPU.kArray[ cellFine ] = kFine + 1;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine + 1;
			IJKFineCPU.jArray[ cellFine ] = jFine;
			IJKFineCPU.kArray[ cellFine ] = kFine + 1;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine;
			IJKFineCPU.jArray[ cellFine ] = jFine + 1;
			IJKFineCPU.kArray[ cellFine ] = kFine + 1;
			cellFine++;
			IJKFineCPU.iArray[ cellFine ] = iFine + 1;
			IJKFineCPU.jArray[ cellFine ] = jFine + 1;
			IJKFineCPU.kArray[ cellFine ] = kFine + 1;
			cellFine++;
		}
	}
	
	grids[level + 1].IJK = IJKArrayStruct( IJKFineCPU );
	
	// MARK KEEP CELLS TO BE ABLE TO DELETE ALL CELLS THAT ARE NOT NEEDED ON OUR COARSE GRID
	keepCellMarkerArray = fluidMarkerArray;
	multiplyBoolArrays( keepCellMarkerArray, refinementInverseMarkerArray, keepCellMarkerArray );
	sumBoolArrays( keepCellMarkerArray, coarseToFineMarkerArray, keepCellMarkerArray );
	sumBoolArrays( keepCellMarkerArray, fineToCoarseMarkerArray, keepCellMarkerArray );
	
	// SORT ALL CELLS BY KEY: keepCell, k, j, i
	sortCoarseGrid( Grid, keepCellMarkerArray, fineToCoarseMarkerArray, coarseToFineMarkerArray );
	Grid.Info.cellCount = countMarkerCells( keepCellMarkerArray );
	Grid.IJK.iArray.resize(Grid.Info.cellCount);
	Grid.IJK.jArray.resize(Grid.Info.cellCount);
	Grid.IJK.kArray.resize(Grid.Info.cellCount);
	fineToCoarseMarkerArray.resize(Grid.Info.cellCount);
	coarseToFineMarkerArray.resize(Grid.Info.cellCount);
	
	// FINAL BOUNCEBACK PASS IN CASE WE ARE NOT REFINING IN SOME WALL AREA
	Grid.bouncebackMarkerArray.setSize(Grid.Info.cellCount);
	Grid.bouncebackMarkerArray.setValue( 0 );
	BoolArrayType markerArraySTL = BoolArrayType( Grid.Info.cellCount );	
	ApplyMarkersFromFunction( Grid.bouncebackMarkerArray, Grid.IJK, Grid.Info );
	for ( int STLIndex = 0; STLIndex < STLs.size(); STLIndex++ )
	{
		const bool insideMarkerValue = 1;
		ApplyMarkersInsideSTL( markerArraySTL, Grid.IJK, STLs[STLIndex], insideMarkerValue, Grid.Info );
		sumBoolArrays( Grid.bouncebackMarkerArray, markerArraySTL, Grid.bouncebackMarkerArray );
	}
	
	// BUILD ESOTWIST CONNECTIONS
	getDIADEsotwistConnections( Grid ); // <- NEED TO FINISH THIS!
	
	// CALL THE FUNCTION RECURSIVELY TO BUILD THE NEXT FINER GRID LEVEL
	grids[level+1].Info.gridID = Grid.Info.gridID + 1;
	grids[level+1].Info.res = Grid.Info.res * 0.5f;
	grids[level+1].Info.dtPhys = Grid.Info.dtPhys * 0.5f;
	grids[level+1].Info.nu = (grids[level+1].Info.dtPhys * nuPhys) / ((grids[level+1].Info.res/1000.f) * (grids[level+1].Info.res/1000.f));
	grids[level+1].Info.ox = Grid.Info.ox - grids[level+1].Info.res * 0.5f;
	grids[level+1].Info.oy = Grid.Info.oy - grids[level+1].Info.res * 0.5f;
	grids[level+1].Info.oz = Grid.Info.oz - grids[level+1].Info.res * 0.5f;
	grids[level+1].Info.cellCountX = Grid.Info.cellCountX * 2;
	grids[level+1].Info.cellCountY = Grid.Info.cellCountY * 2;
	grids[level+1].Info.cellCountZ = Grid.Info.cellCountZ * 2;
	buildDIADGrid( grids, STLs, level + 1 );
	
	// MAP INTERFACE COMMUNICATION PREPARATION
	auto iView = Grid.IJK.iArray.getView();
	auto jView = Grid.IJK.jArray.getView();
	auto kView = Grid.IJK.kArray.getView();
	IJKCPU = IJKArrayStructCPU( Grid.IJK );
	IJKArrayStruct fineIJK;
	const int coarseToFineCount = countMarkerCells( coarseToFineMarkerArray );
	Grid.coarseToFineWriteArray.setSize( coarseToFineCount );
	Grid.coarseToFineReadArray.setSize( coarseToFineCount );
	Grid.coarseToFineWriteArray.setValue( -2 );
	Grid.coarseToFineReadArray.setValue( -2 );
	const int fineToCoarseCount = countMarkerCells( fineToCoarseMarkerArray );
	Grid.fineToCoarseWriteArray.setSize( fineToCoarseCount );
	Grid.fineToCoarseReadArray.setSize( fineToCoarseCount );
	Grid.fineToCoarseWriteArray.setValue( -2 );
	Grid.fineToCoarseReadArray.setValue( -2 );
	int counter;
	
	// MAP INTERFACE COMMUNICATION COARSE TO FINE
	BoolArrayTypeCPU coarseToFineMarkerArrayCPU;
	coarseToFineMarkerArrayCPU = coarseToFineMarkerArray;
	IntArrayTypeCPU coarseToFineReadArrayCPU;
	coarseToFineReadArrayCPU = Grid.coarseToFineReadArray;
	counter = 0;
	for ( int cell = 0; cell < Grid.Info.cellCount; cell++ )
	{
		if ( coarseToFineMarkerArrayCPU[ cell ] )
		{
			coarseToFineReadArrayCPU[ counter ] = cell;
			counter++;
		}
	}
	Grid.coarseToFineReadArray = coarseToFineReadArrayCPU;
	fineIJK.iArray.setSize( coarseToFineCount );
	fineIJK.jArray.setSize( coarseToFineCount );
	fineIJK.kArray.setSize( coarseToFineCount );
	auto coarseToFineReadArrayView = Grid.coarseToFineReadArray.getConstView();
	auto iFineView1 = fineIJK.iArray.getView();
	auto jFineView1 = fineIJK.jArray.getView();
	auto kFineView1 = fineIJK.kArray.getView();
	auto cellLambda1 = [=] __cuda_callable__ ( const int counter ) mutable
	{
		const int cell = coarseToFineReadArrayView[ counter ];
		iFineView1[ counter ] = iView[ cell ] * 2;
		jFineView1[ counter ] = jView[ cell ] * 2;
		kFineView1[ counter ] = kView[ cell ] * 2;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, coarseToFineCount, cellLambda1 );
	IntArrayType kStepArray;
	getIJKStepArray( grids[level+1].IJK, kStepArray );
	findMatchingIJKIndex( fineIJK, grids[level+1].IJK, Grid.coarseToFineWriteArray, kStepArray );
	const int invalidConnections1 = countInvalidIndexes( Grid.coarseToFineWriteArray );
	if ( invalidConnections1 > 0 )
	{
		throw std::runtime_error( 
			"DIAD Grid Generation Error: " + std::to_string(invalidConnections1) + 
			" invalid coarse-to-fine connections found on level " + std::to_string(level) 
		);
	}
	
	// MAP INTERFACE COMMUNICATION FINE TO COARSE
	BoolArrayTypeCPU fineToCoarseMarkerArrayCPU;
	fineToCoarseMarkerArrayCPU = fineToCoarseMarkerArray;
	IntArrayTypeCPU fineToCoarseWriteArrayCPU;
	fineToCoarseWriteArrayCPU = Grid.fineToCoarseWriteArray;
	counter = 0;
	for ( int cell = 0; cell < Grid.Info.cellCount; cell++ )
	{
		if ( fineToCoarseMarkerArrayCPU[ cell ] )
		{
			fineToCoarseWriteArrayCPU[ counter ] = cell;
			counter++;
		}
	}
	Grid.fineToCoarseWriteArray = fineToCoarseWriteArrayCPU;
	fineIJK.iArray.setSize( fineToCoarseCount );
	fineIJK.jArray.setSize( fineToCoarseCount );
	fineIJK.kArray.setSize( fineToCoarseCount );
	auto fineToCoarseWriteArrayView = Grid.fineToCoarseWriteArray.getConstView();
	auto iFineView2 = fineIJK.iArray.getView();
	auto jFineView2 = fineIJK.jArray.getView();
	auto kFineView2 = fineIJK.kArray.getView();
	auto cellLambda2 = [=] __cuda_callable__ ( const int counter ) mutable
	{
		const int cell = fineToCoarseWriteArrayView[ counter ];
		iFineView2[ counter ] = iView[ cell ] * 2;
		jFineView2[ counter ] = jView[ cell ] * 2;
		kFineView2[ counter ] = kView[ cell ] * 2;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, fineToCoarseCount, cellLambda2 );
	getIJKStepArray( grids[level+1].IJK, kStepArray );
	findMatchingIJKIndex( fineIJK, grids[level+1].IJK, Grid.fineToCoarseReadArray, kStepArray );
	const int invalidConnections2 = countInvalidIndexes( Grid.fineToCoarseReadArray );
	if ( invalidConnections2 > 0 )
	{
		throw std::runtime_error( 
			"DIAD Grid Generation Error: " + std::to_string(invalidConnections2) + 
			" invalid coarse-to-fine connections found on level " + std::to_string(level) 
		);
	}
	
	//PLOT
	exportSectionCutPlotXY( Grid, kCut, level );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	std::cout << "Final cell count on level " << level <<" : " << Grid.Info.cellCount << std::endl;
}
