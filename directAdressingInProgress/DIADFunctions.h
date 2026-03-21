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

// Receives sorted IJKSource that !must! be already sorted with key k, j, i so that k index changes the slowest
// For each Wanted, find a matching cell in Source. If its found, write its index to resultArray. If there is no such cell, write -2.		
void findMatchingIJKIndex( 	IJKArrayStruct &Wanted, IJKArrayStruct &Source, IntArrayType &resultArray )
{
	const int wantedCellCount = Wanted.iArray.getSize();
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
	IntArrayType kStepArray;
	kStepArray = kStepArrayCPU;
	
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

	IJKArrayStruct IJKShifted;
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 1, 0, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.iPlusArray );

	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 1, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.jPlusArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 0, 1 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.kPlusArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, -1, 0, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.iMinusArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, -1, 0 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.jMinusArray );
	
	IJKShifted = IJK;
	shiftIJK( IJKShifted, 0, 0, -1 );
	findMatchingIJKIndex( IJKShifted, IJK, Neighbours.kMinusArray );
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
	DIADGridStruct Grid = grids[level];
	InfoStruct Info = Grid.Info;
	
	const int kCut = Grid.Info.cellCountZ / 2;
	
	//PLOT
	Grid.bouncebackMarkerArray = BoolArrayType( Grid.Info.cellCount );
	Grid.bouncebackMarkerArray.setValue( 0 );
	exportSectionCutPlotXY( Grid, kCut, 10*level );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END

	sortIJK( Grid.IJK );
	DIADNeighboursStruct Neighbours;
	getDIADNeighbours( Grid.IJK, Neighbours );	
		
	BoolArrayType fluidMarkerArray = BoolArrayType( Grid.Info.cellCount );
	markWhereFinestFluidIs( fluidMarkerArray, Grid.IJK, STLs, Grid.Info );
	//PLOT
	Grid.bouncebackMarkerArray = fluidMarkerArray;
	exportSectionCutPlotXY( Grid, kCut, 10*level+1 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	
	// BORDER
	BoolArrayType borderMarkerArray = BoolArrayType( Grid.Info.cellCount );
	borderMarkerArray.setValue( 0 );
	markDIADNeighbours( Neighbours, fluidMarkerArray, borderMarkerArray );
	BoolArrayType fluidInverseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	fluidInverseMarkerArray = fluidMarkerArray;
	invertBoolArray( fluidInverseMarkerArray );
	multiplyBoolArrays( fluidInverseMarkerArray, borderMarkerArray, borderMarkerArray );
	//PLOT
	Grid.bouncebackMarkerArray = borderMarkerArray;
	exportSectionCutPlotXY( Grid, kCut, 10*level+2 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	
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
		// getDIADEsotwistConnections( Grid ) <- NEED TO FINISH THIS!
		//PLOT
		Grid.bouncebackMarkerArray = Grid.bouncebackMarkerArray;
		exportSectionCutPlotXY( Grid, kCut, 10*level+3 );
		system("python3 ../include/plotter/plotter.py");
		//PLOT END
		//PLOT
		Grid.bouncebackMarkerArray.setValue( 0 );
		exportSectionCutPlotXY( Grid, kCut, 10*level+4 );
		system("python3 ../include/plotter/plotter.py");
		//PLOT END
		return;
	}
	
	// COARSE GRID BRANCH
	// THICK REFINEMENT REGION
	BoolArrayType refinementMarkerArray = BoolArrayType( Grid.Info.cellCount );
	refinementMarkerArray = borderMarkerArray;
	BoolArrayType newRefinementMarkerArray = BoolArrayType( Grid.Info.cellCount );
	newRefinementMarkerArray.setValue( 0 );
	for ( int layer = 0; layer < wallRefinementSpan; layer++ )
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
	const int cellCountFine = 8 * countMarkerCells( keepCellMarkerArray );
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
	// getDIADEsotwistConnections( Grid ) <- NEED TO FINISH THIS!
	
	//PLOT
	Grid.bouncebackMarkerArray = Grid.bouncebackMarkerArray;
	exportSectionCutPlotXY( Grid, kCut, 10*level+3 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
		
	//PLOT
	Grid.bouncebackMarkerArray = fineToCoarseMarkerArray;
	exportSectionCutPlotXY( Grid, kCut, 10*level+4 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	
	//PLOT
	Grid.bouncebackMarkerArray = coarseToFineMarkerArray;
	exportSectionCutPlotXY( Grid, kCut, 10*level+5 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	
	//PLOT
	Grid.bouncebackMarkerArray.setValue( 0 );
	exportSectionCutPlotXY( Grid, kCut, 10*level+6 );
	system("python3 ../include/plotter/plotter.py");
	//PLOT END
	
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
	
	/*
	for ( int cell = 0; cell < grids[0].Info.cellCount; cell++ )
	{
		int iPlus = Neighbours.iPlusArray.getElement(cell);
		int jPlus = Neighbours.jPlusArray.getElement(cell);
		int kPlus = Neighbours.kPlusArray.getElement(cell);
		int iMinus = Neighbours.iMinusArray.getElement(cell);
		int jMinus = Neighbours.jMinusArray.getElement(cell);
		int kMinus = Neighbours.kMinusArray.getElement(cell);		
		
		std::cout << "cell: " << cell <<
		" i: " << grids[0].IJK.iArray.getElement(cell) <<
		" j: " << grids[0].IJK.jArray.getElement(cell) <<
		" k: " << grids[0].IJK.kArray.getElement(cell) << std::endl;
		std::cout << "neighbours: " <<
		" iPlus: " << iPlus <<
		" jPlus: " << jPlus <<
		" kPlus: " << kPlus <<
		" iMinus: " << iMinus <<
		" jMinus: " << jMinus  <<
		" kMinus: " << kMinus  << std::endl;

		iPlus = std::max({0, iPlus});
		jPlus = std::max({0, jPlus});
		kPlus = std::max({0, kPlus});
		iMinus = std::max({0, iMinus});
		jMinus = std::max({0, jMinus});
		kMinus = std::max({0, kMinus});
		std::cout << "neighbour ijk: " << 
		" iPlus i: " << grids[0].IJK.iArray.getElement(iPlus) <<
		" iPlus j: " << grids[0].IJK.jArray.getElement(iPlus) <<
		" iPlus k: " << grids[0].IJK.kArray.getElement(iPlus) <<
		" jPlus i: " << grids[0].IJK.iArray.getElement(jPlus) <<
		" jPlus j: " << grids[0].IJK.jArray.getElement(jPlus) <<
		" jPlus k: " << grids[0].IJK.kArray.getElement(jPlus) <<
		" kPlus i: " << grids[0].IJK.iArray.getElement(kPlus) <<
		" kPlus j: " << grids[0].IJK.jArray.getElement(kPlus) <<
		" kPlus k: " << grids[0].IJK.kArray.getElement(kPlus) <<
		" iMinus i: " << grids[0].IJK.iArray.getElement(iMinus) <<
		" iMinus j: " << grids[0].IJK.jArray.getElement(iMinus) <<
		" iMinus k: " << grids[0].IJK.kArray.getElement(iMinus) <<
		" jMinus i: " << grids[0].IJK.iArray.getElement(jMinus) <<
		" jMinus j: " << grids[0].IJK.jArray.getElement(jMinus) <<
		" jMinus k: " << grids[0].IJK.kArray.getElement(jMinus) <<
		" kMinus i: " << grids[0].IJK.iArray.getElement(kMinus) <<
		" kMinus j: " << grids[0].IJK.jArray.getElement(kMinus) <<
		" kMinus k: " << grids[0].IJK.kArray.getElement(kMinus) <<
		std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

	}
	*/
}
