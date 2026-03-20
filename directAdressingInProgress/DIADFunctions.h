// DIAD stands for directly adressed

// This holds cell indexes of the neighbour cells in the main positive and negative i,j,k directions for finding geometric neighbours.
struct DIADNeighboursStruct { 	IntArrayType iPlusArray; IntArrayType jPlusArray; IntArrayType kPlusArray; 
								IntArrayType iMinusArray; IntArrayType jMinusArray; IntArrayType kMinusArray; };

// This holds cell indexes of the neighbour cells in positive i,j,k directions for Esotwist streaming.
struct DIADEsotwistConnectionsStruct { 	IntArrayType iNbrArray; IntArrayType jNbrArray; IntArrayType kNbrArray; 
										IntArrayType ijNbrArray; IntArrayType ikNbrArray; IntArrayType jkNbrArray; 
										IntArrayType ijkNbrArray; }; 

struct DIADGridStruct { InfoStruct Info; FloatArray2DType fArray; BoolArrayType bouncebackMarkerArray; 
						IntArrayType iArray; IntArrayType jArray; IntArrayType kArray; DIADEsotwistConnectionsStruct EsotwistConnections; 
						IntArrayType fineToCoarseWriteArray; IntArrayType fineToCoarseReadArray; 
						IntArrayType coarseToFineWriteArray; IntArrayType coarseToFineReadArray; }; 

// Receives sorted iSource, jSource, kSource that !must! be already sorted with key k, j, i -> k index changes the slowest
// For each Wanted, find a matching cell in Source. If its found, write its index to cellSourceIndex. If there is no such cell, write -2.		
void findMatchingIJKIndex( 	IntArrayType &iWantedArray, IntArrayType &jWantedArray, IntArrayType &kWantedArray, 
							IntArrayType &iSourceArray, IntArrayType &jSourceArray, IntArrayType &kSourceArray, 
							IntArrayType &cellSourceIndexArray )
{
	const int wantedCellCount = iWantedArray.getSize();
	const int sourceCellCount = iSourceArray.getSize();
	
	IntArrayTypeCPU kSourceArrayCPU;
	kSourceArrayCPU = kSourceArray;
	IntArrayTypeCPU kChangesArrayCPU = IntArrayTypeCPU( sourceCellCount );
	int counter = 0;
	int lastIndex = -2;
	for ( int cell = 0; cell < sourceCellCount; cell++ )
	{
		int currentIndex = kSourceArrayCPU[ cell ];
		if ( currentIndex != lastIndex )
		{
			kChangesArrayCPU[ counter ] = cell;
			lastIndex = currentIndex;
			counter++;
		}
	}
	const int kChangeCount = counter;
	kChangesArrayCPU.resize( kChangeCount );
	IntArrayType kChangesArray;
	kChangesArray = kChangesArrayCPU;
	
	auto iWantedView = iWantedArray.getConstView();
	auto jWantedView = jWantedArray.getConstView();
	auto kWantedView = kWantedArray.getConstView();
	auto iSourceView = iSourceArray.getConstView();
	auto jSourceView = jSourceArray.getConstView();
	auto kSourceView = kSourceArray.getConstView();
	auto cellSourceIndexView = cellSourceIndexArray.getView();
	auto kChangesView = kChangesArray.getConstView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cellWanted ) mutable
	{
		const int iWanted = iWantedView[ cellWanted ];
		const int jWanted = jWantedView[ cellWanted ];
		const int kWanted = kWantedView[ cellWanted ];
		int start = 0;
		int end = kChangeCount;
		int result = -2;
		// Now by using the kIndexChanges array, we will reduce the search interval
		while ( end > start )
		{
			int half = start + ( end - start ) / 2;
			if ( kSourceView[kChangesView[half]] == kWanted )
			{
				result = half;
				break;
			}
			if ( kSourceView[kChangesView[half]] > kWanted ) end = half;
			else start = half + 1;
		}
		if ( kSourceView[kChangesView[start]] == kWanted )
		{
			result = start;
		}
		if ( result == -2 )
		{
			cellSourceIndexView[ cellWanted ] = -2;
			return;
		}
		
		start = kChangesView[ result ];
		end = sourceCellCount;
		if ( result < kChangeCount - 1 ) end = kChangesView[ result + 1 ];
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
		
		cellSourceIndexView[ cellWanted ] = result;
	};
	
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, wantedCellCount, cellLambda );
	return;
}

void ApplyMarkersInsideSTL( BoolArrayType &markerArray, IntArrayType &iArray, IntArrayType &jArray, IntArrayType &kArray, STLStruct &STL, const bool &insideMarkerValue, InfoStruct &Info )
{
	std::cout << "DIAD Applying markers inside STL" << std::endl;
	auto markerArrayView = markerArray.getView();
	
	auto iArrayView = iArray.getView();
	auto jArrayView = jArray.getView();
	auto kArrayView = kArray.getView();

	auto axArrayView = STL.axArray.getConstView();
	auto ayArrayView = STL.ayArray.getConstView();
	auto azArrayView = STL.azArray.getConstView();
	auto bxArrayView = STL.bxArray.getConstView();
	auto byArrayView = STL.byArray.getConstView();
	auto bzArrayView = STL.bzArray.getConstView();
	auto cxArrayView = STL.cxArray.getConstView();
	auto cyArrayView = STL.cyArray.getConstView();
	auto czArrayView = STL.czArray.getConstView();
	
	IntArray2DType intersectionCounterArray;
	intersectionCounterArray.setSizes( Info.cellCountX, Info.cellCountY );
	intersectionCounterArray.setValue( 0 );
	auto intersectionCounterArrayView = intersectionCounterArray.getView();

    auto counterLambda = [ = ] __cuda_callable__( const int triangleIndex ) mutable
    {
		// transform into the coordinate system of the LBM grid
		const float ax = axArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float ay = ayArrayView[ triangleIndex ] + STL.oy - Info.oy;
		const float bx = bxArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float by = byArrayView[ triangleIndex ] + STL.oy - Info.oy;
		const float cx = cxArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float cy = cyArrayView[ triangleIndex ] + STL.oy - Info.oy;
		// transform STL floats to integer grid that is 100x finer than the LBM grid to prevent float errors
		// make the STL coords odd, rays will be even, this prevents hitting a vortex
		const float scale = 50.0f / Info.res;
		const long long ak = (long long)(round( ax * scale )) * 2 + 1;
		const long long al = (long long)(round( ay * scale )) * 2 + 1;
		const long long bk = (long long)(round( bx * scale )) * 2 + 1;
		const long long bl = (long long)(round( by * scale )) * 2 + 1;
		const long long ck = (long long)(round( cx * scale )) * 2 + 1;
		const long long cl = (long long)(round( cy * scale )) * 2 + 1;
		
		const long long kmin = std::max({ 0LL, std::min({ ak, bk, ck, (long long)(Info.cellCountX-1)*100 }) });
		const long long kmax = std::min({ (long long)(Info.cellCountX-1)*100, std::max({ ak, bk, ck, 0LL }) });
		const long long lmin = std::max({ 0LL, std::min({ al, bl, cl, (long long)(Info.cellCountY-1)*100 }) });
		const long long lmax = std::min({ (long long)(Info.cellCountY-1)*100, std::max({ al, bl, cl, 0LL }) });
		
		const long long imin = (kmin + 99) / 100;
		const long long imax = kmax / 100;
		const long long jmin = (lmin + 99) / 100;
		const long long jmax = lmax / 100;
		
		for ( int j = jmin; j <= jmax; j++ )
		{
			for ( int i = imin; i <= imax; i++ )
			{
				const long long rayK = i * 100;
				const long long rayL = j * 100;
				// transform the triangle into coordinate system where ray is [0, 0]
				const long long ak0 = ak - rayK;
				const long long al0 = al - rayL;
				const long long bk0 = bk - rayK;
				const long long bl0 = bl - rayL;
				const long long ck0 = ck - rayK;
				const long long cl0 = cl - rayL;

				const bool rayHit = getRayHitYesNo( ak0, al0, bk0, bl0, ck0, cl0 );

				if (rayHit) 
				{
					TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(intersectionCounterArrayView(i, j), 1);
				}
			}
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STL.triangleCount, counterLambda );
	
	auto fetch = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int i = singleIndex % Info.cellCountX;
		const int j = singleIndex / Info.cellCountX;
		return intersectionCounterArrayView( i, j );
	};
	auto reduction = [] __cuda_callable__( const int& a, const int& b )
	{
		return TNL::max( a, b );
	};
	const int start = 0;
	const int end = Info.cellCountX * Info.cellCountY;
	int intersectionCountMax = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0 );
	std::cout << "	intersectionCountMax: " << intersectionCountMax << std::endl; 
	
	if ( intersectionCountMax < 1 ) intersectionCountMax = 1;
	// extracting the intersection indexes
	IntArray3DType intersectionIndexArray;
	intersectionIndexArray.setSizes( intersectionCountMax, Info.cellCountX, Info.cellCountY );
	intersectionIndexArray.setValue( Info.cellCountZ );
	auto intersectionIndexArrayView = intersectionIndexArray.getView();
	intersectionCounterArray.setValue(0);
	
	auto rayHitIndexLambda = [ = ] __cuda_callable__( const int triangleIndex ) mutable
    {
		// transform into the coordinate system of the LBM grid
		const float ax = axArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float ay = ayArrayView[ triangleIndex ] + STL.oy - Info.oy;
		const float az = azArrayView[ triangleIndex ] + STL.oz - Info.oz;
		const float bx = bxArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float by = byArrayView[ triangleIndex ] + STL.oy - Info.oy;
		const float bz = bzArrayView[ triangleIndex ] + STL.oz - Info.oz;
		const float cx = cxArrayView[ triangleIndex ] + STL.ox - Info.ox;
		const float cy = cyArrayView[ triangleIndex ] + STL.oy - Info.oy;
		const float cz = czArrayView[ triangleIndex ] + STL.oz - Info.oz;
		// transform STL floats to integer grid that is 100x finer than the LBM grid to prevent float errors
		// make the STL coords odd, rays will be even, this prevents hitting a vortex
		const float scale = 50.0f / Info.res;
		const long long ak = (long long)(round( ax * scale )) * 2 + 1;
		const long long al = (long long)(round( ay * scale )) * 2 + 1;
		const long long bk = (long long)(round( bx * scale )) * 2 + 1;
		const long long bl = (long long)(round( by * scale )) * 2 + 1;
		const long long ck = (long long)(round( cx * scale )) * 2 + 1;
		const long long cl = (long long)(round( cy * scale )) * 2 + 1;
		
		const long long kmin = std::max({ 0LL, std::min({ ak, bk, ck, (long long)(Info.cellCountX-1)*100 }) });
		const long long kmax = std::min({ (long long)(Info.cellCountX-1)*100, std::max({ ak, bk, ck, 0LL }) });
		const long long lmin = std::max({ 0LL, std::min({ al, bl, cl, (long long)(Info.cellCountY-1)*100 }) });
		const long long lmax = std::min({ (long long)(Info.cellCountY-1)*100, std::max({ al, bl, cl, 0LL }) });
		
		const long long imin = (kmin + 99) / 100;
		const long long imax = kmax / 100;
		const long long jmin = (lmin + 99) / 100;
		const long long jmax = lmax / 100;
		
		for ( int j = jmin; j <= jmax; j++ )
		{
			for ( int i = imin; i <= imax; i++ )
			{
				const long long rayK = i * 100;
				const long long rayL = j * 100;
				// transform the triangle into coordinate system where ray is [0, 0]
				const long long ak0 = ak - rayK;
				const long long al0 = al - rayL;
				const long long bk0 = bk - rayK;
				const long long bl0 = bl - rayL;
				const long long ck0 = ck - rayK;
				const long long cl0 = cl - rayL;

				const bool rayHit = getRayHitYesNo( ak0, al0, bk0, bl0, ck0, cl0 );

				if (rayHit) 
				{
					const int writePosition = TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(intersectionCounterArrayView(i, j), 1);
					const float rayX = i * Info.res;
					const float rayY = j * Info.res;
					const float rayZ = getRayHitZCoordinate( ax, ay, az, bx, by, bz, cx, cy, cz, rayX, rayY );
					int k = (int)std::max( 0.0f, ceilf(rayZ / Info.res) );
					k = std::min( k, Info.cellCountZ );
					intersectionIndexArrayView(writePosition, i, j) = k;
				}
			}
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STL.triangleCount, rayHitIndexLambda );
	
	auto rayLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int iCell = doubleIndex.x();
		const int jCell = doubleIndex.y();
		
		for ( int layer = 1; layer < intersectionCountMax; layer++ ) // sort
		{
			int key = intersectionIndexArrayView( layer, iCell, jCell );
			int slider = layer - 1;
			while ( slider >= 0 && intersectionIndexArrayView( slider, iCell, jCell ) > key ) 
			{
				intersectionIndexArrayView( slider + 1, iCell, jCell ) = intersectionIndexArrayView( slider, iCell, jCell );
				slider = slider - 1;
			}
			intersectionIndexArrayView( slider + 1, iCell, jCell ) = key;
		}
	};
	IntPairType startList{ 0, 0 };
	IntPairType endList{ Info.cellCountX, Info.cellCountY };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(startList, endList, rayLambda );	
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		const int iCell = iArrayView[ cell ];
		const int jCell = jArrayView[ cell ];
		const int kCell = kArrayView[ cell ];
		
		bool markerValue = !insideMarkerValue; // we are always starting outside
		for ( int interval = 0; interval <= intersectionCountMax; interval++ )
		{
			int start = 0;
			int end = 0;
			if ( interval == 0 ) end = intersectionIndexArrayView( 0, iCell, jCell );
			else if ( interval == intersectionCountMax ) 
			{
				start = intersectionIndexArrayView( intersectionCountMax-1, iCell, jCell );
				end = Info.cellCountZ;
			}
			else
			{
				start = intersectionIndexArrayView( interval-1, iCell, jCell );
				end = intersectionIndexArrayView( interval, iCell, jCell );
			}
			if ( (kCell >= start) && (kCell < end) ) 
			{
				markerArrayView( cell ) = markerValue;
				return;
			}
			markerValue = !markerValue;
		}		
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );	
	
	std::cout << "	Markers inside STL applied" << std::endl;
}

void ApplyMarkersFromFunction( BoolArrayType &markerArray, IntArrayType &iArray, IntArrayType &jArray, IntArrayType &kArray, InfoStruct &Info )
{
	std::cout << "DIAD Applying markers from function" << std::endl;
	auto markerArrayView = markerArray.getView();
	
	auto iArrayView = iArray.getView();
	auto jArrayView = jArray.getView();
	auto kArrayView = kArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		const int iCell = iArrayView[ cell ];
		const int jCell = jArrayView[ cell ];
		const int kCell = kArrayView[ cell ];
		MarkerStruct Marker;
		getMarkers( iCell, jCell, kCell, Marker, Info );
		markerArrayView( cell ) = Marker.bounceback;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );	
}

// Do a repetetive fluid marking in order to find which of our coarse cells could possibly contain at least one finest (maximum refinement level) fluid cell. 
// To do this we will be repeatedly temporarily shifting the origin of our grid, to simulate being on a finer grid.
void markWhereFinestFluidIs( BoolArrayType &fluidMarkerArray, IntArrayType &iArray, IntArrayType &jArray, IntArrayType &kArray, std::vector<STLStruct> STLs, InfoStruct &Info )
{
	const int shiftCount = std::pow(2, (gridLevelCount - Info.gridID - 1));
	const float finestRes = Info.res / ( std::pow(2, (gridLevelCount - Info.gridID - 1)) );
	const float shiftStart = - 0.5f * Info.res + 0.5f * finestRes;
	const float oxOriginal = Info.ox;
	const float oyOriginal = Info.oy;
	const float ozOriginal = Info.oz;
	BoolArrayType markerArray = BoolArrayType( Info.cellCount );
	BoolArrayType markerArraySTL = BoolArrayType( Info.cellCount );
	for ( int shiftCountX = 0; shiftCountX < shiftCount; shiftCountX++ )
	{
		for ( int shiftCountY = 0; shiftCountY < shiftCount; shiftCountY++ )
		{
			for ( int shiftCountZ = 0; shiftCountZ < shiftCount; shiftCountZ++ )
			{
				Info.ox = oxOriginal + shiftStart + shiftCountX * finestRes;
				Info.oy = oyOriginal + shiftStart + shiftCountY * finestRes;
				Info.oz = ozOriginal + shiftStart + shiftCountZ * finestRes;
				ApplyMarkersFromFunction( markerArray, iArray, jArray, kArray, Info );
				for ( int STLIndex = 0; STLIndex < STLs.size(); STLIndex++ )
				{
					const bool insideMarkerValue = 1;
					ApplyMarkersInsideSTL( markerArraySTL, iArray, jArray, kArray, STLs[STLIndex], insideMarkerValue, Info );
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
