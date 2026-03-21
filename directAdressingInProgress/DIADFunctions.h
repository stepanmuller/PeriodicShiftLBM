
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
		if ( kSourceView[kStepView[start]] == kWanted )
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
	return;
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
