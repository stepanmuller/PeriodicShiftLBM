#include <TNL/Algorithms/sort.h>
#include "../include/STLFunctions.h"

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

// Receives sorted IJKSource that !must! be already sorted with key k, j, i so that k index changes the slowest
// For each Wanted, find a matching cell in Source. If its found, write its index to resultArray. If there is no such cell, write -2.		
void findMatchingIJKIndex( IJKArrayStruct &Wanted, IJKArrayStruct &Source, IntArrayType &resultArray )
{
	const int wantedCellCount = Wanted.iArray.getSize();
	const int sourceCellCount = Source.iArray.getSize();
	
	auto iWantedView = Wanted.iArray.getConstView();
	auto jWantedView = Wanted.jArray.getConstView();
	auto kWantedView = Wanted.kArray.getConstView();
	auto iSourceView = Source.iArray.getConstView();
	auto jSourceView = Source.jArray.getConstView();
	auto kSourceView = Source.kArray.getConstView();
	auto resultView = resultArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cellWanted ) mutable
	{
		if ( resultView[ cellWanted ] >= 0 ) return; // Do not overwrite the cell if its already valid from before
		
		const int iWanted = iWantedView[ cellWanted ];
		const int jWanted = jWantedView[ cellWanted ];
		const int kWanted = kWantedView[ cellWanted ];
		
		int start = 0;
		int end = sourceCellCount;
		int result = -2;
		
		// Search for k, j, and i in a single binary search
		while ( end > start )
		{
			int half = start + ( end - start ) / 2;
			int kHalf = kSourceView[ half ];
			int jHalf = jSourceView[ half ];
			int iHalf = iSourceView[ half ];
			
			if ( kHalf == kWanted && jHalf == jWanted && iHalf == iWanted )
			{
				result = half;
				break;
			}
			
			// Move the end bound if k is too big, OR if k matches but j is too big, OR if k and j match but i is too big
			if ( kHalf > kWanted || 
			   ( kHalf == kWanted && jHalf > jWanted ) || 
			   ( kHalf == kWanted && jHalf == jWanted && iHalf > iWanted ) ) 
			{
				end = half;
			}
			else 
			{
				start = half + 1;
			}
		}
		if ( result == -2 && start < sourceCellCount )
		{
			if ( kSourceView[ start ] == kWanted && jSourceView[ start ] == jWanted && iSourceView[ start ] == iWanted )
			{
				result = start;
			}
		}
		resultView[ cellWanted ] = result;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, wantedCellCount, cellLambda );
}

// Applies the bounceback marker onto an array using the getMarkers function
void ApplyMarkersFromFunction( BoolArrayType &markerArray, IJKArrayStruct &IJK, InfoStruct &Info )
{
	auto markerArrayView = markerArray.getView();
	
	auto iArrayView = IJK.iArray.getView();
	auto jArrayView = IJK.jArray.getView();
	auto kArrayView = IJK.kArray.getView();
	
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
				for ( int STLIndex = 0; STLIndex < (int)STLs.size(); STLIndex++ )
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

void getDIADNeighbours( IJKArrayStruct &IJK, std::vector<IntArrayType> &nbrArrays )
{
	const int cellCount = IJK.iArray.getSize();
	for ( int i = 0; i < 26; i++ ) { nbrArrays[i].setSize( cellCount ); nbrArrays[i].setValue( -2 ); }
	
	const int cxArray[27] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
	const int cyArray[27] = { 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
	const int czArray[27] = { 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

	IJKArrayStruct IJKShifted;
	
	for ( int direction = 1; direction < 27; direction++ )
	{
		IJKShifted = IJK;
		shiftIJK( IJKShifted, cxArray[direction], cyArray[direction], czArray[direction] );
		findMatchingIJKIndex( IJKShifted, IJK, nbrArrays[direction-1] );
	}
}

void getDIADEsotwistNbrArray( DIADGridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	IJKArrayStruct IJK = Grid.IJK;
	const int cellCount = Info.cellCount;
	Grid.EsotwistNbrArray.iNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.jNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.kNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.ijNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.ikNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.jkNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.ijkNbrArray.setSize( cellCount );
	Grid.EsotwistNbrArray.iNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.jNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.kNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.ijNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.ikNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.jkNbrArray.setValue( -2 );
	Grid.EsotwistNbrArray.ijkNbrArray.setValue( -2 );
	
	IJKArrayStruct IJKShifted;
	int invalidNeighbourCount;
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 0, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.iNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.iNbrArray );
	}

	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 1, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.jNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.jNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 0, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.kNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.kNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 1, 0, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.ijNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.ijNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 0, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.ikNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.ikNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 0, 1, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.jkNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.jkNbrArray );
	}
	
	IJKShifted = IJK;
	invalidNeighbourCount = 1;
	while ( invalidNeighbourCount > 0 )
	{
		shiftIJKPeriodic( IJKShifted, 1, 1, 1, Info );
		findMatchingIJKIndex( IJKShifted, IJK, Grid.EsotwistNbrArray.ijkNbrArray );
		invalidNeighbourCount = countInvalidIndexes( Grid.EsotwistNbrArray.ijkNbrArray );
	}
}

// Mark target cell as 1 if at least one of its neighbours in source is 1
void markDIADNeighbours( std::vector<IntArrayType> &nbrArrays, BoolArrayType &sourceArray, BoolArrayType &targetArray )
{
	const int cellCount = nbrArrays[0].getSize();
		
	auto sourceView = sourceArray.getConstView();
	auto targetView = targetArray.getView();
	
	for ( int i = 0; i < 26; i++ )
	{
		auto nbrView = nbrArrays[i].getConstView();
		auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
		{
			for ( int i = 0; i < 26; i++ ) 
			{
				if ( nbrView[cell] >=0 && sourceView[ nbrView[cell] ] ) 
				{ 
					targetView[cell] = true; 
					return; 
				}
			}
		};
		TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, cellCount, cellLambda );
	}
}

void buildDIADGrids( std::vector<DIADGridStruct> &grids, std::vector<STLStruct> STLs, const int level )
{
	DIADGridStruct &Grid = grids[level];
	std::cout << "Initial cell count on level " << level <<" : " << Grid.Info.cellCount << std::endl;

	sortIJK( Grid.IJK );
	std::vector<IntArrayType> nbrArrays( 26 );
	getDIADNeighbours( Grid.IJK, nbrArrays );	
		
	BoolArrayType fluidMarkerArray = BoolArrayType( Grid.Info.cellCount );
	markWhereFinestFluidIs( fluidMarkerArray, Grid.IJK, STLs, Grid.Info );
	
	// BORDER
	BoolArrayType borderMarkerArray = BoolArrayType( Grid.Info.cellCount );
	borderMarkerArray.setValue( 0 );
	markDIADNeighbours( nbrArrays, fluidMarkerArray, borderMarkerArray );
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
		getDIADEsotwistNbrArray( Grid );
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
		markDIADNeighbours( nbrArrays, refinementMarkerArray, newRefinementMarkerArray );
		multiplyBoolArrays( newRefinementMarkerArray, keepCellMarkerArray, refinementMarkerArray );
	}
	
	// FINE TO COARSE COMMUNICATION INTERFACE
	BoolArrayType fineToCoarseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	fineToCoarseMarkerArray.setValue( 0 );
	markDIADNeighbours( nbrArrays, refinementMarkerArray, fineToCoarseMarkerArray );
	multiplyBoolArrays( fluidMarkerArray, fineToCoarseMarkerArray, fineToCoarseMarkerArray );
	BoolArrayType refinementInverseMarkerArray = BoolArrayType( Grid.Info.cellCount );
	refinementInverseMarkerArray = refinementMarkerArray;
	invertBoolArray( refinementInverseMarkerArray );
	multiplyBoolArrays( fineToCoarseMarkerArray, refinementInverseMarkerArray, fineToCoarseMarkerArray );
	
	// COARSE TO FINE COMMUNICATION INTERFACE
	BoolArrayType coarseToFineMarkerArray = BoolArrayType( Grid.Info.cellCount );
	coarseToFineMarkerArray.setValue( 0 );
	markDIADNeighbours( nbrArrays, fineToCoarseMarkerArray, coarseToFineMarkerArray );
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
	std::cout << "Final cell count on level " << level <<" : " << Grid.Info.cellCount << std::endl;
	
	// FINAL BOUNCEBACK PASS IN CASE WE ARE NOT REFINING IN SOME WALL AREA
	Grid.bouncebackMarkerArray.setSize(Grid.Info.cellCount);
	Grid.bouncebackMarkerArray.setValue( 0 );
	BoolArrayType markerArraySTL = BoolArrayType( Grid.Info.cellCount );	
	ApplyMarkersFromFunction( Grid.bouncebackMarkerArray, Grid.IJK, Grid.Info );
	for ( int STLIndex = 0; STLIndex < (int)STLs.size(); STLIndex++ )
	{
		const bool insideMarkerValue = 1;
		ApplyMarkersInsideSTL( markerArraySTL, Grid.IJK, STLs[STLIndex], insideMarkerValue, Grid.Info );
		sumBoolArrays( Grid.bouncebackMarkerArray, markerArraySTL, Grid.bouncebackMarkerArray );
	}
	
	// BUILD ESOTWIST CONNECTIONS
	getDIADEsotwistNbrArray( Grid ); // <- NEED TO FINISH THIS!
	
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
	buildDIADGrids( grids, STLs, level + 1 );
	
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
	findMatchingIJKIndex( fineIJK, grids[level+1].IJK, Grid.coarseToFineWriteArray );
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
	findMatchingIJKIndex( fineIJK, grids[level+1].IJK, Grid.fineToCoarseReadArray );
	const int invalidConnections2 = countInvalidIndexes( Grid.fineToCoarseReadArray );
	if ( invalidConnections2 > 0 )
	{
		throw std::runtime_error( 
			"DIAD Grid Generation Error: " + std::to_string(invalidConnections2) + 
			" invalid coarse-to-fine connections found on level " + std::to_string(level) 
		);
	}
}
