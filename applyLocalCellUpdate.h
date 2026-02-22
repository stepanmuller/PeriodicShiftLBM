// Version with call for bouncebackMarker as an explicit function of i, j, k without report array
void applyLocalCellUpdate( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Grid );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Grid );
		
		bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
		getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Grid );
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		
		if ( bouncebackMarker )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker )
			{
				// do nothing, just skip the else block below
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, periodicMarker, outerNormalX, outerNormalY, outerNormalZ, Grid );
				getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Grid );
				if ( mirrorMarker )
				{
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
				}
				else if ( givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && givenUxUyUzMarker )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Grid );
			applyCollision( f, Grid.nu, SmagorinskyConstant );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Grid.cellCount, cellLambda );
}

// Version with bouncebackMarker loaded from memory without report array
void applyLocalCellUpdate( FStruct &F, BoolArrayType &bouncebackArray, GridStruct &Grid )
{
	auto fArrayView  = F.fArray.getView();
	auto shifterView  = F.shifter.getConstView();
	auto bouncebackArrayView = bouncebackArray.getConstView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Grid );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Grid );
		
		bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
		bouncebackMarker = bouncebackArrayView( cell );
		getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Grid );
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		
		if ( bouncebackMarker )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker )
			{
				// do nothing, just skip the else block below
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, periodicMarker, outerNormalX, outerNormalY, outerNormalZ, Grid );
				getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Grid );
				if ( mirrorMarker )
				{
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
				}
				else if ( givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && givenUxUyUzMarker )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Grid );
			applyCollision( f, SmagorinskyConstant );
		}

		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Grid.cellCount, cellLambda );
}
