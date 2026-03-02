// Version with call for bouncebackMarker as an explicit function of i, j, k
void applyLocalCellUpdate( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	InfoStruct Info = Grid.Info;
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		
		MarkerStruct Marker;
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( Marker.ghost ) return;
		if ( (iCell>Info.iSubgridStart&&iCell<Info.iSubgridEnd-1) 
			&& (jCell>Info.jSubgridStart&&jCell<Info.jSubgridEnd-1) 
			&& (kCell>Info.kSubgridStart&&kCell<Info.kSubgridEnd-1) 
			) return;
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		
		if ( Marker.bounceback )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( Marker.fluid )
			{
				// do nothing, just skip the else block below
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info );  
				if ( Marker.periodic ) outerNormalZ = 0;
				getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Info );
				if ( Marker.mirror )
				{
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
				}
				else if ( Marker.givenRho && !Marker.givenUxUyUz )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !Marker.givenRho && Marker.givenUxUyUz )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !Marker.givenRho && !Marker.givenUxUyUz )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Info );
			applyCollision( f, Info.nu, SmagorinskyConstant );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
}

// Version with bouncebackMarker loaded from memory
void applyLocalCellUpdateBB( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	auto bouncebackArrayView = Grid.bouncebackMarkerArray.getConstView();
	InfoStruct Info = Grid.Info;
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		
		MarkerStruct Marker;
		Marker.bounceback = bouncebackArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( Marker.ghost ) return;
		if ( (iCell>Info.iSubgridStart&&iCell<Info.iSubgridEnd-1) 
			&& (jCell>Info.jSubgridStart&&jCell<Info.jSubgridEnd-1) 
			&& (kCell>Info.kSubgridStart&&kCell<Info.kSubgridEnd-1) 
			) return;
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		
		if ( Marker.bounceback )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( Marker.fluid )
			{
				// do nothing, just skip the else block below
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info );  
				if ( Marker.periodic ) outerNormalZ = 0;
				getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Info );
				if ( Marker.mirror )
				{
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
				}
				else if ( Marker.givenRho && !Marker.givenUxUyUz )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !Marker.givenRho && Marker.givenUxUyUz )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !Marker.givenRho && !Marker.givenUxUyUz )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Info );
			applyCollision( f, Info.nu, SmagorinskyConstant );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
}
