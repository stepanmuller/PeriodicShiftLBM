void applyLocalCellUpdate( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( Marker.ghost ) return;
		
		if ( (iCell>Info.iSubgridStart+1&&iCell<Info.iSubgridEnd-2) 
			&& (jCell>Info.jSubgridStart+1&&jCell<Info.jSubgridEnd-2) 
			&& (kCell>Info.kSubgridStart+1&&kCell<Info.kSubgridEnd-2) 
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
				if ( Marker.periodicX ) outerNormalX = 0;
				if ( Marker.periodicY ) outerNormalY = 0;
				if ( Marker.periodicZ ) outerNormalZ = 0;
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

/// VERSION WITH ADDED SCALAR TRANSPORT ////
/// Raoyang Zhang, Hongli Fan, Hudong Chen - A lattice Boltzmann approach for solving scalar transport equations, 2011
void applyLocalCellUpdate( GridStruct &Grid, ScalarTransportStruct &ScalarTransport )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	
	auto TArrayView = ScalarTransport.TArray.getView();
	const float tauT = ScalarTransport.tauT;
	
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		ScalarTransportMarkerStruct ScalarTransportMarker;
		getScalarTransportMarkers( iCell, jCell, kCell, ScalarTransportMarker, Info );
		
		if ( Marker.ghost ) return;
		
		if ( (iCell>Info.iSubgridStart+1&&iCell<Info.iSubgridEnd-2) 
			&& (jCell>Info.jSubgridStart+1&&jCell<Info.jSubgridEnd-2) 
			&& (kCell>Info.kSubgridStart+1&&kCell<Info.kSubgridEnd-2) 
			) return;
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		
		float T[27];
		for ( int direction = 0; direction < 27; direction++ )	T[direction] = TArrayView(direction, shiftedIndex[direction]);
		
		if ( Marker.bounceback )
		{
			applyBounceback(f);
			applyBounceback(T);
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
				if ( Marker.periodicX ) outerNormalX = 0;
				if ( Marker.periodicY ) outerNormalY = 0;
				if ( Marker.periodicZ ) outerNormalZ = 0;
				getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Info );
				if ( Marker.mirror )
				{
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
					applyMirror( outerNormalX, outerNormalY, outerNormalZ, T );
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
			if ( ScalarTransportMarker.givenT )
			{
				float givenT = 0.f;
				getGivenT( iCell, jCell, kCell, givenT, Info );
				applyScalarTransportGivenT( T, givenT );
			}
			const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Info );
			applyCollision( f, Info.nu, SmagorinskyConstant );
			applyScalarTransportCollision( f, T, tauT );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
		
		for ( int direction = 0; direction < 27; direction++ ) TArrayView( direction, shiftedIndex[direction] ) = T[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
}
