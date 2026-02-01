void applyLocalCellUpdate( FStruct &F, InfoStruct &Info )
{
	auto fArrayView  = F.fArray.getView();
	//auto shifterView  = F.shifter.getConstView();
	
	int shifter[27];
	for ( int direction = 0; direction < 27; direction++ ) shifter[direction] = F.shifter.getElement( direction );
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
			
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifter, Info );
		
		bool fluidMarker, bouncebackMarker, mirrorMarker, givenRhoMarker, givenUxUyUzMarker;
		getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, givenRhoMarker, givenUxUyUzMarker, Info );
		
		float f[27];
		float rho, ux, uy, uz;
		
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, shiftedIndex[direction]);
		//for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(direction, cell);
		
		if ( bouncebackMarker )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker )
			{
				getRhoUxUyUz( rho, ux, uy, uz, f );
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info );
				rho = 1.f;
				ux = 0.f;
				uy = 0.f;
				uz = uzInlet;
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
			applyCollision( rho, ux, uy, uz, f );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = f[direction];
		//for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, cell ) = f[direction];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
}

/*
void applyLocalCellUpdate( 	MarkerStruct& Marker, FStruct& F, InfoStruct& Info )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getConstView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getConstView();
	auto givenRhoMarkerArrayView = Marker.givenRhoArray.getConstView();
	auto givenUxUyUzMarkerArrayView = Marker.givenUxUyUzArray.getConstView();
	
	auto shifterView = F.shifter.getConstView();
	auto fArrayView  = F.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount.n) { shiftedIndex[i] -= cellCount.n; }
		}
		bool fluidMarker = fluidMarkerArrayView[cell];
		bool bouncebackMarker = bouncebackMarkerArrayView[cell];
		bool givenRhoMarker = givenRhoMarkerArrayView[cell];
		bool givenUxUyUzMarker = givenUxUyUzMarkerArrayView[cell];
		
		float f[27];
		float rho, ux, uy, uz;
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		
		if ( bouncebackMarker == 1 )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker == 1 )
			{
				getRhoUxUyUz(rho, ux, uy, uz, f);
			}
			else
			{
				short outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal(cell, outerNormalX, outerNormalY, outerNormalZ, cellCount);
				rho = 1.f;
				ux = 0.f;
				uy = 0.f;
				uz = uzInlet;
				if ( givenRhoMarker == 1 && givenUxUyUzMarker == 0 )
				{
					restoreUxUyUz(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				else if ( givenRhoMarker == 0 && givenUxUyUzMarker == 1 )
				{
					restoreRho(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				else if ( givenRhoMarker == 0 && givenUxUyUzMarker == 0 )
				{
					restoreRhoUxUyUz(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				applyMBBC(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			}
			applyCollision(rho, ux, uy, uz, f);
		}
		for (size_t i = 0; i < 27; i++)	fArrayView(i, shiftedIndex[i]) = f[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount.n, cellLambda );
}
*/
