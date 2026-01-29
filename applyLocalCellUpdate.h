void applyLocalCellUpdate( MarkerStruct& Marker, FloatArray4DType fArray, InfoStruct &Info )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getConstView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getConstView();
	auto givenRhoMarkerArrayView = Marker.givenRhoArray.getConstView();
	auto givenUxUyUzMarkerArrayView = Marker.givenUxUyUzArray.getConstView();
	
	auto fArrayView  = fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const TripleIndexType& tripleIndex ) mutable
	{
		const int iCell = tripleIndex.x();
		const int jCell = tripleIndex.y();
		const int kCell = tripleIndex.z();
		int iStreamed[27], jStreamed[27], kStreamed[27];
		getStreamedIndexes( iCell, jCell, kCell, iStreamed, jStreamed, kStreamed, Info );
		
		bool fluidMarker = fluidMarkerArrayView( iCell, jCell, kCell );
		bool bouncebackMarker = bouncebackMarkerArrayView( iCell, jCell, kCell );
		bool givenRhoMarker = givenRhoMarkerArrayView( iCell, jCell, kCell );
		bool givenUxUyUzMarker = givenUxUyUzMarkerArrayView( iCell, jCell, kCell );
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ ) f[direction] = fArrayView( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] );
		
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
				if ( givenRhoMarker && !givenUxUyUzMarker )
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
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] ) = f[direction];
	};
	TripleIndexType start{ 0, 0, 0 };
	TripleIndexType end{ Info.cellCountX, Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( start, end, cellLambda );
}
