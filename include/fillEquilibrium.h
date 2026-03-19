void fillEquilibriumFromFunction( GridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
		float feq[27];
		float rho = 1.f;
		float ux, uy, uz = 0.f;
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		getInitialRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Marker, Info );
		getFeq(rho, ux, uy, uz, feq);
		for ( int i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Grid.Info.cellCount, cellLambda );
}
