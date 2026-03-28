// DIAD Version that dynamically downsamples finest grids to fit VRAM
void getFlowReportGeneral( std::vector<DIADGridStruct> &grids, const int &cutIndex, XYZBoundsStruct &Bounds, FlowReportStruct &FlowReport, PlaneEnum plane )
{
	const int levelCount = grids.size();
	
	// 1. Find the finest level that fits within the memory limit
	int targetLevelCount = levelCount;
	int targetCellCountHorizontal = 0, targetCellCountVertical = 0;
	
	while ( targetLevelCount > 1 )
	{
		InfoStruct Info = grids[targetLevelCount - 1].Info;
		if ( plane == XY ) 		{ targetCellCountHorizontal = Info.cellCountX; targetCellCountVertical = Info.cellCountY; }
		else if ( plane == ZY ) { targetCellCountHorizontal = Info.cellCountZ; targetCellCountVertical = Info.cellCountY; }
		else 					{ targetCellCountHorizontal = Info.cellCountZ; targetCellCountVertical = Info.cellCountX; }
		
		// Use long long to prevent integer overflow on massive grids
		long long dataSize = (long long)targetCellCountHorizontal * targetCellCountVertical;
		if ( dataSize < 20000000 ) break;
		
		targetLevelCount--;
	}
	
	// How much smaller the output array is compared to the absolute finest grid
	// Using bitshift (1 << X) which is mathematically identical to std::pow(2, X) but faster/safer
	const int targetScale = 1 << (levelCount - targetLevelCount); 

	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	SectionCut.uxArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	SectionCut.uyArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	SectionCut.uzArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	SectionCut.markerArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	SectionCut.gridIDArray.setSizes( targetCellCountVertical, targetCellCountHorizontal );
	
	SectionCut.rhoArray.setValue( 1.f );
	SectionCut.uxArray.setValue( 0.f );
	SectionCut.uyArray.setValue( 0.f );
	SectionCut.uzArray.setValue( 0.f );
	SectionCut.markerArray.setValue( 1 );
	SectionCut.gridIDArray.setValue( 0 );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();
	auto gridIDArrayView = SectionCut.gridIDArray.getView();
	
	// 2. Loop through ALL grids (none are dropped)
	for ( int level = 0; level < levelCount; level++ )
	{
		DIADGridStruct &Grid = grids[level];
		InfoStruct Info = Grid.Info;
		
		// cellScale is relative to the absolute finest grid (levelCount)
		const int cellScale = std::pow(2, (levelCount - Info.gridID - 1) );
		
		auto fArrayView  = Grid.fArray.getConstView();
		bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
		auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
		
		auto iView = Grid.IJK.iArray.getConstView();
		auto jView = Grid.IJK.jArray.getConstView();
		auto kView = Grid.IJK.kArray.getConstView();
		
		bool esotwistFlipper = Grid.esotwistFlipper;
		auto iNbrView = Grid.EsotwistNbrArray.iNbrArray.getConstView();
		auto jNbrView = Grid.EsotwistNbrArray.jNbrArray.getConstView();
		auto kNbrView = Grid.EsotwistNbrArray.kNbrArray.getConstView();
		auto ijNbrView = Grid.EsotwistNbrArray.ijNbrArray.getConstView();
		auto ikNbrView = Grid.EsotwistNbrArray.ikNbrArray.getConstView();
		auto jkNbrView = Grid.EsotwistNbrArray.jkNbrArray.getConstView();
		auto ijkNbrView = Grid.EsotwistNbrArray.ijkNbrArray.getConstView();

		auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
		{
			// Coordinates on the absolute finest grid
			int iCell = iView[ cell ] * cellScale; 
			int jCell = jView[ cell ] * cellScale;
			int kCell = kView[ cell ] * cellScale;
			
			int indexHorizontal = 0;
			int indexVertical = 0;
			
			// cutIndex remains fully accurate because it checks against the absolute finest coordinates
			if ( plane == XY ) 
			{
				if ( cutIndex < kCell || cutIndex >= kCell + cellScale ) return;
				indexHorizontal = iCell; 
				indexVertical = jCell; 
			}
			else if ( plane == ZY ) 
			{ 
				if ( cutIndex < iCell || cutIndex >= iCell + cellScale ) return; 
				indexVertical = jCell; 
				indexHorizontal = kCell; 
			}
			else // ZX plane
			{ 
				if ( cutIndex < jCell || cutIndex >= jCell + cellScale ) return; 
				indexVertical = iCell; 
				indexHorizontal = kCell; 
			}
			
			DIADEsotwistNbrStruct Nbr;
			Nbr.i = iNbrView( cell );
			Nbr.j = jNbrView( cell );
			Nbr.k = kNbrView( cell );
			Nbr.ij = ijNbrView( cell );
			Nbr.ik = ikNbrView( cell );
			Nbr.jk = jkNbrView( cell );
			Nbr.ijk = ijkNbrView( cell );
			
			float f[27];
			int cellReadIndex[27];
			int fReadIndex[27];
			getEsotwistWriteIndex( cell, cellReadIndex, fReadIndex, Nbr, esotwistFlipper, Info ); 
			for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fReadIndex[direction], cellReadIndex[direction]);
			
			float rho, ux, uy, uz;
			getRhoUxUyUz(rho, ux, uy, uz, f);

			MarkerStruct Marker;
			if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
			const float marker = Marker.bounceback;
			
			// 3. Mapping coordinates to the scaled-down output array
			int outYStart = indexVertical / targetScale;
			int outXStart = indexHorizontal / targetScale;
			
			// How many pixels this cell spans on the output array (minimum 1 for downsampled fine grids)
			int spanY = max(1, cellScale / targetScale);
			int spanX = max(1, cellScale / targetScale);
			
			for ( int shiftY = 0; shiftY < spanY; shiftY++ )
			{
				int y = outYStart + shiftY;
				if (y >= targetCellCountVertical) continue; // Memory safety bound
				
				for ( int shiftX = 0; shiftX < spanX; shiftX++ )
				{
					int x = outXStart + shiftX;
					if (x >= targetCellCountHorizontal) continue; // Memory safety bound
					
					rhoArrayView( y, x ) = rho;
					uxArrayView( y, x ) = ux;
					uyArrayView( y, x ) = uy;
					uzArrayView( y, x ) = uz;
					markerArrayView( y, x ) = marker;
					gridIDArrayView( y, x ) = Info.gridID;
				}
			}
		};
		TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
	}
}

void getFlowReportXY( std::vector<DIADGridStruct> &grids, const int &kCell, XYZBoundsStruct &Bounds, FlowReportStruct &FlowReport )
{
	getFlowReportGeneral( grids, kCell, Bounds, FlowReport, XY );
}
void getFlowReportZY( std::vector<DIADGridStruct> &grids, const int &iCell, XYZBoundsStruct &Bounds, FlowReportStruct &FlowReport )
{
	getFlowReportGeneral( grids, iCell, Bounds, FlowReport, ZY );
}
void getFlowReportZX( std::vector<DIADGridStruct> &grids, const int &jCell, XYZBoundsStruct &Bounds, FlowReportStruct &FlowReport )
{
	getFlowReportGeneral( grids, jCell, Bounds, FlowReport, ZX );
}
