// DIAD Version that dynamically downsamples finest grids and clips to specified bounds
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
	const int targetScale = 1 << (levelCount - targetLevelCount); 

	// 2. Map the physical Bounds to index bounds on the absolute finest grid
	InfoStruct finestInfo = grids[levelCount - 1].Info;
	
	int iMin = std::max(0, static_cast<int>(std::floor((Bounds.xmin - finestInfo.ox) / finestInfo.res)));
	int iMax = std::min(finestInfo.cellCountX, static_cast<int>(std::ceil((Bounds.xmax - finestInfo.ox) / finestInfo.res)));
	int jMin = std::max(0, static_cast<int>(std::floor((Bounds.ymin - finestInfo.oy) / finestInfo.res)));
	int jMax = std::min(finestInfo.cellCountY, static_cast<int>(std::ceil((Bounds.ymax - finestInfo.oy) / finestInfo.res)));
	int kMin = std::max(0, static_cast<int>(std::floor((Bounds.zmin - finestInfo.oz) / finestInfo.res)));
	int kMax = std::min(finestInfo.cellCountZ, static_cast<int>(std::ceil((Bounds.zmax - finestInfo.oz) / finestInfo.res)));

	// Assign horizontal and vertical bounds based on the chosen plane
	int hMinFinest = 0, hMaxFinest = 0;
	int vMinFinest = 0, vMaxFinest = 0;

	if ( plane == XY ) 
	{
		hMinFinest = iMin; hMaxFinest = iMax;
		vMinFinest = jMin; vMaxFinest = jMax;
	} 
	else if ( plane == ZY ) 
	{
		hMinFinest = kMin; hMaxFinest = kMax;
		vMinFinest = jMin; vMaxFinest = jMax;
	} 
	else // ZX plane
	{
		hMinFinest = kMin; hMaxFinest = kMax;
		vMinFinest = iMin; vMaxFinest = iMax;
	}

	// Scale the finest index bounds down to the target output resolution
	int hMinTarget = hMinFinest / targetScale;
	int hMaxTarget = (hMaxFinest + targetScale - 1) / targetScale; // Ceil division
	int vMinTarget = vMinFinest / targetScale;
	int vMaxTarget = (vMaxFinest + targetScale - 1) / targetScale; // Ceil division

	int targetWidth = std::max(0, hMaxTarget - hMinTarget);
	int targetHeight = std::max(0, vMaxTarget - vMinTarget);

	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( targetHeight, targetWidth );
	SectionCut.uxArray.setSizes( targetHeight, targetWidth );
	SectionCut.uyArray.setSizes( targetHeight, targetWidth );
	SectionCut.uzArray.setSizes( targetHeight, targetWidth );
	SectionCut.markerArray.setSizes( targetHeight, targetWidth );
	SectionCut.gridIDArray.setSizes( targetHeight, targetWidth );
	
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
	
	// 3. Loop through ALL grids (none are dropped)
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
			
			// 4. Mapping coordinates to the scaled-down and specifically cropped output array
			int outYStart = (indexVertical / targetScale) - vMinTarget;
			int outXStart = (indexHorizontal / targetScale) - hMinTarget;
			
			// How many pixels this cell spans on the output array (minimum 1 for downsampled fine grids)
			int spanY = max(1, cellScale / targetScale);
			int spanX = max(1, cellScale / targetScale);
			
			for ( int shiftY = 0; shiftY < spanY; shiftY++ )
			{
				int y = outYStart + shiftY;
				if ( y < 0 || y >= targetHeight ) continue; // Crop
				
				for ( int shiftX = 0; shiftX < spanX; shiftX++ )
				{
					int x = outXStart + shiftX;
					if ( x < 0 || x >= targetWidth ) continue; // Crop
					
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
	
	// AT THIS POINT THE SECTION CUT IS FILLED, DO REDUCTION ON ITS ARRAYS AND FILL FLOW REPORT
	const int totalCutCells = targetWidth * targetHeight;

	auto fetchCellCount = [=] __cuda_callable__( const int singleIndex )
	{
		const int y = singleIndex / targetWidth;
		const int x = singleIndex % targetWidth;
		if ( markerArrayView( y, x ) ) return 0;
		else return 1;
	};
	auto reductionCellCount = [] __cuda_callable__( const int& a, const int& b )
	{
		return a + b;
	};
	auto fetchUx = [=] __cuda_callable__( const int singleIndex )
	{
		const int y = singleIndex / targetWidth;
		const int x = singleIndex % targetWidth;
		if ( markerArrayView( y, x ) ) return 0.f;
		else return uxArrayView( y, x );
	};
	auto fetchUy = [=] __cuda_callable__( const int singleIndex )
	{
		const int y = singleIndex / targetWidth;
		const int x = singleIndex % targetWidth;
		if ( markerArrayView( y, x ) ) return 0.f;
		else return uyArrayView( y, x );
	};
	auto fetchUz = [=] __cuda_callable__( const int singleIndex )
	{
		const int y = singleIndex / targetWidth;
		const int x = singleIndex % targetWidth;
		if ( markerArrayView( y, x ) ) return 0.f;
		else return uzArrayView( y, x );
	};
	auto fetchRho = [=] __cuda_callable__( const int singleIndex )
	{
		const int y = singleIndex / targetWidth;
		const int x = singleIndex % targetWidth;
		if ( markerArrayView( y, x ) ) return 0.f;
		else return rhoArrayView( y, x );
	};
	auto reductionFloat = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};

	const int cellSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, totalCutCells, fetchCellCount, reductionCellCount, 0 );
	float uxSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, totalCutCells, fetchUx, reductionFloat, 0.f );
	float uySum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, totalCutCells, fetchUy, reductionFloat, 0.f );
	float uzSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, totalCutCells, fetchUz, reductionFloat, 0.f );
	float rhoSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, totalCutCells, fetchRho, reductionFloat, 0.f );
	
	FlowReport.areamm2 = cellSum * ( grids[targetLevelCount-1].Info.res * grids[targetLevelCount-1].Info.res );
	FlowReport.ux = uxSum / (float)cellSum;
	FlowReport.uy = uySum / (float)cellSum;
	FlowReport.uz = uzSum / (float)cellSum;
	FlowReport.rho = rhoSum / (float)cellSum;
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
