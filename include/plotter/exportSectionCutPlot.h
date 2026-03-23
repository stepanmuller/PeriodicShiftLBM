enum PlaneEnum { XY, ZY, ZX };
void exportSectionCutPlotGeneral( GridStruct &Grid, const int &cutIndex, const int &plotNumber, PlaneEnum plane )
{
	InfoStruct Info = Grid.Info;
	int cellCountHorizontal, cellCountVertical;
	if ( plane == XY ) 		{ cellCountHorizontal = Info.cellCountX; cellCountVertical = Info.cellCountY; }
	else if ( plane == ZY ) { cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountY; }
	else 					{ cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountX; }
	
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uxArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uyArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uzArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.markerArray.setSizes( cellCountVertical, cellCountHorizontal );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int indexHorizontal = doubleIndex[0];
		const int indexVertical = doubleIndex[1];
		int iCell, jCell, kCell, cell;
		if ( plane == XY ) 		{ iCell = indexHorizontal; jCell = indexVertical; kCell = cutIndex; }
		else if ( plane == ZY ) { iCell = cutIndex; jCell = indexVertical; kCell = indexHorizontal; }
		else 					{ iCell = indexVertical; jCell = cutIndex; kCell = indexHorizontal; }
		getCellIndex( cell, iCell, jCell, kCell, Info );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		const float marker = Marker.bounceback;
		rhoArrayView( indexVertical, indexHorizontal ) = rho;
		uxArrayView( indexVertical, indexHorizontal ) = ux;
		uyArrayView( indexVertical, indexHorizontal ) = uy;
		uzArrayView( indexVertical, indexHorizontal ) = uz;
		markerArrayView( indexVertical, indexHorizontal ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ cellCountHorizontal, cellCountVertical };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)cellCountVertical, (int)cellCountHorizontal, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int indexVertical = 0; indexVertical < cellCountVertical; indexVertical++)
	{
		for (int indexHorizontal = 0; indexHorizontal < cellCountHorizontal; indexHorizontal++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(indexVertical, indexHorizontal);
			float ux = SectionCutCPU.uxArray.getElement(indexVertical, indexHorizontal);
			float uy = SectionCutCPU.uyArray.getElement(indexVertical, indexHorizontal);
			float uz = SectionCutCPU.uzArray.getElement(indexVertical, indexHorizontal);
			float marker = SectionCutCPU.markerArray.getElement(indexVertical, indexHorizontal);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			float uHorizontal, uVertical, uNormal;
			if ( plane == XY ) 		{ uHorizontal = ux; uVertical = uy; uNormal = uz; }
			else if ( plane == ZY ) { uHorizontal = uz; uVertical = uy; uNormal = ux; }
			else 					{ uHorizontal = uz; uVertical = ux; uNormal = uy; }
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotXY( GridStruct &Grid, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, kCell, plotNumber, XY );
}
void exportSectionCutPlotZY( GridStruct &Grid, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, iCell, plotNumber, ZY );
}
void exportSectionCutPlotZX( GridStruct &Grid, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting ZX section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, jCell, plotNumber, ZX );
}

// Version with scalar transport //
void exportSectionCutPlotGeneral( GridStruct &Grid, ScalarTransportStruct &ScalarTransport, const int &cutIndex, const int &plotNumber, PlaneEnum plane )
{
	InfoStruct Info = Grid.Info;
	int cellCountHorizontal, cellCountVertical;
	if ( plane == XY ) 		{ cellCountHorizontal = Info.cellCountX; cellCountVertical = Info.cellCountY; }
	else if ( plane == ZY ) { cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountY; }
	else 					{ cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountX; }
	
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	auto TArrayView = ScalarTransport.TArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uxArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uyArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uzArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.markerArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.scalarTransportArray.setSizes( cellCountVertical, cellCountHorizontal );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();
	auto scalarTransportArrayView = SectionCut.scalarTransportArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int indexHorizontal = doubleIndex[0];
		const int indexVertical = doubleIndex[1];
		int iCell, jCell, kCell, cell;
		if ( plane == XY ) 		{ iCell = indexHorizontal; jCell = indexVertical; kCell = cutIndex; }
		else if ( plane == ZY ) { iCell = cutIndex; jCell = indexVertical; kCell = indexHorizontal; }
		else 					{ iCell = indexVertical; jCell = cutIndex; kCell = indexHorizontal; }
		getCellIndex( cell, iCell, jCell, kCell, Info );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float T[27];
		for (int direction = 0; direction < 27; direction++) T[direction] = TArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz, scalarTransport;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		getScalarTransport(scalarTransport, f, T);
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		const float marker = Marker.bounceback;
		rhoArrayView( indexVertical, indexHorizontal ) = rho;
		uxArrayView( indexVertical, indexHorizontal ) = ux;
		uyArrayView( indexVertical, indexHorizontal ) = uy;
		uzArrayView( indexVertical, indexHorizontal ) = uz;
		markerArrayView( indexVertical, indexHorizontal ) = marker;
		scalarTransportArrayView( indexVertical, indexHorizontal ) = scalarTransport;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ cellCountHorizontal, cellCountVertical };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	SectionCutCPU.scalarTransportArray = SectionCut.scalarTransportArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)cellCountVertical, (int)cellCountHorizontal, 6};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int indexVertical = 0; indexVertical < cellCountVertical; indexVertical++)
	{
		for (int indexHorizontal = 0; indexHorizontal < cellCountHorizontal; indexHorizontal++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(indexVertical, indexHorizontal);
			float ux = SectionCutCPU.uxArray.getElement(indexVertical, indexHorizontal);
			float uy = SectionCutCPU.uyArray.getElement(indexVertical, indexHorizontal);
			float uz = SectionCutCPU.uzArray.getElement(indexVertical, indexHorizontal);
			float marker = SectionCutCPU.markerArray.getElement(indexVertical, indexHorizontal);
			float scalarTransport = SectionCutCPU.scalarTransportArray.getElement(indexVertical, indexHorizontal);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			float uHorizontal, uVertical, uNormal;
			if ( plane == XY ) 		{ uHorizontal = ux; uVertical = uy; uNormal = uz; }
			else if ( plane == ZY ) { uHorizontal = uz; uVertical = uy; uNormal = ux; }
			else 					{ uHorizontal = uz; uVertical = ux; uNormal = uy; }
			float data[6] = {p, uHorizontal, uVertical, uNormal, marker, scalarTransport};
			fwrite(data, sizeof(float), 6, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotXY( GridStruct &Grid, ScalarTransportStruct &ScalarTransport, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, ScalarTransport, kCell, plotNumber, XY );
}
void exportSectionCutPlotZY( GridStruct &Grid, ScalarTransportStruct &ScalarTransport, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, ScalarTransport, iCell, plotNumber, ZY );
}
void exportSectionCutPlotZX( GridStruct &Grid, ScalarTransportStruct &ScalarTransport, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting ZX section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, ScalarTransport, jCell, plotNumber, ZX );
}

// DIAD version
void exportSectionCutPlotGeneral( DIADGridStruct &Grid, const int &cutIndex, const int &plotNumber, PlaneEnum plane )
{
	InfoStruct Info = Grid.Info;
	int cellCountHorizontal, cellCountVertical;
	if ( plane == XY ) 		{ cellCountHorizontal = Info.cellCountX; cellCountVertical = Info.cellCountY; }
	else if ( plane == ZY ) { cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountY; }
	else 					{ cellCountHorizontal = Info.cellCountZ; cellCountVertical = Info.cellCountX; }
	
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
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uxArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uyArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uzArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.markerArray.setSizes( cellCountVertical, cellCountHorizontal );
	
	SectionCut.rhoArray.setValue( 1.f );
	SectionCut.uxArray.setValue( 0.f );
	SectionCut.uyArray.setValue( 0.f );
	SectionCut.uzArray.setValue( 0.f );
	SectionCut.markerArray.setValue( 1 );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell = iView[ cell ];
		int jCell = jView[ cell ];
		int kCell = kView[ cell ];
		
		DIADEsotwistNbrStruct Nbr;
		Nbr.i = iNbrView( cell );
		Nbr.j = jNbrView( cell );
		Nbr.k = kNbrView( cell );
		Nbr.ij = ijNbrView( cell );
		Nbr.ik = ikNbrView( cell );
		Nbr.jk = jkNbrView( cell );
		Nbr.ijk = ijkNbrView( cell );
		
		int indexHorizontal = 0;
		int indexVertical = 0;
		
		if ( plane == XY ) 
		{
			if (kCell != cutIndex) return;
			indexHorizontal = iCell; 
			indexVertical = jCell; 
		}
		else if ( plane == ZY ) 
		{ 
			if (iCell != cutIndex) return; 
			indexVertical = jCell; 
			indexHorizontal = kCell; 
		}
		else 					
		{ 
			if (jCell != cutIndex) return; 
			indexVertical = iCell; 
			indexHorizontal = kCell; 
		}
		
		float f[27];
		int cellReadIndex[27];
		int fReadIndex[27];
		getEsotwistWriteIndex( cell, cellReadIndex, fReadIndex, Nbr, esotwistFlipper, Info ); 
		// Using the write index because plotting is happening after collision and we want to be consistent with that last write
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fReadIndex[direction], cellReadIndex[direction]);
		
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);

		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		//getMarkers( iCell, jCell, kCell, Marker, Info );
		const float marker = Marker.bounceback;
		rhoArrayView( indexVertical, indexHorizontal ) = rho;
		uxArrayView( indexVertical, indexHorizontal ) = ux;
		uyArrayView( indexVertical, indexHorizontal ) = uy;
		uzArrayView( indexVertical, indexHorizontal ) = uz;
		markerArrayView( indexVertical, indexHorizontal ) = marker;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)cellCountVertical, (int)cellCountHorizontal, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int indexVertical = 0; indexVertical < cellCountVertical; indexVertical++)
	{
		for (int indexHorizontal = 0; indexHorizontal < cellCountHorizontal; indexHorizontal++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(indexVertical, indexHorizontal);
			float ux = SectionCutCPU.uxArray.getElement(indexVertical, indexHorizontal);
			float uy = SectionCutCPU.uyArray.getElement(indexVertical, indexHorizontal);
			float uz = SectionCutCPU.uzArray.getElement(indexVertical, indexHorizontal);
			float marker = SectionCutCPU.markerArray.getElement(indexVertical, indexHorizontal);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			float uHorizontal, uVertical, uNormal;
			if ( plane == XY ) 		{ uHorizontal = ux; uVertical = uy; uNormal = uz; }
			else if ( plane == ZY ) { uHorizontal = uz; uVertical = uy; uNormal = ux; }
			else 					{ uHorizontal = uz; uVertical = ux; uNormal = uy; }
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotXY( DIADGridStruct &Grid, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD XY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, kCell, plotNumber, XY );
}
void exportSectionCutPlotZY( DIADGridStruct &Grid, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD ZY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, iCell, plotNumber, ZY );
}
void exportSectionCutPlotZX( DIADGridStruct &Grid, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD ZX section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( Grid, jCell, plotNumber, ZX );
}

// DIAD Version that plots all grids into one hi res pic
void exportSectionCutPlotGeneral( std::vector<DIADGridStruct> &grids, const int &cutIndex, const int &plotNumber, PlaneEnum plane )
{
	const int levelCount = grids.size();
	
	int cellCountHorizontal, cellCountVertical;
	if ( plane == XY ) 		{ cellCountHorizontal = grids[levelCount-1].Info.cellCountX; cellCountVertical = grids[levelCount-1].Info.cellCountY; }
	else if ( plane == ZY ) { cellCountHorizontal = grids[levelCount-1].Info.cellCountZ; cellCountVertical = grids[levelCount-1].Info.cellCountY; }
	else 					{ cellCountHorizontal = grids[levelCount-1].Info.cellCountZ; cellCountVertical = grids[levelCount-1].Info.cellCountX; }
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uxArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uyArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uzArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.markerArray.setSizes( cellCountVertical, cellCountHorizontal );
	
	SectionCut.rhoArray.setValue( 1.f );
	SectionCut.uxArray.setValue( 0.f );
	SectionCut.uyArray.setValue( 0.f );
	SectionCut.uzArray.setValue( 0.f );
	SectionCut.markerArray.setValue( 1 );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();
	
	for ( int level = 0; level < levelCount; level++ )
	{
		DIADGridStruct &Grid = grids[level];
		InfoStruct Info = Grid.Info;
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
			int iCell = iView[ cell ] * cellScale; // recalculating to index on the finest grid
			int jCell = jView[ cell ] * cellScale;
			int kCell = kView[ cell ] * cellScale;
			
			DIADEsotwistNbrStruct Nbr;
			Nbr.i = iNbrView( cell );
			Nbr.j = jNbrView( cell );
			Nbr.k = kNbrView( cell );
			Nbr.ij = ijNbrView( cell );
			Nbr.ik = ikNbrView( cell );
			Nbr.jk = jkNbrView( cell );
			Nbr.ijk = ijkNbrView( cell );
			
			int indexHorizontal = 0;
			int indexVertical = 0;
			
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
			
			float f[27];
			int cellReadIndex[27];
			int fReadIndex[27];
			getEsotwistWriteIndex( cell, cellReadIndex, fReadIndex, Nbr, esotwistFlipper, Info ); 
			// Using the write index because plotting is happening after collision and we want to be consistent with that last write
			for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fReadIndex[direction], cellReadIndex[direction]);
			
			float rho, ux, uy, uz;
			getRhoUxUyUz(rho, ux, uy, uz, f);

			MarkerStruct Marker;
			if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
			//getMarkers( iCell, jCell, kCell, Marker, Info );
			const float marker = Marker.bounceback;
			
			for ( int shiftVertical = 0; shiftVertical < cellScale; shiftVertical++ )
			{
				for ( int shiftHorizontal = 0; shiftHorizontal < cellScale; shiftHorizontal++ )
				{
					rhoArrayView( indexVertical + shiftVertical, indexHorizontal + shiftHorizontal ) = rho;
					uxArrayView( indexVertical + shiftVertical, indexHorizontal + shiftHorizontal ) = ux;
					uyArrayView( indexVertical + shiftVertical, indexHorizontal + shiftHorizontal ) = uy;
					uzArrayView( indexVertical + shiftVertical, indexHorizontal + shiftHorizontal ) = uz;
					markerArrayView( indexVertical + shiftVertical, indexHorizontal + shiftHorizontal ) = marker;
				}
			}
		};
		TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
	}
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)cellCountVertical, (int)cellCountHorizontal, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	InfoStruct Info = grids[0].Info; // it is needed for the unit conversion below
	
	for (int indexVertical = 0; indexVertical < cellCountVertical; indexVertical++)
	{
		for (int indexHorizontal = 0; indexHorizontal < cellCountHorizontal; indexHorizontal++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(indexVertical, indexHorizontal);
			float ux = SectionCutCPU.uxArray.getElement(indexVertical, indexHorizontal);
			float uy = SectionCutCPU.uyArray.getElement(indexVertical, indexHorizontal);
			float uz = SectionCutCPU.uzArray.getElement(indexVertical, indexHorizontal);
			float marker = SectionCutCPU.markerArray.getElement(indexVertical, indexHorizontal);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			float uHorizontal, uVertical, uNormal;
			if ( plane == XY ) 		{ uHorizontal = ux; uVertical = uy; uNormal = uz; }
			else if ( plane == ZY ) { uHorizontal = uz; uVertical = uy; uNormal = ux; }
			else 					{ uHorizontal = uz; uVertical = ux; uNormal = uy; }
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotXY( std::vector<DIADGridStruct> &grids, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD XY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( grids, kCell, plotNumber, XY );
}
void exportSectionCutPlotZY( std::vector<DIADGridStruct> &grids, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD ZY section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( grids, iCell, plotNumber, ZY );
}
void exportSectionCutPlotZX( std::vector<DIADGridStruct> &grids, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting DIAD ZX section cut plot " << plotNumber << std::endl;
	exportSectionCutPlotGeneral( grids, jCell, plotNumber, ZX );
}

void export3DPlot( GridStruct &Grid, const int &plotNumber )
{
	std::cout << "Exporting Section 3D " << plotNumber << std::endl;
	
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	
	Section3DStruct Section3D;
	Section3D.rhoArray.setSizes( Info.cellCountX, Info.cellCountY, Info.cellCountZ );
	Section3D.uxArray.setSizes( Info.cellCountX, Info.cellCountY, Info.cellCountZ );
	Section3D.uyArray.setSizes( Info.cellCountX, Info.cellCountY, Info.cellCountZ );
	Section3D.uzArray.setSizes( Info.cellCountX, Info.cellCountY, Info.cellCountZ );
	Section3D.markerArray.setSizes( Info.cellCountX, Info.cellCountY, Info.cellCountZ );
		
	auto rhoArrayView = Section3D.rhoArray.getView();
	auto uxArrayView = Section3D.uxArray.getView();
	auto uyArrayView = Section3D.uyArray.getView();
	auto uzArrayView = Section3D.uzArray.getView();
	auto markerArrayView = Section3D.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntTripleType& tripleIndex ) mutable
	{
		const int iCell = tripleIndex[0];
		const int jCell = tripleIndex[1];
		const int kCell = tripleIndex[2];
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		const float marker = Marker.bounceback;
		rhoArrayView( iCell, jCell, kCell ) = rho;
		uxArrayView( iCell, jCell, kCell ) = ux;
		uyArrayView( iCell, jCell, kCell ) = uy;
		uzArrayView( iCell, jCell, kCell ) = uz;
		markerArrayView( iCell, jCell, kCell ) = marker;
	};
	IntTripleType start{ 0, 0, 0 };
	IntTripleType end{ Info.cellCountX, Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	Section3DStructCPU Section3DCPU;
	Section3DCPU.rhoArray = Section3D.rhoArray;
	Section3DCPU.uxArray = Section3D.uxArray;
	Section3DCPU.uyArray = Section3D.uyArray;
	Section3DCPU.uzArray = Section3D.uzArray;
	Section3DCPU.markerArray = Section3D.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[5] = {plotNumber, Info.cellCountX, Info.cellCountY, Info.cellCountZ, 5};
	fwrite(header, sizeof(int), 5, fp);
	
	for ( int iCell = 0; iCell < Info.cellCountX; iCell++ )
	{
		for ( int jCell = 0; jCell < Info.cellCountY; jCell++ )
		{
			for ( int kCell = 0; kCell < Info.cellCountZ; kCell++ )
			{
				float rho = Section3DCPU.rhoArray.getElement( iCell, jCell, kCell );
				float ux = Section3DCPU.uxArray.getElement( iCell, jCell, kCell );
				float uy = Section3DCPU.uyArray.getElement( iCell, jCell, kCell );
				float uz = Section3DCPU.uzArray.getElement( iCell, jCell, kCell );
				float marker = Section3DCPU.markerArray.getElement( iCell, jCell, kCell );
				float p = rho;
				convertToPhysicalVelocity( ux, uy, uz, Info );
				convertToPhysicalPressure( p, Info );
				float data[5] = {p, ux, uy, uz, marker};
				fwrite(data, sizeof(float), 5, fp);
			}
		}
	}
	fclose(fp);
}

void exportSectionCutPlotToiletPaperZ( GridStruct &Grid, const float &r, const int &plotNumber )
{
	std::cout << "Exporting Toilet Paper Z section cut plot " << plotNumber << " at radius " << r << " mm" << std::endl;

	InfoStruct Info = Grid.Info;
	const float PI = 3.14159265359f;

	// Dimensions: Horizontal is Z-axis, Vertical is the unrolled circumference (S)
	int cellCountHorizontal = Info.cellCountZ;
	int cellCountVertical = (int)( (2.0f * PI * r) / Info.res );

	auto fArrayView = Grid.fArray.getConstView();
	auto shifterView = Grid.shifter.getConstView();
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();

	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uxArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uyArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.uzArray.setSizes( cellCountVertical, cellCountHorizontal );
	SectionCut.markerArray.setSizes( cellCountVertical, cellCountHorizontal );

	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int indexVertical = doubleIndex[1]; // Vertical (Circumference)
		const int indexHorizontal = doubleIndex[0]; // Horizontal (Z-axis)

		// 1. Calculate the angle theta based on the unrolled distance S
		float s = (float)indexVertical * Info.res;
		float theta = s / r;

		// 2. Map polar to Cartesian
		const float x = r * cosf(theta);
		const float y = r * sinf(theta);
		const float zTemp = 0.f;
		
		int iCell, jCell, kCell;
		getIJKCellIndexFromXYZ( iCell, jCell, kCell, x, y, zTemp, Info);
		kCell = indexHorizontal;

		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );

		rhoArrayView( indexVertical, indexHorizontal ) = rho;
		uxArrayView( indexVertical, indexHorizontal ) = ux;
		uyArrayView( indexVertical, indexHorizontal ) = uy;
		uzArrayView( indexVertical, indexHorizontal ) = uz;
		markerArrayView( indexVertical, indexHorizontal ) = Marker.bounceback;
	};

	IntPairType start{ 0, 0 };
	IntPairType end{ cellCountHorizontal, cellCountVertical };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );

	// Copy to CPU and save to binary file
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;

	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	int header[4] = {plotNumber, cellCountVertical, cellCountHorizontal, 5};
	fwrite(header, sizeof(int), 4, fp);

	for (int indexVertical = 0; indexVertical < cellCountVertical; indexVertical++)
	{
		for (int indexHorizontal = 0; indexHorizontal < cellCountHorizontal; indexHorizontal++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(indexVertical, indexHorizontal);
			float ux = SectionCutCPU.uxArray.getElement(indexVertical, indexHorizontal);
			float uy = SectionCutCPU.uyArray.getElement(indexVertical, indexHorizontal);
			float uz = SectionCutCPU.uzArray.getElement(indexVertical, indexHorizontal);
			float marker = SectionCutCPU.markerArray.getElement(indexVertical, indexHorizontal);
			
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );

			// Velocity Mapping: 
			// In the "unrolled" plot, Horizontal velocity is Uz.
			// Vertical velocity is the tangential velocity (u_theta).
			float theta = ((float)indexVertical * Info.res) / r;
			float uTangential = -ux * sinf(theta) + uy * cosf(theta);
			float uRadial = ux * cosf(theta) + uy * sinf(theta);

			float data[5] = {p, uz, uTangential, uRadial, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}
