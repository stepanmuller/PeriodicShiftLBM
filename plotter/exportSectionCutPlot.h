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


/*
void exportSectionCutPlotXY( GridStruct &Grid, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.markerArray.setSizes( Info.cellCountY, Info.cellCountX );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int iCell = doubleIndex[0];
		const int jCell = doubleIndex[1];
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
		rhoArrayView( jCell, iCell ) = rho;
		uxArrayView( jCell, iCell ) = ux;
		uyArrayView( jCell, iCell ) = uy;
		uzArrayView( jCell, iCell ) = uz;
		markerArrayView( jCell, iCell ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Info.cellCountX, Info.cellCountY };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Info.cellCountY, (int)Info.cellCountX, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int jCell = 0; jCell < Info.cellCountY; jCell++)
	{
		for (int iCell = 0; iCell < Info.cellCountX; iCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(jCell, iCell);
			float ux = SectionCutCPU.uxArray.getElement(jCell, iCell);
			float uy = SectionCutCPU.uyArray.getElement(jCell, iCell);
			float uz = SectionCutCPU.uzArray.getElement(jCell, iCell);
			float marker = SectionCutCPU.markerArray.getElement(jCell, iCell);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			const float uHorizontal = ux;
			const float uVertical = uy;
			const float uNormal = uz;
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotZY( GridStruct &Grid, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.markerArray.setSizes( Info.cellCountY, Info.cellCountZ );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int kCell = doubleIndex[0];
		const int jCell = doubleIndex[1];
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
		rhoArrayView( jCell, kCell ) = rho;
		uxArrayView( jCell, kCell ) = ux;
		uyArrayView( jCell, kCell ) = uy;
		uzArrayView( jCell, kCell ) = uz;
		markerArrayView( jCell, kCell ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Info.cellCountZ, Info.cellCountY };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Info.cellCountY, (int)Info.cellCountZ, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int jCell = 0; jCell < Info.cellCountY; jCell++)
	{
		for (int kCell = 0; kCell < Info.cellCountZ; kCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(jCell, kCell);
			float ux = SectionCutCPU.uxArray.getElement(jCell, kCell);
			float uy = SectionCutCPU.uyArray.getElement(jCell, kCell);
			float uz = SectionCutCPU.uzArray.getElement(jCell, kCell);
			float marker = SectionCutCPU.markerArray.getElement(jCell, kCell);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			const float uHorizontal = uz;
			const float uVertical = uy;
			const float uNormal = ux;
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

void exportSectionCutPlotZX( GridStruct &Grid, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting ZX section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountX, Info.cellCountZ );
	SectionCut.uxArray.setSizes( Info.cellCountX, Info.cellCountZ );
	SectionCut.uyArray.setSizes( Info.cellCountX, Info.cellCountZ );
	SectionCut.uzArray.setSizes( Info.cellCountX, Info.cellCountZ );
	SectionCut.markerArray.setSizes( Info.cellCountX, Info.cellCountZ );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto markerArrayView = SectionCut.markerArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int kCell = doubleIndex[0];
		const int iCell = doubleIndex[1];
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
		rhoArrayView( iCell, kCell ) = rho;
		uxArrayView( iCell, kCell ) = ux;
		uyArrayView( iCell, kCell ) = uy;
		uzArrayView( iCell, kCell ) = uz;
		markerArrayView( iCell, kCell ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Info.cellCountZ, Info.cellCountX };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Info.cellCountX, (int)Info.cellCountZ, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int iCell = 0; iCell < Info.cellCountX; iCell++)
	{
		for (int kCell = 0; kCell < Info.cellCountZ; kCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(iCell, kCell);
			float ux = SectionCutCPU.uxArray.getElement(iCell, kCell);
			float uy = SectionCutCPU.uyArray.getElement(iCell, kCell);
			float uz = SectionCutCPU.uzArray.getElement(iCell, kCell);
			float marker = SectionCutCPU.markerArray.getElement(iCell, kCell);
			float p = rho;
			convertToPhysicalVelocity( ux, uy, uz, Info );
			convertToPhysicalPressure( p, Info );
			const float uHorizontal = uz;
			const float uVertical = ux;
			const float uNormal = uy;
			float data[5] = {p, uHorizontal, uVertical, uNormal, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
}

*/
