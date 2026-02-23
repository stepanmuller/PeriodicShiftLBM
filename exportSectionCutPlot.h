// Version without marker array
void exportSectionCutPlotXY( GridStruct &Grid, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
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
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}


/*
void exportSectionCutPlotZY( BoolArrayType &inputMarkerArray, GridStruct &Grid, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Grid.Info.cellCountY, Grid.Info.cellCountZ );
	SectionCut.uxArray.setSizes( Grid.Info.cellCountY, Grid.Info.cellCountZ );
	SectionCut.uyArray.setSizes( Grid.Info.cellCountY, Grid.Info.cellCountZ );
	SectionCut.uzArray.setSizes( Grid.Info.cellCountY, Grid.Info.cellCountZ );
	SectionCut.markerArray.setSizes( Grid.Info.cellCountY, Grid.Info.cellCountZ );
		
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
		getCellIndex( cell, iCell, jCell, kCell, Grid );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Grid );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		const float marker = inputMarkerArrayView( cell );
		rhoArrayView( jCell, kCell ) = rho;
		uxArrayView( jCell, kCell ) = ux;
		uyArrayView( jCell, kCell ) = uy;
		uzArrayView( jCell, kCell ) = uz;
		markerArrayView( jCell, kCell ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Grid.Info.cellCountZ, Grid.Info.cellCountY };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Grid.Info.cellCountY, (int)Grid.Info.cellCountZ, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int jCell = 0; jCell < Grid.Info.cellCountY; jCell++)
	{
		for (int kCell = 0; kCell < Grid.Info.cellCountZ; kCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(jCell, kCell);
			float ux = SectionCutCPU.uxArray.getElement(jCell, kCell);
			float uy = SectionCutCPU.uyArray.getElement(jCell, kCell);
			float uz = SectionCutCPU.uzArray.getElement(jCell, kCell);
			float marker = SectionCutCPU.markerArray.getElement(jCell, kCell);
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Grid );
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}

void exportSectionCutPlotZX( BoolArrayType &inputMarkerArray, GridStruct &Grid, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting ZX section cut plot " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Grid.Info.cellCountX, Grid.Info.cellCountZ );
	SectionCut.uxArray.setSizes( Grid.Info.cellCountX, Grid.Info.cellCountZ );
	SectionCut.uyArray.setSizes( Grid.Info.cellCountX, Grid.Info.cellCountZ );
	SectionCut.uzArray.setSizes( Grid.Info.cellCountX, Grid.Info.cellCountZ );
	SectionCut.markerArray.setSizes( Grid.Info.cellCountX, Grid.Info.cellCountZ );
		
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
		getCellIndex( cell, iCell, jCell, kCell, Grid );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Grid );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		const float marker = inputMarkerArrayView( cell );
		rhoArrayView( iCell, kCell ) = rho;
		uxArrayView( iCell, kCell ) = ux;
		uyArrayView( iCell, kCell ) = uy;
		uzArrayView( iCell, kCell ) = uz;
		markerArrayView( iCell, kCell ) = marker;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Grid.Info.cellCountZ, Grid.Info.cellCountX };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.markerArray = SectionCut.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Grid.Info.cellCountX, (int)Grid.Info.cellCountZ, 5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int iCell = 0; iCell < Grid.Info.cellCountX; iCell++)
	{
		for (int kCell = 0; kCell < Grid.Info.cellCountZ; kCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(iCell, kCell);
			float ux = SectionCutCPU.uxArray.getElement(iCell, kCell);
			float uy = SectionCutCPU.uyArray.getElement(iCell, kCell);
			float uz = SectionCutCPU.uzArray.getElement(iCell, kCell);
			float marker = SectionCutCPU.markerArray.getElement(iCell, kCell);
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Grid );
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}

void exportSection3DPlot( BoolArrayType &inputMarkerArray, GridStruct &Grid, const int &iMin, const int &jMin, const int &kMin, const int &iMax, const int &jMax, const int &kMax, const int &plotNumber )
{
	std::cout << "Exporting Section 3D " << plotNumber << std::endl;
	auto fArrayView  = Grid.fArray.getConstView();
	auto shifterView  = Grid.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
	const int sectionSizeX = iMax - iMin + 1;
	const int sectionSizeY = jMax - jMin + 1;
	const int sectionSizeZ = kMax - kMin + 1;
	
	Section3DStruct Section3D;
	Section3D.rhoArray.setSizes( sectionSizeX, sectionSizeY, sectionSizeZ );
	Section3D.uxArray.setSizes( sectionSizeX, sectionSizeY, sectionSizeZ );
	Section3D.uyArray.setSizes( sectionSizeX, sectionSizeY, sectionSizeZ );
	Section3D.uzArray.setSizes( sectionSizeX, sectionSizeY, sectionSizeZ );
	Section3D.markerArray.setSizes( sectionSizeX, sectionSizeY, sectionSizeZ );
		
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
		getCellIndex( cell, iCell, jCell, kCell, Grid );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Grid );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		const float marker = inputMarkerArrayView( cell );
		rhoArrayView( iCell-iMin, jCell-jMin, kCell-kMin ) = rho;
		uxArrayView( iCell-iMin, jCell-jMin, kCell-kMin ) = ux;
		uyArrayView( iCell-iMin, jCell-jMin, kCell-kMin ) = uy;
		uzArrayView( iCell-iMin, jCell-jMin, kCell-kMin ) = uz;
		markerArrayView( iCell-iMin, jCell-jMin, kCell-kMin ) = marker;
	};
	IntTripleType start{ iMin, jMin, kMin };
	IntTripleType end{ iMax+1, jMax+1, kMax+1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	Section3DStructCPU Section3DCPU;
	Section3DCPU.rhoArray = Section3D.rhoArray;
	Section3DCPU.uxArray = Section3D.uxArray;
	Section3DCPU.uyArray = Section3D.uyArray;
	Section3DCPU.uzArray = Section3D.uzArray;
	Section3DCPU.markerArray = Section3D.markerArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[5] = {plotNumber, sectionSizeX, sectionSizeY, sectionSizeZ, 5};
	fwrite(header, sizeof(int), 5, fp);
	
	for ( int iCell = 0; iCell < sectionSizeX; iCell++ )
	{
		for ( int jCell = 0; jCell < sectionSizeY; jCell++ )
		{
			for ( int kCell = 0; kCell < sectionSizeZ; kCell++ )
			{
				float rho = Section3DCPU.rhoArray.getElement( iCell, jCell, kCell );
				float ux = Section3DCPU.uxArray.getElement( iCell, jCell, kCell );
				float uy = Section3DCPU.uyArray.getElement( iCell, jCell, kCell );
				float uz = Section3DCPU.uzArray.getElement( iCell, jCell, kCell );
				float marker = Section3DCPU.markerArray.getElement( iCell, jCell, kCell );
				float p;
				convertToPhysicalUnits( rho, p, ux, uy, uz, Grid );
				float data[5] = {p, ux, uy, uz, marker};
				fwrite(data, sizeof(float), 5, fp);
			}
		}
	}
	fclose(fp);
	system("python3 plotter3D.py");
}

*/
