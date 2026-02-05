// Version without marker array
void exportSectionCutPlotXY( FStruct& F, InfoStruct& Info, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountX );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountX );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();

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
		rhoArrayView( jCell, iCell ) = rho;
		uxArrayView( jCell, iCell ) = ux;
		uyArrayView( jCell, iCell ) = uy;
		uzArrayView( jCell, iCell ) = uz;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Info.cellCountX, Info.cellCountY };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Info.cellCountY, (int)Info.cellCountX, 4};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int jCell = 0; jCell < Info.cellCountY; jCell++)
	{
		for (int iCell = 0; iCell < Info.cellCountX; iCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(jCell, iCell);
			float ux = SectionCutCPU.uxArray.getElement(jCell, iCell);
			float uy = SectionCutCPU.uyArray.getElement(jCell, iCell);
			float uz = SectionCutCPU.uzArray.getElement(jCell, iCell);
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
			float data[4] = {p, ux, uy, uz};
			fwrite(data, sizeof(float), 4, fp);
		}
	}
	fclose(fp);
	system("python3 plotterSimple.py");
}

void exportSectionCutPlotZY( FStruct& F, InfoStruct& Info, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountZ );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const IntPairType& doubleIndex ) mutable
	{
		const int jCell = doubleIndex[0];
		const int kCell = doubleIndex[1];
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		rhoArrayView( jCell, kCell ) = rho;
		uxArrayView( jCell, kCell ) = ux;
		uyArrayView( jCell, kCell ) = uy;
		uzArrayView( jCell, kCell ) = uz;
	};
	IntPairType start{ 0, 0 };
	IntPairType end{ Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb"); 	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	int header[4] = {plotNumber, (int)Info.cellCountY, (int)Info.cellCountZ, 4};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int jCell = 0; jCell < Info.cellCountY; jCell++)
	{
		for (int kCell = 0; kCell < Info.cellCountZ; kCell++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(jCell, kCell);
			float ux = SectionCutCPU.uxArray.getElement(jCell, kCell);
			float uy = SectionCutCPU.uyArray.getElement(jCell, kCell);
			float uz = SectionCutCPU.uzArray.getElement(jCell, kCell);
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
			float data[4] = {p, ux, uy, uz};
			fwrite(data, sizeof(float), 4, fp);
		}
	}
	fclose(fp);
	system("python3 plotterSimple.py");
}


// Version with marker array
void exportSectionCutPlotXY( FStruct &F, BoolArrayType &inputMarkerArray, InfoStruct &Info, const int &kCell, const int &plotNumber )
{
	std::cout << "Exporting XY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
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
		const float marker = inputMarkerArrayView( cell );
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
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}

void exportSectionCutPlotZY( FStruct &F, BoolArrayType &inputMarkerArray, InfoStruct &Info, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting ZY section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
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
		const float marker = inputMarkerArrayView( cell );
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
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}

void exportSectionCutPlotZX( FStruct &F, BoolArrayType &inputMarkerArray, InfoStruct &Info, const int &jCell, const int &plotNumber )
{
	std::cout << "Exporting ZX section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	auto inputMarkerArrayView  = inputMarkerArray.getConstView();
	
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
		const float marker = inputMarkerArrayView( cell );
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
			float p;
			convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
			float data[5] = {p, ux, uy, uz, marker};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}
