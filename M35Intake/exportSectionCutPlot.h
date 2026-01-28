void exportSectionCutPlot( 	MarkerStruct& Marker, DistributionStruct& F, CellCountStruct &cellCount, const size_t &iCut, const int &plotNumber )
{
	std::cout << "Exporting section cut plot " << plotNumber << std::endl;
	auto bouncebackArrayView = Marker.bouncebackArray.getConstView();
	
	auto shifterView = F.shifter.getConstView();
	auto fArrayView  = F.fArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( cellCount.ny, cellCount.nz );
	SectionCut.uxArray.setSizes( cellCount.ny, cellCount.nz );
	SectionCut.uyArray.setSizes( cellCount.ny, cellCount.nz );
	SectionCut.uzArray.setSizes( cellCount.ny, cellCount.nz );
	SectionCut.maskArray.setSizes( cellCount.ny, cellCount.nz );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();
	auto maskArrayView = SectionCut.maskArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 2, int >& doubleIndex) mutable
	{
		const size_t j = doubleIndex.x();
		const size_t k = doubleIndex.y();
		size_t cell = convertIndex(iCut, j, k, cellCount);
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount.n) { shiftedIndex[i] -= cellCount.n; }
		}
		float bouncebackMarker = (float)bouncebackArrayView[cell];
		float f[27];
		float rho, ux, uy, uz;
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		getRhoUxUyUz(rho, ux, uy, uz, f);
		rhoArrayView( j, k ) = rho;
		uxArrayView( j, k ) = ux;
		uyArrayView( j, k ) = uy;
		uzArrayView( j, k ) = uz;
		maskArrayView( j, k ) = bouncebackMarker;
	};
	TNL::Containers::StaticArray< 2, size_t > start{ 0, 0 };
	TNL::Containers::StaticArray< 2, size_t > end{ cellCount.ny, cellCount.nz };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	SectionCutCPU.maskArray = SectionCut.maskArray;
	
	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	// Write metadata first so Python knows the dimensions
	int header[4] = {plotNumber, (int)cellCount.ny, (int)cellCount.nz, (int)5};
	fwrite(header, sizeof(int), 4, fp);
	
	for (size_t j = 0; j < cellCount.ny; j++)
	{
		for (size_t k = 0; k < cellCount.nz; k++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(j, k);
			float ux = SectionCutCPU.uxArray.getElement(j, k);
			float uy = SectionCutCPU.uyArray.getElement(j, k);
			float uz = SectionCutCPU.uzArray.getElement(j, k);
			float mask = SectionCut.maskArray.getElement(j, k);
			float p = 0.f;
			convertToPhysicalUnits(rho, p, ux, uy, uz);
			float data[5] = {p, ux, uy, uz, mask};
			fwrite(data, sizeof(float), 5, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}
