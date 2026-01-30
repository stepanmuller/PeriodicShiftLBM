void exportSectionCutPlot( FloatArray4DType& fArray, InfoStruct &Info, const int &iCell, const int &plotNumber )
{
	std::cout << "Exporting section cut plot " << plotNumber << std::endl;
	auto fArrayView  = fArray.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountZ );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const DoubleIndexType &doubleIndex) mutable
	{
		const int jCell = doubleIndex.x();
		const int kCell = doubleIndex.y();
		int iStreamed[27], jStreamed[27], kStreamed[27];
		getStreamedIndexes( iCell, jCell, kCell, iStreamed, jStreamed, kStreamed, Info );
		float f[27];
		for ( int direction = 0; direction < 27; direction++ ) f[direction] = fArrayView( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] );
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		rhoArrayView( jCell, kCell ) = rho;
		uxArrayView( jCell, kCell ) = ux;
		uyArrayView( jCell, kCell ) = uy;
		uzArrayView( jCell, kCell ) = uz;
	};
	DoubleIndexType start{ 0, 0 };
	DoubleIndexType end{ Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	
	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	// Write metadata first so Python knows the dimensions
	int header[4] = {plotNumber, Info.cellCountY, Info.cellCountZ, 4};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int j = 0; j < Info.cellCountY; j++)
	{
		for (int k = 0; k < Info.cellCountZ; k++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(j, k);
			float ux = SectionCutCPU.uxArray.getElement(j, k);
			float uy = SectionCutCPU.uyArray.getElement(j, k);
			float uz = SectionCutCPU.uzArray.getElement(j, k);
			float data[4] = {rho, ux, uy, uz};
			fwrite(data, sizeof(float), 4, fp);
		}
	}
	fclose(fp);
	system("python3 plotter.py");
}
