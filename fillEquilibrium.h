/*
void fillEquilibriumDefault( GridStruct &Grid )
{
	std::cout << "Filling fArray with default equilibrium, rho = 1, ux, uy, uz = 0" << std::endl;
	auto fArrayView  = Grid.fArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		const float rho = 1.f;
		const float ux = 0.f;
		const float uy = 0.f;
		const float uz = 0.f;
		float feq[27];
		getFeq(rho, ux, uy, uz, feq);
		for ( int i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	int start = 0;
	int end = Grid.cellCountX * Grid.cellCountY * Grid.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

void fillEquilibriumPrescribed( GridStruct &Grid, const float &rho, const float &ux, const float &uy, const float &uz )
{
	std::cout << "Filling fArray with equilibrium as prescribed" << std::endl;
	auto fArrayView  = Grid.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		float feq[27];
		getFeq(rho, ux, uy, uz, feq);
		for ( int i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	int start = 0;
	int end = Grid.cellCountX * Grid.cellCountY * Grid.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
*/

void fillEquilibriumFromFunction( GridStruct &Grid )
{
	std::cout << "Filling fArray with equilibrium per function" << std::endl;
	auto fArrayView  = Grid.fArray.getView();
	InfoStruct Info = Grid.Info;
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		int iCell, jCell, kCell;
		getIJKCellIndex( cell, iCell, jCell, kCell, Info );
		float feq[27];
		float rho = 1.f;
		float ux, uy, uz = 0.f;
		getInitialRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Info );
		getFeq(rho, ux, uy, uz, feq);
		for ( int i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Grid.Info.cellCount, cellLambda );
}
