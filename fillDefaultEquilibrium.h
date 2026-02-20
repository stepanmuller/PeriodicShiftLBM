void fillDefaultEquilibrium( FStruct& F, InfoStruct &Info )
{
	std::cout << "Filling fArray with default equilibrium, rho = 1, ux, uy, uz = 0" << std::endl;
	auto fArrayView  = F.fArray.getView();

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
	int end = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

void fillDefaultEquilibrium( FStruct& F, InfoStruct &Info, const float &rho, const float &ux, const float &uy, const float &uz )
{
	std::cout << "Filling fArray with equilibrium per request" << std::endl;
	auto fArrayView  = F.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		float feq[27];
		getFeq(rho, ux, uy, uz, feq);
		for ( int i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	int start = 0;
	int end = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

void fillDefaultEquilibriumFromFunction( FStruct& F, InfoStruct &Info )
{
	std::cout << "Filling fArray with equilibrium per function" << std::endl;
	auto fArrayView  = F.fArray.getView();

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
	int start = 0;
	int end = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
