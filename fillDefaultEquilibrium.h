void fillDefaultEquilibrium( FStruct& F, InfoStruct &Info )
{
	std::cout << "Filling fArray with default equilibrium, rho = 1, ux, uy, uz = 0" << std::endl;
	auto fArrayView  = F.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( size_t cell ) mutable
	{
		const float rho = 1.f;
		const float ux = 0.f;
		const float uy = 0.f;
		const float uz = 0.f;
		float feq[27];
		getFeq(rho, ux, uy, uz, feq);
		for ( size_t i = 0; i < 27; i++ ) fArrayView( i, cell ) = feq[i];
	};
	size_t start = 0;
	size_t end = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
