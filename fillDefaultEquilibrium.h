void fillDefaultEquilibrium( FloatArray2DType& fArray, InfoStruct &Info )
{
	std::cout << "Filling fArray with default equilibrium, rho = 1, uxyz = 0" << std::endl;
	auto fArrayView  = fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const TripleIndexType& tripleIndex ) mutable
	{
		const int iCell = tripleIndex.x();
		const int jCell = tripleIndex.y();
		const int kCell = tripleIndex.z();	
		size_t shiftedIndex[27];
		getShiftedIndex( iCell, jCell, kCell, shiftedIndex, Info );	
		const float rho = 1.f;
		const float ux = 0.f;
		const float uy = 0.f;
		const float uz = 0.f;
		float feq[27];
		getFeq( rho, ux, uy, uz, feq );
		for ( size_t direction = 0; direction < 27; direction++ ) fArrayView( direction, shiftedIndex[direction] ) = feq[direction];
	};
	TripleIndexType start{ 0, 0, 0 };
	TripleIndexType end{ Info.cellCountX, Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( start, end, cellLambda );
}
