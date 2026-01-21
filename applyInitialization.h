void applyInitialization( 	DistributionStruct& F, CellCountStruct &cellCount )
{
	auto shifterView = F.shifter.getConstView();
	auto fArrayView  = F.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount.n) { shiftedIndex[i] -= cellCount.n; }
		}
		
		const float rho = 1.f;
		const float ux = 0.f;
		const float uy = 0.f;
		const float uz = 0.f;
		
		float feq[27];
		getFeq(rho, ux, uy, uz, feq);
		
		for (size_t i = 0; i < 27; i++)	fArrayView(i, shiftedIndex[i]) = feq[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount.n, cellLambda );
}
