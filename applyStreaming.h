// Perform periodic shift streaming by adjusting F.shifter

// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

using DirectionType = TNL::Containers::Array< long, TNL::Devices::Cuda, size_t >;
DirectionType directionCx{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
DirectionType directionCy{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
DirectionType directionCz{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

void applyStreaming( DistributionStruct& F, CellCountStruct &cellCount )
{
	auto shifterView = F.shifter.getView();
	auto directionCxView = directionCx.getConstView();
	auto directionCyView = directionCy.getConstView();
	auto directionCzView = directionCz.getConstView();
	
	auto streamLambda = [=] __cuda_callable__ (size_t direction) mutable
	{
		long shift = static_cast<long>(shifterView[direction]);
		shift -= directionCxView[direction];
		shift -= static_cast<long>(cellCount.nx) * directionCyView[direction];
		shift -= static_cast<long>(cellCount.ny * cellCount.nx) * (directionCzView[direction]);
		if (shift < 0){ shift += cellCount.n; }
		else if (shift >= static_cast<long>(cellCount.n)) { shift -= cellCount.n; }
		shifterView[direction] = static_cast<size_t>(shift);
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, 27, streamLambda );
}
