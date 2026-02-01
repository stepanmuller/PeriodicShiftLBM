// Perform periodic shift streaming by adjusting F.shifter

// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

IntArrayType cxArray{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
IntArrayType cyArray{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
IntArrayType czArray{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

void applyStreaming( FStruct& F, InfoStruct &Info )
{
	auto shifterView = F.shifter.getView();
	auto cxArrayView = cxArray.getConstView();
	auto cyArrayView = cyArray.getConstView();
	auto czArrayView = czArray.getConstView();
	auto directionLambda = [=] __cuda_callable__ ( const int direction ) mutable
	{
		int shift = shifterView[direction];
		shift -= cxArrayView[direction];
		shift -= Info.cellCountX * cyArrayView[direction];
		shift -= Info.cellCountY * Info.cellCountX * czArrayView[direction];
		if ( shift < 0 ) shift += Info.cellCount;
		else if ( shift >= Info.cellCount ) shift -= Info.cellCount;
		shifterView[direction] = shift;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, 27, directionLambda );
}

