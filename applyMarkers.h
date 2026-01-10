#include "config.h"
#include "cellFunctions.h"

void applyMarkers( MarkerStruct& Marker )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getView();
	auto givenRhoMarkerArrayView = Marker.givenRhoArray.getView();
	auto givenUxUyUzMarkerArrayView = Marker.givenUxUyUzArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 3, int >& tripleIndex) mutable
	{
		const size_t i = tripleIndex.x();
		const size_t j = tripleIndex.y();
		const size_t k = tripleIndex.z();
		size_t cell = convertIndex(i, j, k);
		if ( j>=boxStartJ && j<=boxEndJ && k>=boxStartK && k<=boxEndK ) 
		{
			bouncebackMarkerArrayView[cell] = 1;
			return;
		}
		if ( k==0  || i == 0 || i == cellCountX - 1 || j == 0 || j == cellCountY - 1 ) 
		{
			givenUxUyUzMarkerArrayView[cell] = 1;
		}
		if ( k==cellCountZ - 1 )
		{
			givenRhoMarkerArrayView[cell] = 1;
		}
		if ( i > 0 && i < cellCountX - 1 && j > 0 && j < cellCountY - 1 && k > 0 && k < cellCountZ - 1 )
		{
			fluidMarkerArrayView[cell] = 1;
		}
	};
	TNL::Containers::StaticArray< 3, size_t > start{ 0, 0, 0 };
	TNL::Containers::StaticArray< 3, size_t > end{ cellCountX, cellCountY, cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
