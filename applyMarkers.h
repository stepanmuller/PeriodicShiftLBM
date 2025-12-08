#include "config.h"
#include "convertIndex.h"

void applyMarkers( MarkerStruct& Marker, ArrayType& uzArray )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getView();
	auto inletMarkerArrayView = Marker.inletArray.getView();
	auto outletMarkerArrayView = Marker.outletArray.getView();
	
	auto uzArrayView = uzArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 3, int >& tripleIndex) mutable
	{
		const int i = tripleIndex.x();
		const int j = tripleIndex.y();
		const int k = tripleIndex.z();
		size_t cell = convertIndex(i, j, k);
		if ( j>=350 && j<=450 && k>=250 && k<=350 ) 
		{
			bouncebackMarkerArrayView[cell] = 1;
		}
		else if ( i==0 || i==cellCountX-1 || j==0 || j==cellCountY-1 ) 
		{
			bouncebackMarkerArrayView[cell] = 1;
		}
		else if ( k==0 ) 
		{
			inletMarkerArrayView[cell] = 1;
			uzArrayView[cell] = uzInlet;
		}
		else if ( k==cellCountZ-1 ) 
		{
			outletMarkerArrayView[cell] = 1;
		}
		else
		{
			fluidMarkerArrayView[cell] = 1;
		}
	};
	TNL::Containers::StaticArray< 3, int > start{ 0, 0, 0 };
	TNL::Containers::StaticArray< 3, int > end{ cellCountX, cellCountY, cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}
