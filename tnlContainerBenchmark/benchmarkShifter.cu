#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Timer.h>

using FloatArray2DType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;

using IntArrayType = TNL::Containers::Array< int, TNL::Devices::Cuda, size_t >;

constexpr int cellCountX = 2000;
constexpr int cellCountY = 1000;
constexpr int cellCountZ = 10;
constexpr int cellCount = cellCountX * cellCountY * cellCountZ;

constexpr int iterations = 100;

IntArrayType cxArray{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
IntArrayType cyArray{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
IntArrayType czArray{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

void applyStreaming( IntArrayType& shifter )
{
	auto shifterView = shifter.getView();
	auto cxArrayView = cxArray.getConstView();
	auto cyArrayView = cyArray.getConstView();
	auto czArrayView = czArray.getConstView();
	auto directionLambda = [=] __cuda_callable__ (int direction) mutable
	{
		int shift = shifterView[direction];
		shift -= cxArrayView[direction];
		shift -= cellCountX * cyArrayView[direction];
		shift -= cellCountY * cellCountX * czArrayView[direction];
		if (shift < 0) shift += cellCount;
		else if (shift >= cellCount) shift -= cellCount;
		shifterView[direction] = shift;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, 27, directionLambda );
}

__cuda_callable__ void fTangle( float (&f)[27] )
{
	for ( int direction = 0; direction < 27; direction++ )
	{
		for ( int add = 0; add < 27; add++ )
		{
			f[direction] = f[direction] + f[add] * (( add % 2 ) * 2.f - 1.f) + 1.f;
		}
	}
}

void benchmarkNDArray2DFlat( FloatArray2DType& fArray )
{
	auto fArrayView = fArray.getView();
	
	auto benchmarkLambda = [=] __cuda_callable__ ( int cell ) mutable
	{		
		float f[27];
		for (int i = 0; i < 27; i++)	f[i] = fArrayView( i, cell );
		fTangle( f );
		for (int i = 0; i < 27; i++)	fArrayView( i, cell ) = f[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, benchmarkLambda );
}

void benchmarkNDArray2DShifter( FloatArray2DType& fArray, IntArrayType& shifter )
{
	auto fArrayView = fArray.getView();
	auto shifterView = shifter.getConstView();
	
	auto benchmarkLambda = [=] __cuda_callable__ ( int cell ) mutable
	{		
		int shiftedIndex[27];
		for (int i = 0; i < 27; i++) 
		{
			const int shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
		}
		float f[27];
		for (int i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		fTangle( f );
		for (int i = 0; i < 27; i++)	fArrayView(i, shiftedIndex[i]) = f[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, benchmarkLambda );
}

int main(int argc, char **argv)
{
	// --------------------------------------------------------- //
	// ------------------------NDArray2D------------------------ //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		FloatArray2DType fArray;
		fArray.setSizes( 27, cellCount );
		fArray.setValue( 1.0f ); 
		std::cout << "starting NDArray2D (27 x cellCount) benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
		{
			benchmarkNDArray2DFlat( fArray );
		}
		timer.stop();
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}		
	// --------------------------------------------------------- //
	// ----------------NDArray2D with shifter------------------- //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		FloatArray2DType fArray;
		fArray.setSizes( 27, cellCount );
		fArray.setValue( 1.0f ); 
		IntArrayType shifter = IntArrayType( 27, 0 );
		std::cout << "starting NDArray2D (27 x cellCount) with shifter benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
		{
			applyStreaming( shifter );
			benchmarkNDArray2DShifter( fArray, shifter );
		}
		timer.stop();
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}	

	return EXIT_SUCCESS;
}
