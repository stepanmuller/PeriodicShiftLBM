#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/ndarray/SizesHolder.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>

using ArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using VectorType = TNL::Containers::Vector< float, TNL::Devices::Cuda, size_t >;

using NDArray2DType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;
												
using NDArray4DType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0, 0, 0>,
												std::index_sequence< 0, 1, 2, 3 >,
												TNL::Devices::Cuda >;

const int size = 20000000;
const int iterations = 100;

void benchmarkArray( 
			ArrayType& f0Array, ArrayType& f1Array, ArrayType& f2Array, ArrayType& f3Array, ArrayType& f4Array, 
			ArrayType& f5Array, ArrayType& f6Array, ArrayType& f7Array, ArrayType& f8Array, ArrayType& f9Array, 
			ArrayType& f10Array, ArrayType& f11Array, ArrayType& f12Array, ArrayType& f13Array, ArrayType& f14Array, 
			ArrayType& f15Array, ArrayType& f16Array, ArrayType& f17Array, ArrayType& f18Array, ArrayType& f19Array, 
			ArrayType& f20Array, ArrayType& f21Array, ArrayType& f22Array, ArrayType& f23Array, ArrayType& f24Array, 
			ArrayType& f25Array, ArrayType& f26Array
			)
{
	auto f0ArrayView = f0Array.getView();
	auto f1ArrayView = f1Array.getView();
	auto f2ArrayView = f2Array.getView();
	auto f3ArrayView = f3Array.getView();
	auto f4ArrayView = f4Array.getView();
	auto f5ArrayView = f5Array.getView();
	auto f6ArrayView = f6Array.getView();
	auto f7ArrayView = f7Array.getView();
	auto f8ArrayView = f8Array.getView();
	auto f9ArrayView = f9Array.getView();
	auto f10ArrayView = f10Array.getView();
	auto f11ArrayView = f11Array.getView();
	auto f12ArrayView = f12Array.getView();
	auto f13ArrayView = f13Array.getView();
	auto f14ArrayView = f14Array.getView();
	auto f15ArrayView = f15Array.getView();
	auto f16ArrayView = f16Array.getView();
	auto f17ArrayView = f17Array.getView();
	auto f18ArrayView = f18Array.getView();
	auto f19ArrayView = f19Array.getView();
	auto f20ArrayView = f20Array.getView();
	auto f21ArrayView = f21Array.getView();
	auto f22ArrayView = f22Array.getView();
	auto f23ArrayView = f23Array.getView();
	auto f24ArrayView = f24Array.getView();
	auto f25ArrayView = f25Array.getView();
	auto f26ArrayView = f26Array.getView();
	
	auto benchmarkLambda = [=] __cuda_callable__ (size_t i) mutable
	{
		float f0 = f0ArrayView[i] + 1.f;
		float f1 = f1ArrayView[i] + 1.f;
		float f2 = f2ArrayView[i] + 1.f;
		float f3 = f3ArrayView[i] + 1.f;
		float f4 = f4ArrayView[i] + 1.f;
		float f5 = f5ArrayView[i] + 1.f;
		float f6 = f6ArrayView[i] + 1.f;
		float f7 = f7ArrayView[i] + 1.f;
		float f8 = f8ArrayView[i] + 1.f;
		float f9 = f9ArrayView[i] + 1.f;
		float f10 = f10ArrayView[i] + 1.f;
		float f11 = f11ArrayView[i] + 1.f;
		float f12 = f12ArrayView[i] + 1.f;
		float f13 = f13ArrayView[i] + 1.f;
		float f14 = f14ArrayView[i] + 1.f;
		float f15 = f15ArrayView[i] + 1.f;
		float f16 = f16ArrayView[i] + 1.f;
		float f17 = f17ArrayView[i] + 1.f;
		float f18 = f18ArrayView[i] + 1.f;
		float f19 = f19ArrayView[i] + 1.f;
		float f20 = f20ArrayView[i] + 1.f;
		float f21 = f21ArrayView[i] + 1.f;
		float f22 = f22ArrayView[i] + 1.f;
		float f23 = f23ArrayView[i] + 1.f;
		float f24 = f24ArrayView[i] + 1.f;
		float f25 = f25ArrayView[i] + 1.f;
		float f26 = f26ArrayView[i] + 1.f;
		
		f0ArrayView[i] = f0;
		f1ArrayView[i] = f1;
		f2ArrayView[i] = f2;
		f3ArrayView[i] = f3;
		f4ArrayView[i] = f4;
		f5ArrayView[i] = f5;
		f6ArrayView[i] = f6;
		f7ArrayView[i] = f7;
		f8ArrayView[i] = f8;
		f9ArrayView[i] = f9;
		f10ArrayView[i] = f10;
		f11ArrayView[i] = f11;
		f12ArrayView[i] = f12;
		f13ArrayView[i] = f13;
		f14ArrayView[i] = f14;
		f15ArrayView[i] = f15;
		f16ArrayView[i] = f16;
		f17ArrayView[i] = f17;
		f18ArrayView[i] = f18;
		f19ArrayView[i] = f19;
		f20ArrayView[i] = f20;
		f21ArrayView[i] = f21;
		f22ArrayView[i] = f22;
		f23ArrayView[i] = f23;
		f24ArrayView[i] = f24;
		f25ArrayView[i] = f25;
		f26ArrayView[i] = f26;

	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, size, benchmarkLambda );
}

void benchmarkVector( 
			VectorType& f0Vector, VectorType& f1Vector, VectorType& f2Vector, VectorType& f3Vector, VectorType& f4Vector, 
			VectorType& f5Vector, VectorType& f6Vector, VectorType& f7Vector, VectorType& f8Vector, VectorType& f9Vector, 
			VectorType& f10Vector, VectorType& f11Vector, VectorType& f12Vector, VectorType& f13Vector, VectorType& f14Vector, 
			VectorType& f15Vector, VectorType& f16Vector, VectorType& f17Vector, VectorType& f18Vector, VectorType& f19Vector, 
			VectorType& f20Vector, VectorType& f21Vector, VectorType& f22Vector, VectorType& f23Vector, VectorType& f24Vector, 
			VectorType& f25Vector, VectorType& f26Vector
			)
{
	auto f0VectorView = f0Vector.getView();
	auto f1VectorView = f1Vector.getView();
	auto f2VectorView = f2Vector.getView();
	auto f3VectorView = f3Vector.getView();
	auto f4VectorView = f4Vector.getView();
	auto f5VectorView = f5Vector.getView();
	auto f6VectorView = f6Vector.getView();
	auto f7VectorView = f7Vector.getView();
	auto f8VectorView = f8Vector.getView();
	auto f9VectorView = f9Vector.getView();
	auto f10VectorView = f10Vector.getView();
	auto f11VectorView = f11Vector.getView();
	auto f12VectorView = f12Vector.getView();
	auto f13VectorView = f13Vector.getView();
	auto f14VectorView = f14Vector.getView();
	auto f15VectorView = f15Vector.getView();
	auto f16VectorView = f16Vector.getView();
	auto f17VectorView = f17Vector.getView();
	auto f18VectorView = f18Vector.getView();
	auto f19VectorView = f19Vector.getView();
	auto f20VectorView = f20Vector.getView();
	auto f21VectorView = f21Vector.getView();
	auto f22VectorView = f22Vector.getView();
	auto f23VectorView = f23Vector.getView();
	auto f24VectorView = f24Vector.getView();
	auto f25VectorView = f25Vector.getView();
	auto f26VectorView = f26Vector.getView();
	
	auto benchmarkLambda = [=] __cuda_callable__ (size_t i) mutable
	{
		float f0 = f0VectorView[i] + 1.f;
		float f1 = f1VectorView[i] + 1.f;
		float f2 = f2VectorView[i] + 1.f;
		float f3 = f3VectorView[i] + 1.f;
		float f4 = f4VectorView[i] + 1.f;
		float f5 = f5VectorView[i] + 1.f;
		float f6 = f6VectorView[i] + 1.f;
		float f7 = f7VectorView[i] + 1.f;
		float f8 = f8VectorView[i] + 1.f;
		float f9 = f9VectorView[i] + 1.f;
		float f10 = f10VectorView[i] + 1.f;
		float f11 = f11VectorView[i] + 1.f;
		float f12 = f12VectorView[i] + 1.f;
		float f13 = f13VectorView[i] + 1.f;
		float f14 = f14VectorView[i] + 1.f;
		float f15 = f15VectorView[i] + 1.f;
		float f16 = f16VectorView[i] + 1.f;
		float f17 = f17VectorView[i] + 1.f;
		float f18 = f18VectorView[i] + 1.f;
		float f19 = f19VectorView[i] + 1.f;
		float f20 = f20VectorView[i] + 1.f;
		float f21 = f21VectorView[i] + 1.f;
		float f22 = f22VectorView[i] + 1.f;
		float f23 = f23VectorView[i] + 1.f;
		float f24 = f24VectorView[i] + 1.f;
		float f25 = f25VectorView[i] + 1.f;
		float f26 = f26VectorView[i] + 1.f;
		
		f0VectorView[i] = f0;
		f1VectorView[i] = f1;
		f2VectorView[i] = f2;
		f3VectorView[i] = f3;
		f4VectorView[i] = f4;
		f5VectorView[i] = f5;
		f6VectorView[i] = f6;
		f7VectorView[i] = f7;
		f8VectorView[i] = f8;
		f9VectorView[i] = f9;
		f10VectorView[i] = f10;
		f11VectorView[i] = f11;
		f12VectorView[i] = f12;
		f13VectorView[i] = f13;
		f14VectorView[i] = f14;
		f15VectorView[i] = f15;
		f16VectorView[i] = f16;
		f17VectorView[i] = f17;
		f18VectorView[i] = f18;
		f19VectorView[i] = f19;
		f20VectorView[i] = f20;
		f21VectorView[i] = f21;
		f22VectorView[i] = f22;
		f23VectorView[i] = f23;
		f24VectorView[i] = f24;
		f25VectorView[i] = f25;
		f26VectorView[i] = f26;

	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, size, benchmarkLambda );
}

void benchmarkNDArray2D( 
			NDArray2DType& fNDArray2D
			)
{
	auto fNDArray2DView = fNDArray2D.getView();
	
	auto benchmarkLambda = [=] __cuda_callable__ (size_t i) mutable
	{		
		float f0  = fNDArray2DView(0, i) + 1.f;
		float f1  = fNDArray2DView(1, i) + 1.f;
		float f2  = fNDArray2DView(2, i) + 1.f;
		float f3  = fNDArray2DView(3, i) + 1.f;
		float f4  = fNDArray2DView(4, i) + 1.f;
		float f5  = fNDArray2DView(5, i) + 1.f;
		float f6  = fNDArray2DView(6, i) + 1.f;
		float f7  = fNDArray2DView(7, i) + 1.f;
		float f8  = fNDArray2DView(8, i) + 1.f;
		float f9  = fNDArray2DView(9, i) + 1.f;
		float f10 = fNDArray2DView(10, i) + 1.f;
		float f11 = fNDArray2DView(11, i) + 1.f;
		float f12 = fNDArray2DView(12, i) + 1.f;
		float f13 = fNDArray2DView(13, i) + 1.f;
		float f14 = fNDArray2DView(14, i) + 1.f;
		float f15 = fNDArray2DView(15, i) + 1.f;
		float f16 = fNDArray2DView(16, i) + 1.f;
		float f17 = fNDArray2DView(17, i) + 1.f;
		float f18 = fNDArray2DView(18, i) + 1.f;
		float f19 = fNDArray2DView(19, i) + 1.f;
		float f20 = fNDArray2DView(20, i) + 1.f;
		float f21 = fNDArray2DView(21, i) + 1.f;
		float f22 = fNDArray2DView(22, i) + 1.f;
		float f23 = fNDArray2DView(23, i) + 1.f;
		float f24 = fNDArray2DView(24, i) + 1.f;
		float f25 = fNDArray2DView(25, i) + 1.f;
		float f26 = fNDArray2DView(26, i) + 1.f;
		
		fNDArray2DView(0, i)  = f0;
		fNDArray2DView(1, i)  = f1;
		fNDArray2DView(2, i)  = f2;
		fNDArray2DView(3, i)  = f3;
		fNDArray2DView(4, i)  = f4;
		fNDArray2DView(5, i)  = f5;
		fNDArray2DView(6, i)  = f6;
		fNDArray2DView(7, i)  = f7;
		fNDArray2DView(8, i)  = f8;
		fNDArray2DView(9, i)  = f9;
		fNDArray2DView(10, i) = f10;
		fNDArray2DView(11, i) = f11;
		fNDArray2DView(12, i) = f12;
		fNDArray2DView(13, i) = f13;
		fNDArray2DView(14, i) = f14;
		fNDArray2DView(15, i) = f15;
		fNDArray2DView(16, i) = f16;
		fNDArray2DView(17, i) = f17;
		fNDArray2DView(18, i) = f18;
		fNDArray2DView(19, i) = f19;
		fNDArray2DView(20, i) = f20;
		fNDArray2DView(21, i) = f21;
		fNDArray2DView(22, i) = f22;
		fNDArray2DView(23, i) = f23;
		fNDArray2DView(24, i) = f24;
		fNDArray2DView(25, i) = f25;
		fNDArray2DView(26, i) = f26;

	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, size, benchmarkLambda );
}

void benchmarkNDArray4D( 
			NDArray4DType& fNDArray4D
			)
{
	auto fNDArray4DView = fNDArray4D.getView();
	
	auto benchmarkLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 3, int >& tripleIndex) mutable
	{	
		const int i = tripleIndex.x();
		const int j = tripleIndex.y();
		const int k = tripleIndex.z();
		
		float f0  = fNDArray4DView(0, i, j, k) + 1.f;
		float f1  = fNDArray4DView(1, i, j, k) + 1.f;
		float f2  = fNDArray4DView(2, i, j, k) + 1.f;
		float f3  = fNDArray4DView(3, i, j, k) + 1.f;
		float f4  = fNDArray4DView(4, i, j, k) + 1.f;
		float f5  = fNDArray4DView(5, i, j, k) + 1.f;
		float f6  = fNDArray4DView(6, i, j, k) + 1.f;
		float f7  = fNDArray4DView(7, i, j, k) + 1.f;
		float f8  = fNDArray4DView(8, i, j, k) + 1.f;
		float f9  = fNDArray4DView(9, i, j, k) + 1.f;
		float f10 = fNDArray4DView(10, i, j, k) + 1.f;
		float f11 = fNDArray4DView(11, i, j, k) + 1.f;
		float f12 = fNDArray4DView(12, i, j, k) + 1.f;
		float f13 = fNDArray4DView(13, i, j, k) + 1.f;
		float f14 = fNDArray4DView(14, i, j, k) + 1.f;
		float f15 = fNDArray4DView(15, i, j, k) + 1.f;
		float f16 = fNDArray4DView(16, i, j, k) + 1.f;
		float f17 = fNDArray4DView(17, i, j, k) + 1.f;
		float f18 = fNDArray4DView(18, i, j, k) + 1.f;
		float f19 = fNDArray4DView(19, i, j, k) + 1.f;
		float f20 = fNDArray4DView(20, i, j, k) + 1.f;
		float f21 = fNDArray4DView(21, i, j, k) + 1.f;
		float f22 = fNDArray4DView(22, i, j, k) + 1.f;
		float f23 = fNDArray4DView(23, i, j, k) + 1.f;
		float f24 = fNDArray4DView(24, i, j, k) + 1.f;
		float f25 = fNDArray4DView(25, i, j, k) + 1.f;
		float f26 = fNDArray4DView(26, i, j, k) + 1.f;
		
		fNDArray4DView(0, i, j, k)  = f0;
		fNDArray4DView(1, i, j, k)  = f1;
		fNDArray4DView(2, i, j, k)  = f2;
		fNDArray4DView(3, i, j, k)  = f3;
		fNDArray4DView(4, i, j, k)  = f4;
		fNDArray4DView(5, i, j, k)  = f5;
		fNDArray4DView(6, i, j, k)  = f6;
		fNDArray4DView(7, i, j, k)  = f7;
		fNDArray4DView(8, i, j, k)  = f8;
		fNDArray4DView(9, i, j, k)  = f9;
		fNDArray4DView(10, i, j, k) = f10;
		fNDArray4DView(11, i, j, k) = f11;
		fNDArray4DView(12, i, j, k) = f12;
		fNDArray4DView(13, i, j, k) = f13;
		fNDArray4DView(14, i, j, k) = f14;
		fNDArray4DView(15, i, j, k) = f15;
		fNDArray4DView(16, i, j, k) = f16;
		fNDArray4DView(17, i, j, k) = f17;
		fNDArray4DView(18, i, j, k) = f18;
		fNDArray4DView(19, i, j, k) = f19;
		fNDArray4DView(20, i, j, k) = f20;
		fNDArray4DView(21, i, j, k) = f21;
		fNDArray4DView(22, i, j, k) = f22;
		fNDArray4DView(23, i, j, k) = f23;
		fNDArray4DView(24, i, j, k) = f24;
		fNDArray4DView(25, i, j, k) = f25;
		fNDArray4DView(26, i, j, k) = f26;

	};
	TNL::Containers::StaticArray< 3, int > start{ 0, 0, 0 };
	TNL::Containers::StaticArray< 3, int > end{ size, 1, 1 };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, benchmarkLambda );
}

int main(int argc, char **argv)
{
	// --------------------------------------------------------- //
	// ------------------------Array---------------------------- //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		ArrayType f0Array( size, 1.f );
		ArrayType f1Array( size, 1.f );
		ArrayType f2Array( size, 1.f );
		ArrayType f3Array( size, 1.f );
		ArrayType f4Array( size, 1.f );
		ArrayType f5Array( size, 1.f );
		ArrayType f6Array( size, 1.f );
		ArrayType f7Array( size, 1.f );
		ArrayType f8Array( size, 1.f );
		ArrayType f9Array( size, 1.f );
		ArrayType f10Array( size, 1.f );
		ArrayType f11Array( size, 1.f );
		ArrayType f12Array( size, 1.f );
		ArrayType f13Array( size, 1.f );
		ArrayType f14Array( size, 1.f );
		ArrayType f15Array( size, 1.f );
		ArrayType f16Array( size, 1.f );
		ArrayType f17Array( size, 1.f );
		ArrayType f18Array( size, 1.f );
		ArrayType f19Array( size, 1.f );
		ArrayType f20Array( size, 1.f );
		ArrayType f21Array( size, 1.f );
		ArrayType f22Array( size, 1.f );
		ArrayType f23Array( size, 1.f );
		ArrayType f24Array( size, 1.f );
		ArrayType f25Array( size, 1.f );
		ArrayType f26Array( size, 1.f );  
		#ifdef __CUDACC__
		std::cout << "starting Array benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
			benchmarkArray(
							f0Array, f1Array, f2Array, f3Array, f4Array, f5Array, f6Array, f7Array, f8Array, 
							f9Array, f10Array, f11Array, f12Array, f13Array, f14Array, f15Array, f16Array, 
							f17Array, f18Array, f19Array, f20Array, f21Array, f22Array, f23Array, f24Array, 
							f25Array, f26Array
							);
		timer.stop();
		#endif
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}
	
	// --------------------------------------------------------- //
	// ------------------------Vector--------------------------- //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		VectorType f0Vector( size, 1.f );
		VectorType f1Vector( size, 1.f );
		VectorType f2Vector( size, 1.f );
		VectorType f3Vector( size, 1.f );
		VectorType f4Vector( size, 1.f );
		VectorType f5Vector( size, 1.f );
		VectorType f6Vector( size, 1.f );
		VectorType f7Vector( size, 1.f );
		VectorType f8Vector( size, 1.f );
		VectorType f9Vector( size, 1.f );
		VectorType f10Vector( size, 1.f );
		VectorType f11Vector( size, 1.f );
		VectorType f12Vector( size, 1.f );
		VectorType f13Vector( size, 1.f );
		VectorType f14Vector( size, 1.f );
		VectorType f15Vector( size, 1.f );
		VectorType f16Vector( size, 1.f );
		VectorType f17Vector( size, 1.f );
		VectorType f18Vector( size, 1.f );
		VectorType f19Vector( size, 1.f );
		VectorType f20Vector( size, 1.f );
		VectorType f21Vector( size, 1.f );
		VectorType f22Vector( size, 1.f );
		VectorType f23Vector( size, 1.f );
		VectorType f24Vector( size, 1.f );
		VectorType f25Vector( size, 1.f );
		VectorType f26Vector( size, 1.f );  
		#ifdef __CUDACC__
		std::cout << "starting Vector benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
			benchmarkVector(
							f0Vector, f1Vector, f2Vector, f3Vector, f4Vector, f5Vector, f6Vector, f7Vector, f8Vector, 
							f9Vector, f10Vector, f11Vector, f12Vector, f13Vector, f14Vector, f15Vector, f16Vector, 
							f17Vector, f18Vector, f19Vector, f20Vector, f21Vector, f22Vector, f23Vector, f24Vector, 
							f25Vector, f26Vector
							);
		timer.stop();
		#endif
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}
	
	// --------------------------------------------------------- //
	// ------------------------NDArray2D------------------------ //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		NDArray2DType fNDArray2D;
		fNDArray2D.setSizes( 27, size );
		fNDArray2D.setValue( 1.0f ); 
		#ifdef __CUDACC__
		std::cout << "starting NDArray2D (27 x size) benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
			benchmarkNDArray2D( fNDArray2D );
		timer.stop();
		#endif
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}	
	
	// --------------------------------------------------------- //
	// ------------------------NDArray4D------------------------ //
	// --------------------------------------------------------- //
	for (int scope=0; scope<1; scope++)
	{
		NDArray4DType fNDArray4D;
		fNDArray4D.setSizes( 27, size, 1, 1 );
		fNDArray4D.setValue( 1.0f ); 
		#ifdef __CUDACC__
		std::cout << "starting NDArray4D (27 x size x 1 x 1) benchmark" << std::endl;
		TNL::Timer timer;
		timer.start();
		for (int i=0; i<iterations; i++)
			benchmarkNDArray4D( fNDArray4D );
		timer.stop();
		#endif
		auto totalTime = timer.getRealTime();
		std::cout << "this took " << totalTime << " s" << std::endl;
		std::cout << std::endl;
	}	

	return EXIT_SUCCESS;
}
