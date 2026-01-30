#include <iostream>
#include <cmath>
#include <fstream> 
#include <cstdlib>
#include <limits>
#include <string>

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/AtomicOperations.h>
#include <TNL/Algorithms/reduce.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>

using BoolArray3DType = TNL::Containers::NDArray< bool, 
												TNL::Containers::SizesHolder< size_t, 0, 0, 0 >,
												std::index_sequence< 0, 1, 2 >,
												TNL::Devices::Cuda >;
												
using IntArrayType = TNL::Containers::Array< int, TNL::Devices::Cuda, size_t >;

using IntArray2DType = TNL::Containers::NDArray< int, 
												TNL::Containers::SizesHolder< size_t, 0, 0 >,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;

using FloatArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using FloatArrayTypeCPU = TNL::Containers::Array< float, TNL::Devices::Host, size_t >;

using FloatArray2DType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< size_t, 0, 0 >,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;
using FloatArray2DTypeCPU = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< size_t, 0, 0 >,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Host >;
												
using LongIntArrayType = TNL::Containers::Array< long int, TNL::Devices::Cuda, size_t >;
using LongIntArrayTypeCPU = TNL::Containers::Array< long int, TNL::Devices::Host, size_t >;

using LongIntPairType = TNL::Containers::StaticArray< 2, long int >;											
using LongIntTripleType = TNL::Containers::StaticArray< 3, long int >;

struct FStruct { LongIntArrayType shifter; FloatArray2DType fArray; };
struct FStructCPU { LongIntArrayTypeCPU shifter; FloatArray2DTypeCPU fArray; };

struct SectionCutStruct { FloatArray2DType rhoArray; FloatArray2DType uxArray; FloatArray2DType uyArray; FloatArray2DType uzArray; FloatArray2DType maskArray; };
struct SectionCutStructCPU { FloatArray2DTypeCPU rhoArray; FloatArray2DTypeCPU uxArray; FloatArray2DTypeCPU uyArray; FloatArray2DTypeCPU uzArray; FloatArray2DTypeCPU maskArray; };

struct MarkerStruct { BoolArray3DType fluidArray; BoolArray3DType bouncebackArray; BoolArray3DType givenUxUyUzArray; BoolArray3DType givenRhoArray;  };

struct STLStruct { 	FloatArrayType axArray; FloatArrayType ayArray; FloatArrayType azArray; 
					FloatArrayType bxArray; FloatArrayType byArray; FloatArrayType bzArray; 
					FloatArrayType cxArray; FloatArrayType cyArray; FloatArrayType czArray; 
					float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; 
					float ox; float oy; float oz; int triangleCount; }; // ox, oy, oz is the position of STL origin in global coordinates

struct STLStructCPU { 	FloatArrayTypeCPU axArray; FloatArrayTypeCPU ayArray; FloatArrayTypeCPU azArray; 
						FloatArrayTypeCPU bxArray; FloatArrayTypeCPU byArray; FloatArrayTypeCPU bzArray; 
						FloatArrayTypeCPU cxArray; FloatArrayTypeCPU cyArray; FloatArrayTypeCPU czArray; 
						float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; 
						float ox; float oy; float oz; int triangleCount; }; // ox, oy, oz is the position of STL origin in global coordinates
								
struct InfoStruct { float res; int cellCountX; int cellCountY; int cellCountZ; int iterationsFinished; 
					float rhoNominalPhys; float soundspeedPhys; float dtPhys; }; 
