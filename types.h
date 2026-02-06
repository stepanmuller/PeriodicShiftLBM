#include <iostream>
#include <cmath>
#include <fstream> 
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/AtomicOperations.h>
#include <TNL/Algorithms/reduce.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>

using BoolArrayType = TNL::Containers::Array< bool, TNL::Devices::Cuda, size_t >;
												
using IntArrayType = TNL::Containers::Array< int, TNL::Devices::Cuda, size_t >;
using IntArrayConstViewType = TNL::Containers::ArrayView< const int, TNL::Devices::Cuda, size_t >;
using IntArrayTypeCPU = TNL::Containers::Array< int, TNL::Devices::Host, size_t >;

using IntArray2DType = TNL::Containers::NDArray< int, 
												TNL::Containers::SizesHolder< size_t, 0, 0 >,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;
												
using IntArray3DType = TNL::Containers::NDArray< int, 
												TNL::Containers::SizesHolder< size_t, 0, 0, 0 >,
												std::index_sequence< 0, 1, 2 >,
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

using IntPairType = TNL::Containers::StaticArray< 2, int >;											
using IntTripleType = TNL::Containers::StaticArray< 3, int >;

struct FStruct { IntArrayType shifter; FloatArray2DType fArray; };
struct FStructCPU { IntArrayTypeCPU shifter; FloatArray2DTypeCPU fArray; };

struct SectionCutStruct { FloatArray2DType rhoArray; FloatArray2DType uxArray; FloatArray2DType uyArray; FloatArray2DType uzArray; FloatArray2DType markerArray; };
struct SectionCutStructCPU { FloatArray2DTypeCPU rhoArray; FloatArray2DTypeCPU uxArray; FloatArray2DTypeCPU uyArray; FloatArray2DTypeCPU uzArray; FloatArray2DTypeCPU markerArray; };		
								
struct InfoStruct { float res = 1.f; int cellCountX; int cellCountY; int cellCountZ; int cellCount;
					float rhoNominalPhys = 1.f; float soundspeedPhys = 1.f; float dtPhys = 1.f; 
					float iRegulator = 0.f; }; 
					
struct STLStructCPU { 	FloatArrayTypeCPU axArray; FloatArrayTypeCPU ayArray; FloatArrayTypeCPU azArray; 
						FloatArrayTypeCPU bxArray; FloatArrayTypeCPU byArray; FloatArrayTypeCPU bzArray; 
						FloatArrayTypeCPU cxArray; FloatArrayTypeCPU cyArray; FloatArrayTypeCPU czArray; 
						float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; 
						float ox; float oy; float oz; int triangleCount; }; // ox, oy, oz is the position of STL origin in global coordinates

struct STLStruct { 	FloatArrayType axArray; FloatArrayType ayArray; FloatArrayType azArray; 
					FloatArrayType bxArray; FloatArrayType byArray; FloatArrayType bzArray; 
					FloatArrayType cxArray; FloatArrayType cyArray; FloatArrayType czArray; 
					float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; 
					float ox; float oy; float oz; int triangleCount; // ox, oy, oz is the position of STL origin in global coordinates
					STLStruct() = default;
					// Constructor copies data from STLStructCPU
					STLStruct( const STLStructCPU& STLCPU )
					{
						axArray = STLCPU.axArray;
						ayArray = STLCPU.ayArray;
						azArray = STLCPU.azArray; 
						bxArray = STLCPU.bxArray;
						byArray = STLCPU.byArray;
						bzArray = STLCPU.bzArray; 
						cxArray = STLCPU.cxArray;
						cyArray = STLCPU.cyArray;
						czArray = STLCPU.czArray; 
						xmin = STLCPU.xmin;
						ymin = STLCPU.ymin;
						zmin = STLCPU.zmin; 
						xmax = STLCPU.xmax;
						ymax = STLCPU.ymax;
						zmax = STLCPU.zmax; 
						ox = STLCPU.ox;
						oy = STLCPU.oy;
						oz = STLCPU.oz;
						triangleCount = STLCPU.triangleCount;
					}
				};	
