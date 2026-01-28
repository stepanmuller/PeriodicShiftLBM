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

using IndexArrayType = TNL::Containers::Array< size_t, TNL::Devices::Cuda, size_t >;
using IndexArrayTypeCPU = TNL::Containers::Array< size_t, TNL::Devices::Host, size_t >;

using IndexArray3DType = TNL::Containers::NDArray< size_t, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0, 0>,
												std::index_sequence< 0, 1, 2 >,
												TNL::Devices::Cuda >;

using FloatArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using FloatArrayTypeCPU = TNL::Containers::Array< float, TNL::Devices::Host, size_t >;

using FloatArray2DType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;
using FloatArray2DTypeCPU = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Host >;

using MarkerArrayType = TNL::Containers::Array< bool, TNL::Devices::Cuda, size_t >;

using CounterArray2DType = TNL::Containers::NDArray< int, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;

struct DistributionStruct { IndexArrayType shifter; FloatArray2DType fArray; };
struct DistributionStructCPU { IndexArrayTypeCPU shifter; FloatArray2DTypeCPU fArray; };

struct SectionCutStruct { FloatArray2DType rhoArray; FloatArray2DType uxArray; FloatArray2DType uyArray; FloatArray2DType uzArray; FloatArray2DType maskArray; };
struct SectionCutStructCPU { FloatArray2DTypeCPU rhoArray; FloatArray2DTypeCPU uxArray; FloatArray2DTypeCPU uyArray; FloatArray2DTypeCPU uzArray; FloatArray2DTypeCPU maskArray; };

struct MarkerStruct { MarkerArrayType fluidArray; MarkerArrayType bouncebackArray; MarkerArrayType givenUxUyUzArray; MarkerArrayType givenRhoArray;  };

struct STLArbeiterStruct { 	FloatArrayType axArray; FloatArrayType ayArray; FloatArrayType azArray; 
							FloatArrayType bxArray; FloatArrayType byArray; FloatArrayType bzArray; 
							FloatArrayType cxArray; FloatArrayType cyArray; FloatArrayType czArray; 
							float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; size_t triangleCount; };

struct STLArbeiterStructCPU { 	FloatArrayTypeCPU axArray; FloatArrayTypeCPU ayArray; FloatArrayTypeCPU azArray; 
								FloatArrayTypeCPU bxArray; FloatArrayTypeCPU byArray; FloatArrayTypeCPU bzArray; 
								FloatArrayTypeCPU cxArray; FloatArrayTypeCPU cyArray; FloatArrayTypeCPU czArray; 
								float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; size_t triangleCount; };
								
struct CellCountStruct { float res; size_t nx; size_t ny; size_t nz; size_t n; float ox; float oy; float oz; }; // ox, oy, oz is the position of cell i,j,k = 0 in global coordinates

#include "cellFunctions.h"
#include "applyStreaming.h"
#include "applyCollision.h"

#include "STLFunctions.h"

#include "applyInitialization.h"

#include "boundaryConditions/applyBounceback.h"
#include "boundaryConditions/restoreRho.h"
#include "boundaryConditions/restoreUxUyUz.h"
#include "boundaryConditions/restoreRhoUxUyUz.h"
#include "boundaryConditions/applyMBBC.h"

#include "applyLocalCellUpdate.h"
