#include <iostream>
#include <cmath>
#include <fstream> 
#include <cstdlib>
#include <limits>

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>

using IndexArrayType = TNL::Containers::Array< size_t, TNL::Devices::Cuda, size_t >;
using IndexArrayTypeCPU = TNL::Containers::Array< size_t, TNL::Devices::Host, size_t >;

using FloatArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using FloatArrayTypeCPU = TNL::Containers::Array< float, TNL::Devices::Host, size_t >;

using DistributionArrayType = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Cuda >;
using DistributionArrayTypeCPU = TNL::Containers::NDArray< float, 
												TNL::Containers::SizesHolder< std::size_t, 0, 0>,
												std::index_sequence< 0, 1 >,
												TNL::Devices::Host >;

using MarkerArrayType = TNL::Containers::Array< bool, TNL::Devices::Cuda, size_t >;

struct DistributionStruct { IndexArrayType shifter; DistributionArrayType fArray; };
struct DistributionStructCPU { IndexArrayTypeCPU shifter; DistributionArrayTypeCPU fArray; };

struct MarkerStruct { MarkerArrayType fluidArray; MarkerArrayType bouncebackArray; MarkerArrayType givenUxUyUzArray; MarkerArrayType givenRhoArray;  };

struct STLArbeiterStructCPU { 	FloatArrayTypeCPU nxArray; FloatArrayTypeCPU nyArray; FloatArrayTypeCPU nzArray; 
								FloatArrayTypeCPU axArray; FloatArrayTypeCPU ayArray; FloatArrayTypeCPU azArray; 
								FloatArrayTypeCPU bxArray; FloatArrayTypeCPU byArray; FloatArrayTypeCPU bzArray; 
								FloatArrayTypeCPU cxArray; FloatArrayTypeCPU cyArray; FloatArrayTypeCPU czArray; 
								float xmin; float ymin; float zmin; float xmax; float ymax; float zmax; };
								
struct CellCountStruct { size_t nx; size_t ny; size_t nz; size_t n; float ox; float oy; float oz; }; // ox, oy, oz is the position of cell i,j,k = 0 in global coordinates

#include "applyMarkers.h"
#include "applyInitialization.h"

#include "STL/STLFunctions.h"

#include "applyStreaming.h"
#include "applyCollision.h"
#include "cellFunctions.h"

#include "boundaryConditions/applyBounceback.h"
#include "boundaryConditions/restoreRho.h"
#include "boundaryConditions/restoreUxUyUz.h"
#include "boundaryConditions/restoreRhoUxUyUz.h"
#include "boundaryConditions/applyMBBC.h"

#include "applyLocalCellUpdate.h"
