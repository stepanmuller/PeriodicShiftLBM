#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>
#include <cmath>
#include <fstream> 
#include <cstdlib>
#include <limits>

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

#include "config.h"
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

int main(int argc, char **argv)
{
	std::cout << "Starting Periodic Shift LBM" << std::endl;
	
	STLArbeiterStructCPU STLArbeiterCPU;
	std::cout << "Initialization: Reading STL" << std::endl;
	readSTL( STLArbeiterCPU );
	
	DistributionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	F.fArray.setSizes( 27, cellCount );
	F.fArray.setValue( 1.0f );
	
	MarkerStruct Marker;
	Marker.fluidArray = MarkerArrayType( cellCount, 0);
	Marker.bouncebackArray = MarkerArrayType( cellCount, 0);
	Marker.givenRhoArray = MarkerArrayType( cellCount, 0);
	Marker.givenUxUyUzArray = MarkerArrayType( cellCount, 0);
	
	std::cout << "Initialization: Marking cells" << std::endl;
	applyMarkers(Marker);
	
	std::cout << "Initialization: Filling F" << std::endl;
	applyInitialization( F );
	
	#ifdef __CUDACC__
	std::cout << "Starting simulation" << std::endl;
	TNL::Timer timer;
	TNL::Timer lapTimer;
	timer.start();
	lapTimer.start();
	for (int i=0; i<iterationCount; i++)
	{
		applyStreaming( F );
		applyLocalCellUpdate( Marker, F );
		
		if (i%100 == 0 && i!=0)
		{
			lapTimer.stop();
			std::cout << "Finished iteration " << i << std::endl;
			auto lapTime = lapTimer.getRealTime();
			float glups = (cellCount * 100) / lapTime / 1000000000;
			std::cout << "GLUPS: " << glups << std::endl;
			lapTimer.reset();
			lapTimer.start();
		}
	}
	timer.stop();
	#endif
	auto totalTime = timer.getRealTime();
	std::cout << "This took " << totalTime << " s" << std::endl;
	float glups = (cellCount * iterationCount) / totalTime / 1000000000;
	std::cout << "Total average GLUPS: " << glups << std::endl;
	
	std::cout << "Saving result" << std::endl;
	
	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	// Write metadata first so Python knows the dimensions
	int dims[2] = {(int)cellCountY, (int)cellCountZ};
	fwrite(dims, sizeof(int), 2, fp);
	
	DistributionStructCPU FCPU;
	FCPU.shifter = IndexArrayType( 27, 0 );
	FCPU.fArray.setSizes( 27, cellCount );
	FCPU.shifter = F.shifter;
	FCPU.fArray = F.fArray;
	
	for (size_t j = 0; j < cellCountY; j++)
	{
		for (size_t k = 0; k < cellCountZ; k++)
		{
			size_t i = cellCountX / 2; 
			size_t cell = convertIndex(i, j, k);
			size_t shiftedIndex[27];
			for (size_t i = 0; i < 27; i++) 
			{
				const size_t shift = FCPU.shifter[i];
				shiftedIndex[i] = cell + shift;
				if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
			}
			float f[27];
			float rho, ux, uy, uz;
			for (size_t i = 0; i < 27; i++)	f[i] = FCPU.fArray.getElement(i, shiftedIndex[i]);
			getRhoUxUyUz(rho, ux, uy, uz, f);
			float uMag = sqrt(uy * uy + uz * uz);
			fwrite(&uMag, sizeof(float), 1, fp);
		}
	}
	fclose(fp);
	system("python3 visualizer.py");
	return EXIT_SUCCESS;
}
