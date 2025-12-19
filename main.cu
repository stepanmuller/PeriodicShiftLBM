#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Timer.h>
#include <cmath>
#include <fstream> 
#include <cstdlib>

using ArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using IndexArrayType = TNL::Containers::Array< size_t, TNL::Devices::Cuda, size_t >;

using MarkerArrayType = TNL::Containers::Array< bool, TNL::Devices::Cuda, size_t >;
using MarkerArrayTypeCPU = TNL::Containers::Array< bool, TNL::Devices::Host, size_t >;

struct DistributionFunctionStruct { IndexArrayType shifter; ArrayType fArray[27]; };

struct MarkerStruct { MarkerArrayType fluidArray; MarkerArrayType bouncebackArray; MarkerArrayType inletArray; MarkerArrayType outletArray;  };

#include "config.h"

#include "convertIndex.h"

#include "applyMarkers.h"
#include "applyInitialization.h"
#include "applyStreaming.h"
#include "applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	DistributionFunctionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	for (size_t i = 0; i < 27; i++) { F.fArray[i] = ArrayType(cellCount, 1.f); }
	
	MarkerStruct Marker;
	Marker.fluidArray = MarkerArrayType( cellCount, 0);
	Marker.bouncebackArray = MarkerArrayType( cellCount, 0);
	Marker.inletArray = MarkerArrayType( cellCount, 0);
	Marker.outletArray = MarkerArrayType( cellCount, 0);
	
	std::cout << "Periodic Shift LBM" << std::endl;
	
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
	
	std::ofstream out("result.csv");
	for (size_t j = 0; j < cellCountY; j++)
	{
		for (size_t k = 0; k < cellCountZ; k++)
		{
			size_t i = cellCountX / 2; 
			size_t cell = convertIndex(i, j, k);
			float uy = F.fArray[0].getElement(cell);
			float uz = F.fArray[0].getElement(cell);
			float uMag = sqrt(uy * uy + uz * uz);
			float rho = F.fArray[0].getElement(cell);
			out << uMag;
			if (k + 1 < cellCountZ)
				out << ",";
		}
		out << "\n";   // new row
	}
	out.close();
	int pythonResult = system("python3 visualizer.py");
	return EXIT_SUCCESS;
}
