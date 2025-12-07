#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Timer.h>
#include <cmath>
#include <fstream> 

#include "config.h"

using ArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using IndexArrayType = TNL::Containers::Array< size_t, TNL::Devices::Cuda, size_t >;
using SolidmaskType = TNL::Containers::Array< bool, TNL::Devices::Cuda, size_t >;

struct DistributionFunctionStruct { IndexArrayType shifter; ArrayType fArray[27]; };

#include "applyInitialization.h"
#include "applyCollision.h"
#include "applyStreaming.h"
#include "applyBCInlet.h"
#include "applyBCOutlet.h"
#include "convertIndex.h"

int main(int argc, char **argv)
{
	DistributionFunctionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	for (size_t i = 0; i < 27; i++) { F.fArray[i] = ArrayType(cellCount, 1.f); }
	
	ArrayType rhoArray( cellCount, 1.f );
	ArrayType uxArray( cellCount, 0.f );
	ArrayType uyArray( cellCount, 0.f );
	ArrayType uzArray( cellCount, 0.f );
	
	ArrayType gxArray( cellCount, 0.f );
	ArrayType gyArray( cellCount, 0.f );
	ArrayType gzArray( cellCount, 0.f );
	
	SolidmaskType solidmask( cellCount, 0 );
	for (size_t k = 0; k < cellCountZ; k++)
	{
		for (size_t i = 0; i < cellCountX; i++)
		{
			size_t j = 0;
			size_t cell = convertIndex(i, j, k);
			solidmask.setElement(cell, 1);
			j = cellCountY-1;
			cell = convertIndex(i, j, k);
			solidmask.setElement(cell, 1);
		}
		for (size_t j = 0; j < cellCountY; j++)
		{
			size_t i = 0;
			size_t cell = convertIndex(i, j, k);
			solidmask.setElement(cell, 1);
			i = cellCountX-1;
			cell = convertIndex(i, j, k);
			solidmask.setElement(cell, 1);
		}
	}
	for (size_t i = 0; i < cellCountX; i++)
	{
		for (size_t j = 350; j < 451; j++)
		{
			for (size_t k = 250; k < 351; k++)
			{
				size_t cell = convertIndex(i, j, k);
				solidmask.setElement(cell, 1);
			}
		}
	}
	
	IndexArrayType inletIndexArray( cellCountX * cellCountY );
	size_t counter = 0;
	for (size_t i = 1; i < cellCountX-1; i++)
	{
		for (size_t j = 1; j < cellCountY-1; j++)
		{
			size_t k = 0; //first cell layer in Z direction
			size_t cell = convertIndex(i, j, k);
			inletIndexArray.setElement(counter, cell);
			counter++;
		}
	}
	
	IndexArrayType outletIndexArray( (cellCountX-2) * (cellCountY-2) );
	counter = 0;
	for (size_t i = 1; i < cellCountX-1; i++)
	{
		for (size_t j = 1; j < cellCountY-1; j++)
		{
			size_t k = cellCountZ - 1; //last cell layer in Z direction
			size_t cell = convertIndex(i, j, k);
			outletIndexArray.setElement(counter, cell);
			counter++;
		}
	}
	
	applyInitialization( F, rhoArray, uxArray, uyArray, uzArray );
	
	#ifdef __CUDACC__
	std::cout << "Tnl periodic shift streaming + cummulant collision" << std::endl;
	TNL::Timer timer;
	timer.start();
	for (int i=0; i<iterationCount; i++)
	{
		applyStreaming( F );
		applyBCInlet( F, rhoArray, uxArray, uyArray, uzArray, inletIndexArray );
		applyBCOutlet( F, rhoArray, uxArray, uyArray, uzArray, outletIndexArray );
		applyCollision(	F, rhoArray, uxArray, uyArray, uzArray, gxArray, gyArray, gzArray, solidmask );
		if (i%100 == 0){std::cout << "Finished iteration " << i << std::endl;}
	}
	timer.stop();
	#endif
	auto totalTime = timer.getRealTime();
	std::cout << "This took " << totalTime << " s" << std::endl;
	float glups = (cellCount * iterationCount) / totalTime / 1000000000;
	std::cout << "GLUPS: " << glups << std::endl;
	
	std::cout << "Saving result" << std::endl;
	
	std::ofstream out("result.csv");
	for (size_t j = 0; j < cellCountY; j++)
	{
		for (size_t k = 0; k < cellCountZ; k++)
		{
			size_t i = cellCountX / 2; 
			size_t cell = convertIndex(i, j, k);
			float uy = uyArray.getElement(cell);
			float uz = uzArray.getElement(cell);
			float uMag = sqrt(uy * uy + uz * uz);
			out << uMag;
			if (k + 1 < cellCountZ)
				out << ",";
		}
		out << "\n";   // new row
	}
	out.close();
	
	return EXIT_SUCCESS;
}
