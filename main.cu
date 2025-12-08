#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Timer.h>
#include <cmath>
#include <fstream> 

using ArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using FlagArrayType = TNL::Containers::Array< int, TNL::Devices::Cuda, size_t >;
using CPUFlagArrayType = TNL::Containers::Array< int, TNL::Devices::Host, size_t >;

struct DistributionFunctionStruct { IndexArrayType shifter; ArrayType fArray[27]; };

#include "config.h"

#include "convertIndex.h"

#include "applyInitialization.h"
#include "applyStreaming.h"
#include "applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	DistributionFunctionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	for (size_t i = 0; i < 27; i++) { F.fArray[i] = ArrayType(cellCount, 1.f); }
	
	ArrayType rhoArray = ArrayType( cellCount, 1.f );
	ArrayType uxArray = ArrayType( cellCount, 0.f );
	ArrayType uyArray = ArrayType( cellCount, 0.f );
	ArrayType uzArray = ArrayType( cellCount, 0.f );
	ArrayType gxArray = ArrayType( cellCount, 0.f );
	ArrayType gyArray = ArrayType( cellCount, 0.f );
	ArrayType gzArray = ArrayType( cellCount, 0.f );
	
	FlagArrayType flagArray = FlagArrayType( cellCount, 0);	// 1 = fluid, 2 = bounceback, 1NNN = velocity inlet, 2NNN = pressure outlet
	CPUFlagArrayType CPUFlagArray = CPUFlagArrayType( cellCount, 0);
	
	std::cout << "Periodic Shift LBM" << std::endl;
	std::cout << "Initialization: Flagging cells" << std::endl;
	
	for (size_t k = 0; k < cellCountZ; k++)
	{
		for (size_t j = 0; j < cellCountY; j++)
		{
			for (size_t i = 0; i < cellCountX; i++)
			{
				size_t cell = convertIndex(i, j, k);
				
				if (cell%1000000 == 0){std::cout << "Classifying cell " << cell << " out of " << cellCount << std::endl;}
				
				if ( j>=350 && j<=450 && k>=250 && k<=350 ) 
				{
					CPUFlagArray.setElement( cell, 2 ); // 2 = bounceback
				}
				else if ( i==0 || i==cellCountX-1 || j==0 || j==cellCountY-1 ) 
				{
					CPUFlagArray.setElement( cell, 2 ); // 2 = bounceback
				}
				else if ( k==0 ) 
				{
					CPUFlagArray.setElement( cell, 1554 ); // 1NNN = velocity inlet, 554 -> outer normal (0, 0, -1)
				}
				else if ( k==cellCountZ-1 ) 
				{
					CPUFlagArray.setElement( cell, 2556 ); // 2NNN = pressure outlet, 556 -> outer normal (0, 0, 1)
				}
				else
				{
					CPUFlagArray.setElement( cell, 1 ); // 1 = fluid
				}
			}
		}
	}
	
	std::cout << "Passing cells to GPU" << std::endl;
	flagArray = CPUflagArray;
	
	applyInitialization( F, rho, ux, uy, uz );
	
	#ifdef __CUDACC__
	std::cout << "Starting simulation" << std::endl;
	TNL::Timer timer;
	TNL::Timer lapTimer;
	timer.start();
	lapTimer.start();
	for (int i=0; i<iterationCount; i++)
	{
		applyStreaming( F );
		applyLocalCellUpdate( flagArray, F, rhoArray, uxArray, uyArray, uzArray, gxArray, gyArray, gzArray );
		
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
