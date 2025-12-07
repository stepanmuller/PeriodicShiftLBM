#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>
#include <TNL/Timer.h>
#include <cmath>
#include <fstream> 

#include "config.h"

#include "types.h"

#include "applyInitialization.h"
#include "applyStreaming.h"
#include "convertIndex.h"
#include "convertNormal.h"

#include "updateCells/updateFluidCells.h"
#include "updateCells/updateBouncebackCells.h"
#include "updateCells/updateInletCells.h"
#include "updateCells/updateOutletCells.h"
#include "updateCells/updatePeriodicCells.h"

int main(int argc, char **argv)
{
	DistributionFunctionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	for (size_t i = 0; i < 27; i++) { F.fArray[i] = ArrayType(cellCount, 1.f); }
	
	RhoUGStruct rhoUG;
	rhoUG.rhoArray = ArrayType( cellCount, 1.f );
	rhoUG.uxArray = ArrayType( cellCount, 0.f );
	rhoUG.uyArray = ArrayType( cellCount, 0.f );
	rhoUG.uzArray = ArrayType( cellCount, 0.f );
	rhoUG.gxArray = ArrayType( cellCount, 0.f );
	rhoUG.gyArray = ArrayType( cellCount, 0.f );
	rhoUG.gzArray = ArrayType( cellCount, 0.f );
	
	CellGroupStruct fluidCells;
	CellGroupStruct bouncebackCells;
	CellGroupStruct inletCells;
	CellGroupStruct outletCells;
	CellGroupStruct periodicCells;
	
	std::cout << "Periodic Shift LBM" << std::endl;
	std::cout << "Initialization: Classifying cells" << std::endl;
	
	for (size_t k = 0; k < cellCountZ; k++)
	{
		for (size_t j = 0; j < cellCountY; j++)
		{
			for (size_t i = 0; i < cellCountX; i++)
			{
				size_t cell = convertIndex(i, j, k);
				if ( j>=350 && j<=450 && k>=250 && k<=350 ) bouncebackCells.groupSize++;
				else if ( i==0 || i==cellCountX-1 || j==0 || j==cellCountY-1 ) periodicCells.groupSize++;
				else if ( k==0 ) inletCells.groupSize++;
				else if ( k==cellCountZ-1 ) outletCells.groupSize++;
				else fluidCells.groupSize++;
			}
		}
	}
	
	fluidCells.indexArray.setSize(fluidCells.groupSize);
	bouncebackCells.indexArray.setSize(bouncebackCells.groupSize);
	inletCells.indexArray.setSize(inletCells.groupSize);
	inletCells.normalArray.setSize(inletCells.groupSize);
	outletCells.indexArray.setSize(outletCells.groupSize);
	outletCells.normalArray.setSize(outletCells.groupSize);
	periodicCells.indexArray.setSize(periodicCells.groupSize);
	periodicCells.sourceIndexArray.setSize(periodicCells.groupSize);
	periodicCells.normalArray.setSize(periodicCells.groupSize);
	
	fluidCells.CPUindexArray.setSize(fluidCells.groupSize);
	bouncebackCells.CPUindexArray.setSize(bouncebackCells.groupSize);
	inletCells.CPUindexArray.setSize(inletCells.groupSize);
	inletCells.CPUnormalArray.setSize(inletCells.groupSize);
	outletCells.CPUindexArray.setSize(outletCells.groupSize);
	outletCells.CPUnormalArray.setSize(outletCells.groupSize);
	periodicCells.CPUindexArray.setSize(periodicCells.groupSize);
	periodicCells.CPUsourceIndexArray.setSize(periodicCells.groupSize);
	periodicCells.CPUnormalArray.setSize(periodicCells.groupSize);
	
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
					bouncebackCells.CPUindexArray.setElement( bouncebackCells.temp, cell );
					bouncebackCells.temp++;
				}
				else if ( i==0 || i==cellCountX-1 || j==0 || j==cellCountY-1 ) 
				{
					periodicCells.CPUindexArray.setElement( periodicCells.temp, cell );
					if ( i==0 )
					{
						uint8_t normalCode;
						encodeNormal( normalCode, -1, 0, 0 );
						periodicCells.CPUnormalArray.setElement( periodicCells.temp, normalCode );
						size_t sourceI = cellCountX - 2;
						size_t sourceCell = convertIndex(sourceI, j, k);
						periodicCells.CPUsourceIndexArray.setElement( periodicCells.temp, sourceCell );
					}
					else if ( i==cellCountX-1 )
					{
						uint8_t normalCode;
						encodeNormal( normalCode, 1, 0, 0 );
						periodicCells.CPUnormalArray.setElement( periodicCells.temp, normalCode );
						size_t sourceI = 1;
						size_t sourceCell = convertIndex(sourceI, j, k);
						periodicCells.CPUsourceIndexArray.setElement( periodicCells.temp, sourceCell );
					}
					else if ( j==0 )
					{
						uint8_t normalCode;
						encodeNormal( normalCode, 0, -1, 0 );
						periodicCells.CPUnormalArray.setElement( periodicCells.temp, normalCode );
						size_t sourceJ = cellCountY - 2;
						size_t sourceCell = convertIndex(i, sourceJ, k);
						periodicCells.CPUsourceIndexArray.setElement( periodicCells.temp, sourceCell );
					}
					else if ( j==cellCountY-1 )
					{
						uint8_t normalCode;
						encodeNormal( normalCode, 0, 1, 0 );
						periodicCells.CPUnormalArray.setElement( periodicCells.temp, normalCode );
						size_t sourceJ = 1;
						size_t sourceCell = convertIndex(i, sourceJ, k);
						periodicCells.CPUsourceIndexArray.setElement( periodicCells.temp, sourceCell );
					}
					periodicCells.temp++;
				}
				else if ( k==0 ) 
				{
					inletCells.CPUindexArray.setElement( inletCells.temp, cell );
					uint8_t normalCode;
					encodeNormal( normalCode, 0, 0, -1 );
					inletCells.CPUnormalArray.setElement( inletCells.temp, normalCode );
					rhoUG.uzArray.setElement( cell, uzInlet );
					inletCells.temp++;
				}
				else if ( k==cellCountZ-1 ) 
				{
					outletCells.CPUindexArray.setElement( outletCells.temp, cell );
					outletCells.temp++;
				}
				else
				{
					fluidCells.CPUindexArray.setElement( fluidCells.temp, cell );
					uint8_t normalCode;
					encodeNormal( normalCode, 0, 0, 1 );
					outletCells.CPUnormalArray.setElement( outletCells.temp, normalCode );
					fluidCells.temp++;
				}
			}
		}
	}
	
	std::cout << "Passing cells to GPU" << std::endl;
	fluidCells.indexArray = fluidCells.CPUindexArray;
	bouncebackCells.indexArray = bouncebackCells.CPUindexArray;
	inletCells.indexArray = inletCells.CPUindexArray;
	inletCells.normalArray = inletCells.CPUnormalArray;
	outletCells.indexArray = outletCells.CPUindexArray;
	outletCells.normalArray = outletCells.CPUnormalArray;
	periodicCells.indexArray = periodicCells.CPUindexArray;
	periodicCells.sourceIndexArray = periodicCells.CPUsourceIndexArray;
	periodicCells.normalArray = periodicCells.CPUnormalArray;
	
	applyInitialization( F, rhoUG );
	
	#ifdef __CUDACC__
	std::cout << "Starting simulation" << std::endl;
	TNL::Timer timer;
	TNL::Timer timerLap;
	timer.start();
	timerLap.start();
	for (int i=0; i<iterationCount; i++)
	{
		applyStreaming( F );
		
		updateFluidCells( fluidCells, F, rhoUG );
		updateBouncebackCells( bouncebackCells, F, rhoUG );
		updateInletCells( inletCells, F, rhoUG );
		updateOutletCells( outletCells, F, rhoUG );
		updatePeriodicCells( periodicCells, F, rhoUG );
		
		if (i%100 == 0 && i!=0)
		{
			timerLap.stop();
			std::cout << "Finished iteration " << i << std::endl;
			auto lapTime = timerLap.getRealTime();
			float glups = (cellCount * 100) / lapTime / 1000000000;
			std::cout << "GLUPS: " << glups << std::endl;
			timerLap.reset();
			timerLap.start();
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
			float uy = rhoUG.uyArray.getElement(cell);
			float uz = rhoUG.uzArray.getElement(cell);
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
