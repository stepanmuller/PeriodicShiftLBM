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
	
	RhoUGStruct RhoUG;
	RhoUG.rhoArray = ArrayType( cellCount, 1.f );
	RhoUG.uxArray = ArrayType( cellCount, 0.f );
	RhoUG.uyArray = ArrayType( cellCount, 0.f );
	RhoUG.uzArray = ArrayType( cellCount, 0.f );
	RhoUG.gxArray = ArrayType( cellCount, 0.f );
	RhoUG.gyArray = ArrayType( cellCount, 0.f );
	RhoUG.gzArray = ArrayType( cellCount, 0.f );
	
	CellGroupStruct fluidCells;
	CellGroupStruct bouncebackCells;
	CellGroupStruct inletCells;
	CellGroupStruct outletCells;
	CellGroupStruct periodicCells;
	
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
	std::cout << fluidCells.groupSize << std::endl;
	
	/*
	applyInitialization( F, rhoUG );
	
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
	*/
	return EXIT_SUCCESS;
}
