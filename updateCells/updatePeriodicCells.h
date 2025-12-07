#include "../config.h"
#include "../convertNormal.h"
#include "../types.h"

void updatePeriodicCells( 	CellGroupStruct& cells,
							DistributionFunctionStruct& F, 
							RhoUGStruct& rhoUG )
{
	size_t groupSize = cells.groupSize;
	auto indexArrayView = cells.indexArray.getConstView();
	auto sourceIndexArrayView = cells.sourceIndexArray.getConstView();
	auto normalArrayView = cells.normalArray.getConstView();
	
	auto shifterView = F.shifter.getConstView();
	
	auto f0ArrayView  = F.fArray[0].getView();
	auto f1ArrayView  = F.fArray[1].getView();
	auto f2ArrayView  = F.fArray[2].getView();
	auto f3ArrayView  = F.fArray[3].getView();
	auto f4ArrayView  = F.fArray[4].getView();
	auto f5ArrayView  = F.fArray[5].getView();
	auto f6ArrayView  = F.fArray[6].getView();
	auto f7ArrayView  = F.fArray[7].getView();
	auto f8ArrayView  = F.fArray[8].getView();
	auto f9ArrayView  = F.fArray[9].getView();
	auto f10ArrayView = F.fArray[10].getView();
	auto f11ArrayView = F.fArray[11].getView();
	auto f12ArrayView = F.fArray[12].getView();
	auto f13ArrayView = F.fArray[13].getView();
	auto f14ArrayView = F.fArray[14].getView();
	auto f15ArrayView = F.fArray[15].getView();
	auto f16ArrayView = F.fArray[16].getView();
	auto f17ArrayView = F.fArray[17].getView();
	auto f18ArrayView = F.fArray[18].getView();
	auto f19ArrayView = F.fArray[19].getView();
	auto f20ArrayView = F.fArray[20].getView();
	auto f21ArrayView = F.fArray[21].getView();
	auto f22ArrayView = F.fArray[22].getView();
	auto f23ArrayView = F.fArray[23].getView();
	auto f24ArrayView = F.fArray[24].getView();
	auto f25ArrayView = F.fArray[25].getView();
	auto f26ArrayView = F.fArray[26].getView();

	auto cellLambda = [=] __cuda_callable__ (size_t index) mutable
	{
		size_t cell = indexArrayView[index];
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
		}
		
		size_t sourceCell = sourceIndexArrayView[index];
		size_t shiftedSourceIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedSourceIndex[i] = cell + shift;
			if (shiftedSourceIndex[i] >= cellCount) { shiftedSourceIndex[i] -= cellCount; }
		}
		
		uint8_t normalCode = normalArrayView[index];
		int outerNormalX, outerNormalY, outerNormalZ;
		decodeNormal( normalCode, outerNormalX, outerNormalY, outerNormalZ );
		
		#include "../includeInPlace/applyPeriodicBC.hpp"
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, groupSize, cellLambda );
}
