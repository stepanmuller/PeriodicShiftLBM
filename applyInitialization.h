#include "config.h"

void applyInitialization( 	DistributionFunctionStruct& F, 
							ArrayType& rhoArray, 
							ArrayType& uxArray, ArrayType& uyArray, ArrayType& uzArray )
{
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
	
	auto rhoArrayView = rhoArray.getView();
	auto uxArrayView = uxArray.getView();
	auto uyArrayView = uyArray.getView();
	auto uzArrayView = uzArray.getView();

	auto initializeLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
		}
		
		const float rho = rhoArrayView[cell];
		const float ux = uxArrayView[cell];
		const float uy = uyArrayView[cell];
		const float uz = uzArrayView[cell];
		
		#include "inPlaceInclude/getFeq.hpp"
		
		f0ArrayView[shiftedIndex[0]]   = feq0;
		f1ArrayView[shiftedIndex[1]]   = feq1;
		f2ArrayView[shiftedIndex[2]]   = feq2;
		f3ArrayView[shiftedIndex[3]]   = feq3;
		f4ArrayView[shiftedIndex[4]]   = feq4;
		f5ArrayView[shiftedIndex[5]]   = feq5;
		f6ArrayView[shiftedIndex[6]]   = feq6;
		f7ArrayView[shiftedIndex[7]]   = feq7;
		f8ArrayView[shiftedIndex[8]]   = feq8;
		f9ArrayView[shiftedIndex[9]]   = feq9;
		f10ArrayView[shiftedIndex[10]] = feq10;
		f11ArrayView[shiftedIndex[11]] = feq11;
		f12ArrayView[shiftedIndex[12]] = feq12;
		f13ArrayView[shiftedIndex[13]] = feq13;
		f14ArrayView[shiftedIndex[14]] = feq14;
		f15ArrayView[shiftedIndex[15]] = feq15;
		f16ArrayView[shiftedIndex[16]] = feq16;
		f17ArrayView[shiftedIndex[17]] = feq17;
		f18ArrayView[shiftedIndex[18]] = feq18;
		f19ArrayView[shiftedIndex[19]] = feq19;
		f20ArrayView[shiftedIndex[20]] = feq20;
		f21ArrayView[shiftedIndex[21]] = feq21;
		f22ArrayView[shiftedIndex[22]] = feq22;
		f23ArrayView[shiftedIndex[23]] = feq23;
		f24ArrayView[shiftedIndex[24]] = feq24;
		f25ArrayView[shiftedIndex[25]] = feq25;
		f26ArrayView[shiftedIndex[26]] = feq26;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, initializeLambda );
}
