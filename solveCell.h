#include "config.h"

void applyCollision( 	DistributionFunctionStruct& F, 
						ArrayType& rhoArray, 
						ArrayType& uxArray, ArrayType& uyArray, ArrayType& uzArray, 
						ArrayType& gxArray, ArrayType& gyArray, ArrayType& gzArray, 
						SolidmaskType& solidmask )
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
	auto gxArrayView = gxArray.getView();
	auto gyArrayView = gyArray.getView();
	auto gzArrayView = gzArray.getView();
	
	auto solidmaskView = solidmask.getConstView();

	auto collideLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
		}
		
		float f0  = f0ArrayView[shiftedIndex[0]];
		float f1  = f1ArrayView[shiftedIndex[1]];
		float f2  = f2ArrayView[shiftedIndex[2]];
		float f3  = f3ArrayView[shiftedIndex[3]];
		float f4  = f4ArrayView[shiftedIndex[4]];
		float f5  = f5ArrayView[shiftedIndex[5]];
		float f6  = f6ArrayView[shiftedIndex[6]];
		float f7  = f7ArrayView[shiftedIndex[7]];
		float f8  = f8ArrayView[shiftedIndex[8]];
		float f9  = f9ArrayView[shiftedIndex[9]];
		float f10 = f10ArrayView[shiftedIndex[10]];
		float f11 = f11ArrayView[shiftedIndex[11]];
		float f12 = f12ArrayView[shiftedIndex[12]];
		float f13 = f13ArrayView[shiftedIndex[13]];
		float f14 = f14ArrayView[shiftedIndex[14]];
		float f15 = f15ArrayView[shiftedIndex[15]];
		float f16 = f16ArrayView[shiftedIndex[16]];
		float f17 = f17ArrayView[shiftedIndex[17]];
		float f18 = f18ArrayView[shiftedIndex[18]];
		float f19 = f19ArrayView[shiftedIndex[19]];
		float f20 = f20ArrayView[shiftedIndex[20]];
		float f21 = f21ArrayView[shiftedIndex[21]];
		float f22 = f22ArrayView[shiftedIndex[22]];
		float f23 = f23ArrayView[shiftedIndex[23]];
		float f24 = f24ArrayView[shiftedIndex[24]];
		float f25 = f25ArrayView[shiftedIndex[25]];
		float f26 = f26ArrayView[shiftedIndex[26]];
		float rho, ux, uy, uz;
		float gx, gy, gz;
		bool solidCell = solidmaskView[cell];
		
		if (solidCell == 0)
		{ 
		gx  = gxArrayView[cell];
		gy  = gyArrayView[cell];
		gz  = gzArrayView[cell];
		#include "applyCollisionEquations.hpp"
		}
		else
		{
		#include "applyBounceback.hpp"
		}
		
		f0ArrayView[shiftedIndex[0]]   = f0;
		f1ArrayView[shiftedIndex[1]]   = f1;
		f2ArrayView[shiftedIndex[2]]   = f2;
		f3ArrayView[shiftedIndex[3]]   = f3;
		f4ArrayView[shiftedIndex[4]]   = f4;
		f5ArrayView[shiftedIndex[5]]   = f5;
		f6ArrayView[shiftedIndex[6]]   = f6;
		f7ArrayView[shiftedIndex[7]]   = f7;
		f8ArrayView[shiftedIndex[8]]   = f8;
		f9ArrayView[shiftedIndex[9]]   = f9;
		f10ArrayView[shiftedIndex[10]] = f10;
		f11ArrayView[shiftedIndex[11]] = f11;
		f12ArrayView[shiftedIndex[12]] = f12;
		f13ArrayView[shiftedIndex[13]] = f13;
		f14ArrayView[shiftedIndex[14]] = f14;
		f15ArrayView[shiftedIndex[15]] = f15;
		f16ArrayView[shiftedIndex[16]] = f16;
		f17ArrayView[shiftedIndex[17]] = f17;
		f18ArrayView[shiftedIndex[18]] = f18;
		f19ArrayView[shiftedIndex[19]] = f19;
		f20ArrayView[shiftedIndex[20]] = f20;
		f21ArrayView[shiftedIndex[21]] = f21;
		f22ArrayView[shiftedIndex[22]] = f22;
		f23ArrayView[shiftedIndex[23]] = f23;
		f24ArrayView[shiftedIndex[24]] = f24;
		f25ArrayView[shiftedIndex[25]] = f25;
		f26ArrayView[shiftedIndex[26]] = f26;
		
		if (solidCell == 0)
		{
		rhoArrayView[cell] = rho;
		uxArrayView[cell] = ux;
		uyArrayView[cell] = uy;
		uzArrayView[cell] = uz;
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, collideLambda );
}
