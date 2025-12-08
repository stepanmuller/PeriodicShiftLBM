#include "config.h"

void applyLocalCellUpdate( 	FlagArrayType& flagArray,
							MarkerArrayType& fluidMarkerArray,
							MarkerArrayType& bouncebackMarkerArray,
							DistributionFunctionStruct& F, 
							ArrayType& rhoArray, 
							ArrayType& uxArray, ArrayType& uyArray, ArrayType& uzArray,
							ArrayType& gxArray, ArrayType& gyArray, ArrayType& gzArray )
{
	auto flagArrayView = flagArray.getConstView();
	
	auto fluidMarkerArrayView = fluidMarkerArray.getConstView();
	auto bouncebackMarkerArrayView = bouncebackMarkerArray.getConstView();
	
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

	auto cellLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) { shiftedIndex[i] -= cellCount; }
		}
		bool fluidMarker = fluidMarkerArrayView[cell];
		if ( fluidMarker == 1 ) // fluid
		{
			#include "inPlaceInclude/readF.hpp"
			#include "inPlaceInclude/getRhoUxUyUz.hpp"
			#include "inPlaceInclude/applyCollision.hpp"
			#include "inPlaceInclude/writeF.hpp"
			#include "inPlaceInclude/writeRho.hpp"
			#include "inPlaceInclude/writeUxUyUz.hpp"
			return;
		}
		bool bouncebackMarker = bouncebackMarkerArrayView[cell];
		if ( bouncebackMarker == 1 ) // bounceback
		{
			#include "inPlaceInclude/readF.hpp"
			#include "inPlaceInclude/applyBounceback.hpp"
			#include "inPlaceInclude/writeF.hpp"
			return;
		}
		short flag = flagArrayView[cell];
		if ( flag == 0 ) return; // ignore cell
		else if ( flag > 1000 && flag < 2000 ) // velocity inlet
		{
			#include "inPlaceInclude/applyVelocityInlet.hpp"
			#include "inPlaceInclude/applyCollision.hpp"
			#include "inPlaceInclude/writeF.hpp"
			#include "inPlaceInclude/writeRho.hpp"
			return;
		}
		else if ( flag > 2000 && flag < 3000 ) // pressure outlet
		{
			#include "inPlaceInclude/applyPressureOutlet.hpp"
			#include "inPlaceInclude/applyCollision.hpp"
			#include "inPlaceInclude/writeF.hpp"
			#include "inPlaceInclude/writeUxUyUz.hpp"
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, cellLambda );
}
