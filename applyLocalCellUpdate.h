#include "config.h"
#include "boundaryConditions/applyBounceback.h"

void applyLocalCellUpdate( 	MarkerStruct& Marker, DistributionStruct& F )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getConstView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getConstView();
	auto inletMarkerArrayView = Marker.inletArray.getConstView();
	auto outletMarkerArrayView = Marker.outletArray.getConstView();
	
	auto shifterView = F.shifter.getConstView();
	auto fArrayView  = F.fArray.getView();

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
		bool bouncebackMarker = bouncebackMarkerArrayView[cell];
		bool inletMarker = inletMarkerArrayView[cell];
		bool outletMarker = outletMarkerArrayView[cell];
		
		float f[27];
		float rho, ux, uy, uz;
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		
		if ( fluidMarker == 1 )
		{
			getRhoUxUyUz(rho, ux, uy, uz, f);
			applyCollision(rho, ux, uy, uz, f);
		}
		else if ( bouncebackMarker == 1 )
		{
			applyBounceback(f);
		}
		else if ( inletMarker == 1 )
		{
			short outerNormalX, outerNormalY, outerNormalZ;
			getOuterNormal(cell, outerNormalX, outerNormalY, outerNormalZ);
			restoreRho(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			applyMBBC(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			applyCollision(rho, ux, uy, uz, f);
		}
		else if ( outletMarker == 1 )
		{
			short outerNormalX, outerNormalY, outerNormalZ;
			getOuterNormal(cell, outerNormalX, outerNormalY, outerNormalZ);
			restoreUxUyUz(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			applyMBBC(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			applyCollision(rho, ux, uy, uz, f);
		}
		for (size_t i = 0; i < 27; i++)	fArrayView(i, shiftedIndex[i]) = f[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount, cellLambda );
}
