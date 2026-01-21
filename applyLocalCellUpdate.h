#include "boundaryConditions/applyBounceback.h"

void applyLocalCellUpdate( 	MarkerStruct& Marker, DistributionStruct& F, CellCountStruct &cellCount )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getConstView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getConstView();
	auto givenRhoMarkerArrayView = Marker.givenRhoArray.getConstView();
	auto givenUxUyUzMarkerArrayView = Marker.givenUxUyUzArray.getConstView();
	
	auto shifterView = F.shifter.getConstView();
	auto fArrayView  = F.fArray.getView();

	auto cellLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount.n) { shiftedIndex[i] -= cellCount.n; }
		}
		bool fluidMarker = fluidMarkerArrayView[cell];
		bool bouncebackMarker = bouncebackMarkerArrayView[cell];
		bool givenRhoMarker = givenRhoMarkerArrayView[cell];
		bool givenUxUyUzMarker = givenUxUyUzMarkerArrayView[cell];
		
		float f[27];
		float rho, ux, uy, uz;
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		
		if ( bouncebackMarker == 1 )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker == 1 )
			{
				getRhoUxUyUz(rho, ux, uy, uz, f);
			}
			else
			{
				short outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal(cell, outerNormalX, outerNormalY, outerNormalZ, cellCount);
				rho = 1.f;
				ux = 0.f;
				uy = 0.f;
				uz = uzInlet;
				if ( givenRhoMarker == 1 && givenUxUyUzMarker == 0 )
				{
					restoreUxUyUz(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				else if ( givenRhoMarker == 0 && givenUxUyUzMarker == 1 )
				{
					restoreRho(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				else if ( givenRhoMarker == 0 && givenUxUyUzMarker == 0 )
				{
					restoreRhoUxUyUz(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
				}
				applyMBBC(outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f);
			}
			applyCollision(rho, ux, uy, uz, f);
		}
		for (size_t i = 0; i < 27; i++)	fArrayView(i, shiftedIndex[i]) = f[i];
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, cellCount.n, cellLambda );
}
