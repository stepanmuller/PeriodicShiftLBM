constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES
constexpr float nu = 1e-6;												// LBM nu
constexpr float tau = 3.f * nu + 0.5f;									// LBM tau

constexpr int cellCountX = 5;
constexpr int cellCountY = 100;
constexpr int cellCountZ = 200;

constexpr int boxStartY = 400;
constexpr int boxStartZ = 400;
constexpr int boxEndY = 600;
constexpr int boxEndZ = 600;

constexpr int iterationCount = 10000;

#include "../includesTypes.h"
#include "exportSectionCutPlot.h"

void applyLocalCellUpdate( FloatArray4DType fArray, InfoStruct &Info )
{
	auto fArrayView  = fArray.getView();

	auto cellLambda = [=] __cuda_callable__ ( const TripleIndexType& tripleIndex ) mutable
	{
		const int iCell = tripleIndex.x();
		const int jCell = tripleIndex.y();
		const int kCell = tripleIndex.z();
		int iStreamed[27], jStreamed[27], kStreamed[27];
		getStreamedIndexes( iCell, jCell, kCell, iStreamed, jStreamed, kStreamed, Info );
		
		bool fluidMarker = 0;
		bool bouncebackMarker = 0;
		bool givenRhoMarker = 0;
		bool givenUxUyUzMarker = 0;
		
		if ( jCell >= boxStartY && jCell < boxEndY && kCell >= boxStartZ && kCell < boxEndZ ) bouncebackMarker = 1;
		else if ( kCell == 0 ) givenUxUyUzMarker = 1;
		else if ( kCell == Info.cellCountZ-1 ) givenRhoMarker = 1;
		else fluidMarker = 1;
		
		float f[27];
		float rho, ux, uy, uz;
		for ( int direction = 0; direction < 27; direction++ ) f[direction] = fArrayView( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] );
		
		if ( bouncebackMarker )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker )
			{
				getRhoUxUyUz( rho, ux, uy, uz, f );
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info );
				rho = 1.f;
				ux = 0.f;
				uy = 0.f;
				uz = uzInlet;
				if ( givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && givenUxUyUzMarker )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			applyCollision( rho, ux, uy, uz, f );
		}
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] ) = f[direction];
	};
	TripleIndexType start{ 0, 0, 0 };
	TripleIndexType end{ Info.cellCountX, Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( start, end, cellLambda );
}

int main(int argc, char **argv)
{
	InfoStruct Info;
	Info.cellCountX = cellCountX;
	Info.cellCountY = cellCountY;
	Info.cellCountZ = cellCountZ;
	Info.iterationsFinished = 0;
	
	FloatArray4DType fArray;
	fArray.setSizes( 27, Info.cellCountX, Info.cellCountY, Info.cellCountZ );	
	fillDefaultEquilibrium( fArray, Info);
	
	std::cout << "Starting simulation" << std::endl;
	
	const int iCut = Info.cellCountX / 2;
	int plotNumber = 0;
	
	TNL::Timer lapTimer;
	lapTimer.start();
	for (int iteration=0; iteration<iterationCount; iteration++)
	{
		applyLocalCellUpdate( fArray, Info );
		Info.iterationsFinished++;
		
		if (iteration%1000 == 0 && iteration!=0)
		{
			lapTimer.stop();
			std::cout << "Finished iteration " << iteration << std::endl;
			auto lapTime = lapTimer.getRealTime();
			float glups = (Info.cellCountX * Info.cellCountY * Info.cellCountZ * 1000) / lapTime / 1000000000;
			std::cout << "GLUPS: " << glups << std::endl;
			lapTimer.reset();
			
			exportSectionCutPlot( fArray, Info, iCut, plotNumber );
			plotNumber++;
			
			//check
			
			const int iCell = 0;
			const int jCell = 0;
			const int kCell = 0;
			int iStreamed[27], jStreamed[27], kStreamed[27];
			getStreamedIndexes( iCell, jCell, kCell, iStreamed, jStreamed, kStreamed, Info );
			float f[27];
			float rho, ux, uy, uz;
			for ( int direction = 0; direction < 27; direction++ ) f[direction] = fArray.getElement( direction, iStreamed[direction], jStreamed[direction], kStreamed[direction] );
			getRhoUxUyUz( rho, ux, uy, uz, f );
			
			std::cout << "uz" << uz << std::endl;
			
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
