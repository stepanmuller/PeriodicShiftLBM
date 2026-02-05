constexpr float rhoInlet = 1.01f; 										
constexpr float SmagorinskyConstantGlobal = 0.0f; 						// set to zero to turn off LES
constexpr float nu = 1e-3;												// LBM nu
constexpr float tau = 3.f * nu + 0.5f;									// LBM tau

constexpr int cellCountX = 256;
constexpr int cellCountY = 256;
constexpr int cellCountZ = 256;

constexpr int iterationCount = 10000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"
#include "../fillDefaultEquilibrium.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/applyMirror.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

#include "../exportSectionCutPlot.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& mirrorMarker, bool& periodicMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									InfoStruct& Info )
{
    fluidMarker = 0;
	bouncebackMarker = 0;
	mirrorMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	if ( jCell == 0 || jCell == Info.cellCountY-1 ) bouncebackMarker = 1;
	else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) bouncebackMarker = 1;
	else if ( iCell == 0 ) givenRhoMarker = 1;
	else if ( iCell == Info.cellCountX-1 ) givenRhoMarker = 1;
	else fluidMarker = 1;
	
	periodicMarker = 0;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	if ( iCell == 0 ) rho = rhoInlet;
	else if ( iCell == Info.cellCountX-1 ) rho = 1.f;
}

__cuda_callable__ float getSmagorinskyConstant( const int& iCell, const int& jCell, const int& kCell )
{
	return SmagorinskyConstantGlobal;
}

#include "../applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	InfoStruct Info;
	Info.cellCountX = cellCountX;
	Info.cellCountY = cellCountY;
	Info.cellCountZ = cellCountZ;
	Info.cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	
	FStruct F;
	FloatArray2DType fArray;
	F.fArray.setSizes( 27, Info.cellCountX * Info.cellCountY * Info.cellCountZ );	
	F.shifter = IntArrayType( 27, 0 );
	
	fillDefaultEquilibrium( F, Info);
	
	std::cout << "Starting simulation" << std::endl;
	
	const int kCut = Info.cellCountZ / 2;
	int plotNumber = 0;
	
	const int iterationChunk = 1000;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<iterationCount; iteration++)
	{
		applyStreaming( F, Info );
		applyLocalCellUpdate( F, Info );
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			const int cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
			float glups = ((float)cellCount * (float)iterationChunk) / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlotXY( F, Info, kCut, plotNumber );
			plotNumber++;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
