constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES
constexpr float nu = 1e-6;												// LBM nu
constexpr float tau = 3.f * nu + 0.5f;									// LBM tau

constexpr int cellCountX = 50;
constexpr int cellCountY = 1000;
constexpr int cellCountZ = 2000;

constexpr int boxStartY = (int)(cellCountY * 0.4f);
constexpr int boxStartZ = (int)(cellCountY * 0.3f);
constexpr int boxEndY = (int)(cellCountY * 0.6f);
constexpr int boxEndZ = (int)(cellCountY * 0.5f);

constexpr int iterationCount = 10000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"
#include "../fillDefaultEquilibrium.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

#include "../exportSectionCutPlot.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									InfoStruct& Info )
{
    fluidMarker = 0;
	bouncebackMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	if ( jCell >= boxStartY && jCell < boxEndY && kCell >= boxStartZ && kCell < boxEndZ ) bouncebackMarker = 1;
	else if ( iCell == 0 || iCell == Info.cellCountX-1 || jCell == 0 || jCell == Info.cellCountY-1 ) bouncebackMarker = 1;
	else if ( kCell == 0 ) givenUxUyUzMarker = 1;
	else if ( kCell == Info.cellCountZ-1 ) givenRhoMarker = 1;
	else fluidMarker = 1;
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
	
	const int iCut = Info.cellCountX / 2;
	int plotNumber = 0;
	
	const int iterationChunk = 100;
	
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
			
			exportSectionCutPlotZY( F, Info, iCut, plotNumber );
			plotNumber++;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
