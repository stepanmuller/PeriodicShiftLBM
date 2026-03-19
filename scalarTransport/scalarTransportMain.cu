constexpr float channelLengthPhys = 2000.f;												// mm
constexpr float channelWidthPhys = 1000.f;												// mm
constexpr float resGlobal = 10.f; 														// mm
constexpr float uxInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 100000.f;
constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float nuPhys = 1.50e-5;														// m2/s air
constexpr float uxInletPhys = nuPhys * reynoldsNumber / (channelWidthPhys / 1000.f); 	// m/s
constexpr float rhoNominalPhys = 1.225f;												// kg/m3 air
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float scalarTransportDiffusivity = 0.1f; 										// dimensionless for now
constexpr float TInlet = 1.f;

const int cellCountX = static_cast<int>(std::ceil(channelLengthPhys / resGlobal));
const int cellCountY = static_cast<int>(std::ceil(channelWidthPhys / resGlobal));
const int cellCountZ = static_cast<int>(std::ceil(channelWidthPhys / resGlobal));

constexpr int iterationCount = 20000;
constexpr int iterationChunk = 1000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"

#include "../scalarTransportFunctions.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/applyMirror.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.bounceback = 1; // top bottom wall
	else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.bounceback = 1; // front back wall
	else if ( iCell <= Info.cellCountX / 4 && jCell <= Info.cellCountY / 2 ) Marker.bounceback = 1; // step wall
	else if ( iCell == 0 ) Marker.givenUxUyUz = 1; // velocity inlet
	else if ( iCell == Info.cellCountX-1 ) Marker.givenRho = 1; // pressure outlet (weakly compressible, rho is proportional to pressure)
	else Marker.fluid = 1;
}

__cuda_callable__ void getScalarTransportMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
													ScalarTransportMarkerStruct &ScalarTransportMarker, const InfoStruct& Info )
{
	if ( iCell == 0 ) ScalarTransportMarker.givenT = 1; // prescribed concentration inlet
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											const InfoStruct& Info )
{
    rho = 1.f;
	ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
}

__cuda_callable__ void getGivenT( 	const int& iCell, const int& jCell, const int& kCell, 
											float& T,
											const InfoStruct& Info )
{
    T = 0.f;
	if ( jCell <= Info.cellCountY / 2 + Info.cellCountY / 10 ) T = TInlet;
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	if ( Marker.bounceback ) ux = 0;
	else ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"

int main(int argc, char **argv)
{
	GridStruct Grid;
	Grid.Info.res = resGlobal;
	Grid.Info.dtPhys = dtPhysGlobal;
	Grid.Info.nu = (Grid.Info.dtPhys * nuPhys) / ((Grid.Info.res/1000) * (Grid.Info.res/1000));
	Grid.Info.cellCountX = cellCountX;
	Grid.Info.cellCountY = cellCountY;
	Grid.Info.cellCountZ = cellCountZ;
	Grid.Info.cellCount = Grid.Info.cellCountX * Grid.Info.cellCountY * Grid.Info.cellCountZ;
	Grid.fArray.setSizes( 27, Grid.Info.cellCount );
	Grid.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid );
	
	ScalarTransportStruct ScalarTransport;	
	ScalarTransport.TArray.setSizes( 27, Grid.Info.cellCount );
	ScalarTransport.TArray.setValue( 0.f );
	ScalarTransport.tauT = 3.f * scalarTransportDiffusivity + 0.5f;
	
	std::cout << "Cell count: " << Grid.Info.cellCount << std::endl;
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyStreaming(Grid);
		applyLocalCellUpdate(Grid, ScalarTransport);
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;

			const float updateCount = (float)Grid.Info.cellCount * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			const int kCut = Grid.Info.cellCountZ / 2;
			exportSectionCutPlotXY( Grid, kCut, iteration );
			system("python3 ../plotter/plotter.py");
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
