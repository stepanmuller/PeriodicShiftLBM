// input physics
constexpr float res = 0.35f; 											// mm
constexpr float uzInlet = 0.1f; 										// also works as nominal LBM Mach number

constexpr float nuPhys = 1e-6;											// m2/s water
constexpr float rhoNominalPhys = 1000.0f;								// kg/m3 water
constexpr float SmagorinskyConstantGlobal = 0.1f; 						// set to zero to turn off LES

// calculated from input
constexpr float uzInletPhys = 17.f; 								// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * (res/1000); 		// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (res/1000) / dtPhys; 	// m/s
constexpr float nu = (dtPhys * nuPhys) / ((res/1000) * (res/1000));	// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau

constexpr int iterationCount = 500001;

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

#include "../STLFunctions.h"

std::string STLPath = "M35IntakeSTL.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& mirrorMarker, bool& periodicMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									InfoStruct& Info )
{
    fluidMarker = 0;
	mirrorMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	if ( bouncebackMarker )	{ } // do nothing but skip all else
	else if ( kCell == 0 || jCell == 0 ) givenUxUyUzMarker = 1;
	else if ( iCell == 0 || iCell == Info.cellCountX-1 ) givenUxUyUzMarker = 1;
	else if ( kCell == Info.cellCountZ-1 ) givenRhoMarker = 1;
	else fluidMarker = 1;
	
	periodicMarker = 0;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
    rho = 1.f;
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
}

__cuda_callable__ float getSmagorinskyConstant( const int& iCell, const int& jCell, const int& kCell, const InfoStruct &Info )
{
	const float xToEnd = (Info.cellCountX-1 - iCell) * Info.res;
	if ( xToEnd > 30.f ) return SmagorinskyConstantGlobal;
	else
	{
		const float interpoler = xToEnd / 30.f;
		return SmagorinskyConstantGlobal * interpoler + 1.f * (1.f - interpoler); // set Smagorinsky high to dampen vortices before reaching the outlet
	}
	
}

#include "../applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	STLStructCPU STLCPU;
	readSTL( STLCPU, STLPath );
	
	InfoStruct Info;
	Info.res = res;
	Info.rhoNominalPhys = rhoNominalPhys;
	Info.soundspeedPhys = soundspeedPhys;
	Info.dtPhys = dtPhys;
	
	std::cout << "dtPhys: " << dtPhys << " s" << std::endl;
	
	std::cout << "Sizing domain around the STL" << std::endl;
	
	Info.cellCountX = static_cast<int>( std::ceil(( STLCPU.xmax - STLCPU.xmin - 1e-9 ) / Info.res ));
	Info.cellCountY = static_cast<int>( std::ceil(( STLCPU.ymax - STLCPU.ymin - 1e-9 ) / Info.res ));
	Info.cellCountZ = static_cast<int>( std::ceil(( STLCPU.zmax - STLCPU.zmin - 1e-9 ) / Info.res ));
	
	STLCPU.ox = - STLCPU.xmin - ( 0.5f * ( ( STLCPU.xmax - STLCPU.xmin ) - Info.res * ( Info.cellCountX-1 ) ) );
	STLCPU.oy = - STLCPU.ymin - ( 0.5f * ( ( STLCPU.ymax - STLCPU.ymin ) - Info.res * ( Info.cellCountY-1 ) ) );
	STLCPU.oz = - STLCPU.zmin - ( 0.5f * ( ( STLCPU.zmax - STLCPU.zmin ) - Info.res * ( Info.cellCountZ-1 ) ) );
	
	std::cout << "	STL.ox: " << STLCPU.ox << std::endl;
	std::cout << "	STL.oy: " << STLCPU.oy << std::endl;
	std::cout << "	STL.oz: " << STLCPU.oz << std::endl;
	
	Info.cellCountY = Info.cellCountY + 1; // adding one more "wall" layer on top
	
	Info.cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	std::cout << "	cellCountX: " << Info.cellCountX << std::endl;
	std::cout << "	cellCountY: " << Info.cellCountY << std::endl;
	std::cout << "	cellCountZ: " << Info.cellCountZ << std::endl;
	std::cout << "	cellCount: " << Info.cellCount << std::endl;
	
	STLStruct STL( STLCPU );
	
	checkSTLEdges( STL );
	
	BoolArrayType bouncebackArray = BoolArrayType( Info.cellCount, 0 );
	
	const bool insideMarkerValue = 0;
	applyMarkersInsideSTL( bouncebackArray, STL, insideMarkerValue, Info );
	
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
		applyLocalCellUpdate( F, bouncebackArray, Info );
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			const int cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
			float glups = ((float)cellCount * (float)iterationChunk) / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlotZY( F, bouncebackArray, Info, iCut, plotNumber );
			plotNumber++;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
