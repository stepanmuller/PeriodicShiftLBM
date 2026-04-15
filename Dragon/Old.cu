constexpr float res = 1.f; 												// mm
constexpr float uxInlet = 0.05f; 										// also works as nominal LBM Mach number
constexpr float angleOfAttack = -10.f;									// deg

constexpr float nuPhys = 1.5e-5;										// m2/s air
constexpr float rhoNominalPhys = 1.225f;								// kg/m3 water
constexpr float SmagorinskyConstantGlobal = 0.0f; 						// set to zero to turn off LES

constexpr float uxInletPhys = 90.f; 									// m/s
constexpr float dtPhys = (uxInlet / uxInletPhys) * (res/1000); 			// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (res/1000) / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / ((res/1000) * (res/1000));		// LBM nu
constexpr float tau = 3.f * nu + 0.5f;									// LBM tau

constexpr float domainSizeX = 3000.f;									// mm
constexpr float domainSizeY = 3000.f;									// mm
constexpr float domainSizeZ = 10.f;										// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizeX / res));
const int cellCountY = static_cast<int>(std::ceil(domainSizeY / res));
const int cellCountZ = static_cast<int>(std::ceil(domainSizeZ / res));

constexpr int iterationCount = 10000001;

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

std::string STLPath = "NACA0012.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& mirrorMarker, bool& periodicMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									InfoStruct& Info )
{
    fluidMarker = 0;
	mirrorMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	if ( iCell == 0 || jCell == 0 || jCell == Info.cellCountY-1 ) givenUxUyUzMarker = 1;
	else if ( iCell == Info.cellCountX-1 ) givenRhoMarker = 1;
	else if ( !bouncebackMarker ) fluidMarker = 1;
	
	periodicMarker = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
    rho = 1.f;
	ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
}

__cuda_callable__ float getSmagorinskyConstant( const int& iCell, const int& jCell, const int& kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

#include "../applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	InfoStruct Info;
	Info.res = res;
	Info.cellCountX = cellCountX;
	Info.cellCountY = cellCountY;
	Info.cellCountZ = cellCountZ;
	Info.cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	std::cout << "Cell count: " << Info.cellCount << std::endl;
	Info.rhoNominalPhys = rhoNominalPhys;
	Info.soundspeedPhys = soundspeedPhys;
	Info.dtPhys = dtPhys;
	
	STLStructCPU STLCPU;
	readSTL( STLCPU, STLPath );
	STLCPU.ox = ((Info.cellCountX - 1) * Info.res) * 0.5f;
	STLCPU.oy = ((Info.cellCountY - 1) * Info.res) * 0.5f;
	STLCPU.oz = ((Info.cellCountZ - 1) * Info.res) * 0.5f;
	
	STLStruct STL( STLCPU );
	
	checkSTLEdges( STL );
	
	float radians = angleOfAttack * (M_PI / 180.0f);
	rotateSTLAlongZ( STL, radians );
	
	BoolArrayType bouncebackArray = BoolArrayType( Info.cellCount, 0 );
	
	const bool insideMarkerValue = 1;
	applyMarkersInsideSTL( bouncebackArray, STL, insideMarkerValue, Info );
	
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
		applyLocalCellUpdate( F, bouncebackArray, Info );
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			const int cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
			float glups = ((float)cellCount * (float)iterationChunk) / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlotXY( F, bouncebackArray, Info, kCut, plotNumber+100000 );
			exportSectionCutPlotXY( F, Info, kCut, plotNumber );
			plotNumber++;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
