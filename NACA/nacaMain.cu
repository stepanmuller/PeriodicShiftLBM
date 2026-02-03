constexpr float res = 1.2f; 											// mm
constexpr float uxInlet = 0.05f; 										// also works as nominal LBM Mach number

constexpr float nuPhys = 1.5e-5;										// m2/s air
constexpr float rhoNominalPhys = 1.225f;								// kg/m3 water
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES

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

constexpr int iterationCount = 50000;

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

std::string STLPath = "NACA0012Repaired.STL";

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

#include "../applyLocalCellUpdate.h"

int main(int argc, char **argv)
{
	InfoStruct Info;
	Info.res = res;
	Info.cellCountX = cellCountX;
	Info.cellCountY = cellCountY;
	Info.cellCountZ = cellCountZ;
	Info.cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	
	STLStructCPU STLCPU;
	readSTL( STLCPU, STLPath );
	STLCPU.ox = ((Info.cellCountX - 1) * Info.res) * 0.5f;
	STLCPU.oy = ((Info.cellCountY - 1) * Info.res) * 0.5f;
	STLCPU.oz = ((Info.cellCountZ - 1) * Info.res) * 0.5f;
	
	STLStruct STL( STLCPU );
	
	checkSTLEdges( STL );
	
	/*
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
	*/
}


/*
int main(int argc, char **argv)
{
	STLArbeiterStructCPU STLArbeiterCPU;
	
	readSTL( STLArbeiterCPU, STLPath );
	
	std::cout << "Sizing domain around the STL" << std::endl;
	CellCountStruct cellCount;
	cellCount.res = res;
	cellCount.nx = static_cast<size_t>(std::ceil((STLArbeiterCPU.xmax - STLArbeiterCPU.xmin) / cellCount.res));
	cellCount.ny = static_cast<size_t>(std::ceil((STLArbeiterCPU.ymax - STLArbeiterCPU.ymin) / cellCount.res));
	cellCount.nz = static_cast<size_t>(std::ceil((STLArbeiterCPU.zmax - STLArbeiterCPU.zmin) / cellCount.res));
	
	cellCount.ox = STLArbeiterCPU.xmin + ( 0.5f * ( (STLArbeiterCPU.xmax - STLArbeiterCPU.xmin) - cellCount.res * (cellCount.nx-1) ) );
	cellCount.oy = STLArbeiterCPU.ymin + ( 0.5f * ( (STLArbeiterCPU.ymax - STLArbeiterCPU.ymin) - cellCount.res * (cellCount.ny-1) ) );
	cellCount.oz = STLArbeiterCPU.zmin + ( 0.5f * ( (STLArbeiterCPU.zmax - STLArbeiterCPU.zmin) - cellCount.res * (cellCount.nz-1) ) );
	
	cellCount.n = cellCount.nx * cellCount.ny * cellCount.nz;
	std::cout << "	nx: " << cellCount.nx << "\n";
    std::cout << "	ny: " << cellCount.ny << "\n";
    std::cout << "	nz: " << cellCount.nz << "\n";
    std::cout << "	n: " << cellCount.n << "\n";
	std::cout << "	ox: " << cellCount.ox << " mm \n";
    std::cout << "	oy: " << cellCount.oy << " mm \n";
    std::cout << "	oz: " << cellCount.oz << " mm \n";
    
    STLArbeiterStruct STLArbeiter; 
    STLArbeiter.axArray = STLArbeiterCPU.axArray;
    STLArbeiter.ayArray = STLArbeiterCPU.ayArray;
    STLArbeiter.azArray = STLArbeiterCPU.azArray; 
	STLArbeiter.bxArray = STLArbeiterCPU.bxArray;
    STLArbeiter.byArray = STLArbeiterCPU.byArray;
    STLArbeiter.bzArray = STLArbeiterCPU.bzArray; 
    STLArbeiter.cxArray = STLArbeiterCPU.cxArray;
    STLArbeiter.cyArray = STLArbeiterCPU.cyArray;
    STLArbeiter.czArray = STLArbeiterCPU.czArray; 
    STLArbeiter.xmin = STLArbeiterCPU.xmin;
    STLArbeiter.ymin = STLArbeiterCPU.ymin;
    STLArbeiter.zmin = STLArbeiterCPU.zmin; 
    STLArbeiter.xmax = STLArbeiterCPU.xmax;
    STLArbeiter.ymax = STLArbeiterCPU.ymax;
    STLArbeiter.zmax = STLArbeiterCPU.zmax; 
    STLArbeiter.triangleCount = STLArbeiterCPU.triangleCount;
    
	checkSTLEdges( STLArbeiter );
    
    MarkerStruct Marker;
	Marker.fluidArray = MarkerArrayType( cellCount.n, 0);
	Marker.bouncebackArray = MarkerArrayType( cellCount.n, 0);
	Marker.givenRhoArray = MarkerArrayType( cellCount.n, 0);
	Marker.givenUxUyUzArray = MarkerArrayType( cellCount.n, 0);
	
	applyMarkersFromSTL( Marker, STLArbeiter, cellCount );
	
	DistributionStruct F;
	F.shifter = IndexArrayType( 27, 0 );
	F.fArray.setSizes( 27, cellCount.n );
	F.fArray.setValue( 1.0f );
	
	std::cout << "Marking cells" << std::endl;
	applyMarkers(Marker, cellCount);
	
	std::cout << "Filling F" << std::endl;
	applyInitialization( F, cellCount);
	
	std::cout << "Starting simulation" << std::endl;
	
	const size_t iCut = (size_t)cellCount.nx / 2;
	int plotNumber = 0;
	
	TNL::Timer lapTimer;
	lapTimer.start();
	for (int i=0; i<iterationCount; i++)
	{
		applyStreaming( F, cellCount );
		applyLocalCellUpdate( Marker, F, cellCount );
		
		if (i%100 == 0 && i!=0)
		{
			lapTimer.stop();
			std::cout << "Finished iteration " << i << std::endl;
			auto lapTime = lapTimer.getRealTime();
			float glups = (cellCount.n * 100) / lapTime / 1000000000;
			std::cout << "GLUPS: " << glups << std::endl;
			lapTimer.reset();
			
			exportSectionCutPlot( Marker, F, cellCount, iCut, plotNumber );
			plotNumber++;
			
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
*/
