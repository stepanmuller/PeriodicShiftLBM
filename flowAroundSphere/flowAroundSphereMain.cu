constexpr float sphereDiameterPhys = 2.f;										// m
constexpr float res = 50.f; 													// mm
constexpr float uxInlet = 0.07f; 												// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 1000.f;
constexpr float SmagorinskyConstantGlobal = 0.0f; 								// set to zero to turn off LES

constexpr float sphereRadiusPhys = 0.5 * sphereDiameterPhys;					// m
constexpr float uxInletPhys = uxInlet; 											// m/s, physical velocity set to same as LBM velocity
constexpr float nuPhys = uxInletPhys * sphereDiameterPhys / reynoldsNumber;		// m2/s
constexpr float rhoNominalPhys = 1.225f;										// kg/m3 air
constexpr float dtPhys = (uxInlet / uxInletPhys) * (res/1000); 					// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (res/1000) / dtPhys; 				// m/s
constexpr float nu = (dtPhys * nuPhys) / ((res/1000) * (res/1000));				// LBM nu
constexpr float tau = 3.f * nu + 0.5f;											// LBM tau

constexpr float domainSizePhys = 11.f * sphereDiameterPhys;						// m
constexpr float sphereXPhys = 2.f * sphereDiameterPhys;							// m
constexpr float sphereYPhys = 5.5f * sphereDiameterPhys;						// m
constexpr float sphereZPhys = 5.5f * sphereDiameterPhys;						// m

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys * 1000.f / res));
const int cellCountY = cellCountX;
const int cellCountZ = cellCountX;

constexpr int iterationCount = 100000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/applyMirror.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& mirrorMarker, bool& periodicMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									const InfoStruct& Info )
{
    fluidMarker = 0;
	bouncebackMarker = 0;
	mirrorMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	const float xPhys = iCell * Info.res * 0.001f;
	const float yPhys = jCell * Info.res * 0.001f;
	const float zPhys = kCell * Info.res * 0.001f;
	
	const float r2 = (xPhys-sphereXPhys) * (xPhys - sphereXPhys) + (yPhys - sphereYPhys) * (yPhys - sphereYPhys) + (zPhys - sphereZPhys) * (zPhys - sphereZPhys);
		
	if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) bouncebackMarker = 1;
	else if ( jCell == 0 || jCell == Info.cellCountY-1 ) givenUxUyUzMarker = 1;
	else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) givenUxUyUzMarker = 1;
	else if ( iCell == 0 ) givenUxUyUzMarker = 1;
	else if ( iCell == Info.cellCountX-1 ) givenRhoMarker = 1;
	else fluidMarker = 1;
	
	periodicMarker = 0;
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

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const InfoStruct &Info )
{
	rho = 1.f;
	ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
	bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
	getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Info );
	if ( bouncebackMarker ) ux = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../exportSectionCutPlot.h"
#include "../fillDefaultEquilibrium.h"

float getSphereDrag( FStruct &F, const InfoStruct &Info )
{
	auto fArrayView  = F.fArray.getView();
	auto shifterView  = F.shifter.getConstView();
	
	const int iStart = (int)((sphereXPhys - sphereRadiusPhys) / (Info.res / 1000));
	const int iEnd = (int)((sphereXPhys + sphereRadiusPhys) / (Info.res / 1000)) + 2;
	const int jStart = (int)((sphereYPhys - sphereRadiusPhys) / (Info.res / 1000));
	const int jEnd = (int)((sphereYPhys + sphereRadiusPhys) / (Info.res / 1000)) + 2;
	const int kStart = (int)((sphereZPhys - sphereRadiusPhys) / (Info.res / 1000));
	const int kEnd = (int)((sphereZPhys + sphereRadiusPhys) / (Info.res / 1000)) + 2;
	
	const int start = 0;
	const int end = (iEnd - iStart) * (jEnd - jStart) * (kEnd - kStart);
	
	auto fetch = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iSpan = iEnd - iStart;
		const int jSpan = jEnd - jStart;
		const int kRelative = singleIndex / (iSpan * jSpan);
		const int remainder = singleIndex % (iSpan * jSpan);
		const int jRelative = remainder / iSpan;
		const int iRelative = remainder % iSpan;
		
		const int iCell = iRelative + iStart;
		const int jCell = jRelative + jStart;
		const int kCell = kRelative + kStart;
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
		getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Info );
		
		if ( !bouncebackMarker ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz( rho, ux, uy, uz, f );
		applyBounceback(f);
		float rhoPrev, uxPrev, uyPrev, uzPrev;
		getRhoUxUyUz( rhoPrev, uxPrev, uyPrev, uzPrev, f );
		const float gx = rho * (ux -  uxPrev);
		return gx;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float gxSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.f );
	float gy, gz = 0;
	convertToPhysicalForce( gxSum, gy, gz, Info );
	return gxSum;
}

void exportHistoryData( const std::vector<float>& historyDragCoefficient, const int &currentIteration ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyDragCoefficient.data(), sizeof(float), count, fp);
    fclose(fp);
    system("python3 historyPlotter.py &"); 
}

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
	
	FStruct F;
	FloatArray2DType fArray;
	F.fArray.setSizes( 27, Info.cellCount );	
	F.shifter = IntArrayType( 27, 0 );
	
	fillDefaultEquilibriumFromFunction( F, Info );
	
	std::cout << "Starting simulation" << std::endl;
	
	const int kCut = Info.cellCountZ / 2;
	
	const int iterationChunk = 1000;
	
	std::vector<float> historyDragCoefficient( iterationCount, 0.f );
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyStreaming( F, Info );
		applyLocalCellUpdate( F, Info );
		const float drag = getSphereDrag( F, Info );
		const float dragCoefficient = - (8 * drag) / (Info.rhoNominalPhys * uxInletPhys * uxInletPhys * 3.14159f * sphereDiameterPhys * sphereDiameterPhys);
		
		historyDragCoefficient[iteration] = dragCoefficient;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			std::cout << "Drag " << drag << std::endl;
			std::cout << "Drag coefficient " << dragCoefficient << std::endl;
			
			const int cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
			float glups = ((float)cellCount * (float)iterationChunk) / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlotXY( F, Info, kCut, iteration );
			
			exportHistoryData( historyDragCoefficient, iteration );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
