// input physics
constexpr float res = 0.35f; 											// mm
constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number

constexpr float nuPhys = 1e-6;											// m2/s water
constexpr float rhoNominalPhys = 1000.0f;								// kg/m3 water
constexpr float SmagorinskyConstantGlobal = 0.1f; 						// set to zero to turn off LES

// calculated from input
constexpr float uzInletPhys = 17.f; 									// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * (res/1000); 			// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (res/1000) / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / ((res/1000) * (res/1000));		// LBM nu
constexpr float tau = 3.f * nu + 0.5f;		

constexpr float massFlowDesired = 4.0f; 								// kg/s in the intake
constexpr float iRegulatorSensitivity = 2e-8;							// how fast outlet pressure BC gets adjusted to achieve correct mass flow
constexpr float pRegulatorSensitivity = 1e-4;

constexpr int iterationCount = 1000001;
constexpr int iterationChunk = 100;

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

std::string STLPath = "M35IntakeLoaderSTL.STL";

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
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
	const int cellsInRadius = (int)(17.5f / Info.res) + 1;
	const int iMin = ((int)Info.cellCountX / 2) - cellsInRadius;
	const int iMax = ((int)Info.cellCountX / 2) + cellsInRadius; // both included
	const int jMin = Info.cellCountY-1 - 2*cellsInRadius;
	const int jMax = Info.cellCountY-1; // both included
	if ( iCell >= iMin && iCell <= iMax && jCell >= jMin && jCell <= jMax ) rho = 1.f + Info.iRegulator + Info.pRegulator; // intake
	else rho = 1.f; // lake
}

__cuda_callable__ float getSmagorinskyConstant( const int& iCell, const int& jCell, const int& kCell, const InfoStruct &Info )
{
	const float xToEnd = (Info.cellCountX-1 - iCell) * Info.res;
	if ( xToEnd > 60.f ) return SmagorinskyConstantGlobal;
	else
	{
		return 1.f; // set Smagorinsky high to dampen vortices before reaching the outlet
	}
}

#include "../applyLocalCellUpdate.h"

float getMassFlowOut( FStruct &F, BoolArrayType &bouncebackArray, InfoStruct &Info )
{
	auto fArrayView  = F.fArray.getView();
	auto shifterView  = F.shifter.getConstView();
	auto bouncebackArrayView = bouncebackArray.getConstView();
	
	const int cellsInRadius = (int)(17.5f / Info.res) + 1;
	const int iMin = ((int)Info.cellCountX / 2) - cellsInRadius;
	const int iMax = ((int)Info.cellCountX / 2) + cellsInRadius; // both included
	const int jMin = Info.cellCountY-1 - 2*cellsInRadius;
	const int jMax = Info.cellCountY-1; // both included
	
	const int spanX = iMax - iMin + 1;
	const int spanY = jMax - jMin + 1;
	
	const int start = 0;
	const int end = spanX * spanY;
	
	auto fetch = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = iMin + singleIndex % spanX;
		const int jCell = jMin + singleIndex / spanX;
		const int kCell = Info.cellCountZ-1;
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		bool bouncebackMarker = bouncebackArrayView( cell );
		if ( bouncebackMarker ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float uz = 0;
		float rho, ux, uy;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		float p;
		convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
		return uz;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float uzPhysSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.f );
	float volumetricFlow = uzPhysSum * (Info.res / 1000.f) * (Info.res / 1000.f);
	float massFlow = volumetricFlow * 1000.f;
	return massFlow;
}

float getIntakePower( FStruct &F, BoolArrayType &bouncebackArray, InfoStruct &Info )
{
	auto fArrayView  = F.fArray.getView();
	auto shifterView  = F.shifter.getConstView();
	auto bouncebackArrayView = bouncebackArray.getConstView();
	
	const int cellsInRadius = (int)(17.5f / Info.res) + 1;
	const int iMin = ((int)Info.cellCountX / 2) - cellsInRadius;
	const int iMax = ((int)Info.cellCountX / 2) + cellsInRadius; // both included
	const int jMin = Info.cellCountY-1 - 2*cellsInRadius;
	const int jMax = Info.cellCountY-1; // both included
	
	const int spanX = iMax - iMin + 1;
	const int spanY = jMax - jMin + 1;
	
	const int start = 0;
	const int end = spanX * spanY;
	
	auto fetch = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = iMin + singleIndex % spanX;
		const int jCell = jMin + singleIndex / spanX;
		const int kCell = Info.cellCountZ-1;
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		bool bouncebackMarker = bouncebackArrayView( cell );
		if ( bouncebackMarker ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float uz = 0;
		float rho, ux, uy;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		float p;
		convertToPhysicalUnits( rho, p, ux, uy, uz, Info );
		float volumetricFlow = uz * (Info.res / 1000.f) * (Info.res / 1000.f);
		float volumetricPower = p * volumetricFlow;
		float massFlow = volumetricFlow * 1000.f;
		float kineticPower = 0.5f * massFlow * uz * std::abs( uz );
		float totalPower = volumetricPower + kineticPower;
		return totalPower;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float intakePower = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.f );
	return intakePower;
}

void exportRegulatorData(const std::vector<float>& regulator, const std::vector<float>& massFlow, const std::vector<float>& eta, int &currentIteration) {
    FILE* fp = fopen("/dev/shm/regulator_history.bin", "wb");
    if (!fp) return;
    
    // Header: [Number of entries recorded so far]
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    
    // Write the three arrays
    fwrite(regulator.data(), sizeof(float), count, fp);
    fwrite(massFlow.data(), sizeof(float), count, fp);
    fwrite(eta.data(), sizeof(float), count, fp);
    
    fclose(fp);
    // Call the python plotter (non-blocking if possible, but system() is fine for now)
    system("python3 regulatorPlotter.py &"); 
}

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
	
	fillDefaultEquilibrium( F, Info, 1.0f, 0.f, 0.f, uzInlet );
	
	std::cout << "Starting simulation" << std::endl;
	
	const int iCut = Info.cellCountX / 2;
	int plotNumber = 0;
	
	std::vector<float> historyRegulator( iterationCount, 0.f );
	std::vector<float> historyMassFlow( iterationCount, 0.f );
	std::vector<float> historyIntakeEta( iterationCount, 0.f );
	float movingError = 0.f;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<iterationCount; iteration++)
	{
		applyStreaming( F, Info );
		applyLocalCellUpdate( F, bouncebackArray, Info );
		
		float rhoPrevious = 1.0f + Info.iRegulator;
		
		float massFlow = getMassFlowOut( F, bouncebackArray, Info );
		float error = massFlowDesired - massFlow;
		Info.iRegulator = Info.iRegulator - iRegulatorSensitivity * error;
		Info.iRegulator = std::clamp( Info.iRegulator, -0.01f, 0.01f );
		movingError = movingError * 0.999f + 0.001f * error;
		Info.pRegulator = - pRegulatorSensitivity * movingError;

		float intakePower = getIntakePower( F, bouncebackArray, Info );
		float lakePower = 0.5f * massFlow * uzInletPhys * uzInletPhys;
		float intakeEta = std::max({0.f, intakePower}) / lakePower;
		
		historyRegulator[iteration] = Info.iRegulator + Info.pRegulator;
		historyMassFlow[iteration] = massFlow;
		historyIntakeEta[iteration] = intakeEta;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			std::cout << "Outlet rho adjusted to: " << 1.0f + Info.iRegulator + Info.pRegulator << std::endl;
			std::cout << "Mass flow out: " << massFlow << " kg/s" << std::endl; 
			std::cout << "Intake eta: " << intakeEta << std::endl;
			
			float glups = ((float)Info.cellCount * (float)iterationChunk) / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlotZY( F, bouncebackArray, Info, iCut, plotNumber );
			plotNumber++;
			
			exportRegulatorData( historyRegulator, historyMassFlow, historyIntakeEta, iteration );
			
			std::cout << std::endl;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
