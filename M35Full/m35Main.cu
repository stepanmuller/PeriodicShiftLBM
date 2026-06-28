constexpr float resGlobal = 0.4f; 														// mm
constexpr int gridLevelCount = 3;
constexpr int wallRefinementSpan = 1;

constexpr int iterationCount = 500000;
constexpr int iterationChunk = 1000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number	
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 4.5986f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float RIn = 3.75f;															// mm
constexpr float ROut = 16.5f;															// mm
constexpr float angularVelocity = 2700.f;												// rad/s
const float boundaryLayerThickness = 0.2f;												// mm

#include "../include/types.h"

#include <sstream>

#include "../include/cellFunctions.h"
#include "../include/applyStreaming.h"
//#include "../include/applyCollision.h"
#include "./applyCollisionCustomized.h"

#include "../include/boundaryConditions/applyBounceback.h"
#include "../include/boundaryConditions/applyMovingBounceback.h"
#include "../include/boundaryConditions/restoreRho.h"
#include "../include/boundaryConditions/restoreUxUyUz.h"
#include "../include/boundaryConditions/restoreRhoUxUyUz.h"
#include "../include/boundaryConditions/applyMBBC.h"
#include "../include/boundaryConditions/applyNonReflectiveOutlet.h"

#include "../include/STLFunctions.h"
std::string STLPathMain = "M-Jet_35_pump_main.STL";
std::string STLPathImpeller = "M-Jet_35_impeller.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	if ( Marker.forcedVelocity ) return;
	if ( kCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info, MarkerStruct &Marker )
{
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r = std::sqrt( x * x + y * y );
	const float vtPhys = angularVelocity * (r / 1000.f);
	const float vt = vtPhys * ( uzInlet / uzInletPhys );
	const float wallDistancePhys = std::max(0.f, std::min(r - RIn, ROut - r));
	const float delta = std::max( 0.f, std::min( 1.f, wallDistancePhys / boundaryLayerThickness ));
	const float velocityMultiplier = delta * delta * (3.0f - 2.0f * delta);
	if ( Marker.forcedVelocity )
	{
		ux = - vt * (y / r);
		uy = vt * (x / r);
		uz = 0.f;
	}
	else
	{
		ux = 0.f;
		uy = 0.f;
		uz = uzInlet * velocityMultiplier;
	}
	rho = 1.f;
}

__cuda_callable__ void getForcing( 	const int& iCell, const int& jCell, const int& kCell, 
									const float (&f)[27], 
									float& gx, float& gy, float& gz,
									InfoStruct& Info )
{
	float rho, ux, uy, uz;
	getRhoUxUyUz( rho, ux, uy, uz, f );
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r = std::sqrt( x * x + y * y );
	const float vtPhys = angularVelocity * (r / 1000.f);
	const float vt = vtPhys * ( uzInlet / uzInletPhys );
	const float uxTarget = - vt * (y / r);
	const float uyTarget = vt * (x / r);
	const float uzTarget = 0.f;
	gx = rho * ( uxTarget - ux );
	gy = rho * ( uyTarget - uy );
	gz = rho * ( uzTarget - uz );
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
	rho = 1.f;
}

//#include "../include/applyLocalCellUpdate.h"
#include "./applyLocalCellUpdateCustomized.h"
#include "./exportSectionCutPlotCustomized.h"
#include "../include/fillEquilibrium.h"
#include "../include/gridRefinementFunctions.h"
#include "../include/DIADFunctions.h"
#include "../include/flowReportFunctions.h"

void updateGrid( std::vector<DIADGridStruct>& grids, int level ) 
{
	//applyNonReflectiveOutletZ(grids[level]);
    applyStreaming(grids[level]);
    applyLocalCellUpdate(grids[level]);
    if (level < gridLevelCount - 1) 
    {
        for (int i = 0; i < 2; i++) updateGrid(grids, level + 1);
        applyCoarseFineGridCommunication(grids[level], grids[level + 1]);
    }
}

void exportHistoryData( const std::vector<float>& historyVector, const int &currentIteration, int fileNumber ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyVector.data(), sizeof(float), count, fp);
    fclose(fp);
    // Construct the command string to pass the number as an argument
    std::string cmd = "python3 historyPlotter.py " + std::to_string(fileNumber) + " &";
    system(cmd.c_str());
}

int main(int argc, char **argv)
{
	int caseID = 0;
    if (argc > 1) 
    {
		caseID = std::stoi(argv[1]);
	}

	STLStructCPU STLCPUMain;
	readSTL( STLCPUMain, STLPathMain );
	STLStruct STLMain( STLCPUMain );
	checkSTLEdges( STLMain );
	
	std::vector<STLStruct> STLs = { STLMain };
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.ox = STLMain.xmin + 0.5f * grids[0].Info.res;
	grids[0].Info.oy = STLMain.ymin + 0.5f * grids[0].Info.res;
	grids[0].Info.oz = STLMain.zmin + 0.5f * grids[0].Info.res;
	grids[0].Info.cellCountX = (int)( (STLMain.xmax - STLMain.xmin) / grids[0].Info.res );
	grids[0].Info.cellCountY = (int)( (STLMain.ymax - STLMain.ymin) / grids[0].Info.res );
	grids[0].Info.cellCountZ = (int)( (STLMain.zmax - STLMain.zmin) / grids[0].Info.res );
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	buildDIADGrids( grids, STLs, 0 );
	
	int cellCountTotal = 0;
	long int cellUpdatesPerIteration = 0;
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		const int cellCountLevel = grids[level].Info.cellCount;
		cellCountTotal += cellCountLevel; 
		cellUpdatesPerIteration += cellCountLevel * std::pow(2, level);
	}
	std::cout << "Cell count total: " << cellCountTotal << std::endl;
	std::cout << "Cell updates per iteration: " << cellUpdatesPerIteration << std::endl;	
	
	std::vector<float> historyVector( iterationCount, 0.f );
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		if (iteration%20 == 0 && iteration != 0)
		{
			XYZBoundsStruct Bounds;
			Bounds.xmin = -ROut;
			Bounds.xmax = ROut;
			Bounds.ymin = -ROut;
			Bounds.ymax = ROut;
			
			FlowReportStruct FlowReportIn;
			int kCut = 0;
			getFlowReportXY( grids, kCut, Bounds, FlowReportIn );
			float pIn = FlowReportIn.rho;
			float uzIn = FlowReportIn.uz;
			convertToPhysicalPressure( pIn );
			float uTemp = 0.f;
			convertToPhysicalVelocity( uzIn, uTemp, uTemp, grids[0].Info );
			
			float pTotalIn = 0.5f * rhoNominalPhys * uzIn * uzIn + pIn;
			
			for ( int shifter = 0; shifter < 20; shifter++ )
			{
				historyVector[iteration-shifter] = pTotalIn;
			}
			
		}
		
		if ( iteration % iterationChunk == 0) 
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportHistoryData( historyVector, iteration, 0 );
			
			int counter = 1;
			for (float r = RIn + 1.f; r < ROut; r = r + 1.f) 
			{
				exportSectionCutPlotToiletPaperZ( grids, r, iteration + counter );
				system("python3 ../include/plotter/plotterGridID.py");
				counter++;
			}	
			counter = 1;
			for (float z = grids[gridLevelCount-1].Info.oz; z < ((grids[gridLevelCount-1].Info.cellCountZ-1) * grids[gridLevelCount-1].Info.res + grids[gridLevelCount-1].Info.oz); z = z + 5.f) 
			{
				int iTemp, jTemp, kCut;
				const float xTemp = 0;
				const float yTemp = 0;
				getIJKCellIndexFromXYZ( iTemp, jTemp, kCut, xTemp, yTemp, z, grids[gridLevelCount-1].Info);
				exportSectionCutPlotXY( grids, kCut, iteration + 30 + counter );
				system("python3 ../include/plotter/plotterGridID.py");
				counter++;
			}	
			const int kCut = grids[gridLevelCount-1].Info.cellCountZ / 2;
			exportSectionCutPlotXY( grids, kCut, iteration + 30 + counter + 1 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			const int iCut = grids[gridLevelCount-1].Info.cellCountX / 2;
			exportSectionCutPlotZY( grids, iCut, iteration + 60 );
			system("python3 ../include/plotter/plotterGridID.py");

			lapTimer.reset();
			lapTimer.start();	
		}
	}
	return EXIT_SUCCESS;
}
