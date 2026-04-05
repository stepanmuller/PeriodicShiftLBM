constexpr int caseID = 1;

constexpr float widthGlobal = 150.f; 													// mm
constexpr float zMinGlobal = -150.f; 													// mm
constexpr float zMaxGlobal = 60.f; 														// mm
constexpr float yMinGlobal = -80.f; 													// mm
constexpr float yMaxGlobal = 20.f; 														// mm

constexpr float resGlobal = 1.6f; 														// mm
constexpr int gridLevelCount = 5;
constexpr int wallRefinementSpan = 0;

constexpr int iterationCount = 1000000;
constexpr int iterationChunk = 5000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uzInlet = 0.02f; 														// also works as nominal LBM Mach number
constexpr float hullAngle = 0.f;														// degrees
constexpr float uyInlet = (2.f * 3.14159f * hullAngle / 360.f) * uzInlet;				// this gives hull angle		
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 20.f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float pOutletPhys = 30000.f;													// Pa
constexpr float rhoOutlet = pOutletPhys / ( rhoNominalPhys * soundspeedPhys * soundspeedPhys ) + 1.f;

#include "../include/types.h"

#include "../include/cellFunctions.h"
#include "../include/applyStreaming.h"
#include "../include/applyCollision.h"

#include "../include/boundaryConditions/applyBounceback.h"
#include "../include/boundaryConditions/applyMirror.h"
#include "../include/boundaryConditions/restoreRho.h"
#include "../include/boundaryConditions/restoreUxUyUz.h"
#include "../include/boundaryConditions/restoreRhoUxUyUz.h"
#include "../include/boundaryConditions/applyMBBC.h"

#include "../include/STLFunctions.h"
std::string STLPathIntake = "IntakeSTL3.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float zPhys = kCell * Info.res + Info.oz;
	
	if ( Info.gridID == 0 )
	{
		if ( fabs(xPhys) < 25.f && yPhys > -35.f ) Marker.refinement = 1;
	}
	// Reduce the refinement area
	if ( Info.gridID == 1 )
	{
		if ( fabs(xPhys) > 35.f ) Marker.refinement = 0;
		if ( fabs(xPhys) < 20.f && yPhys > -34.f ) Marker.refinement = 1;
	}
	if ( Info.gridID == 2 )
	{
		if ( fabs(xPhys) > 25.f ) Marker.refinement = 0;
		if ( fabs(xPhys) < 18.f && yPhys > -30.f ) Marker.refinement = 1;
	}
	if ( Info.gridID == 3 )
	{
		if ( fabs(xPhys) > 17.f ) Marker.refinement = 0;
	}
	
	if ( Marker.bounceback ) return;
	if ( kCell == 0 || jCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.givenUxUyUz = 1;
	else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	ux = 0.f;
	uy = uyInlet;
	uz = uzInlet;
	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float r2 = xPhys * xPhys + yPhys * yPhys;
	if ( r2 < 26.f * 26.f ) rho = rhoOutlet; // intake
	else rho = 1.f; // lake
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
	if (Marker.bounceback) uz = 0.f;
}

#include "../include/applyLocalCellUpdate.h"
#include "../include/plotter/exportSectionCutPlot.h"
#include "../include/fillEquilibrium.h"
#include "../include/gridRefinementFunctions.h"
#include "../include/DIADFunctions.h"
#include "../include/flowReportFunctions.h"

void updateGrid( std::vector<DIADGridStruct>& grids, int level ) 
{
    applyStreaming(grids[level]);
    applyLocalCellUpdate(grids[level]);
    if (level < gridLevelCount - 1) 
    {
        for (int i = 0; i < 2; i++) updateGrid(grids, level + 1);
        applyCoarseFineGridCommunication(grids[level], grids[level + 1]);
    }
}

void exportHistoryData( const std::vector<float>& historyMassFlow, const std::vector<float>& historyEta, const int &currentIteration, int fileNumber ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyMassFlow.data(), sizeof(float), count, fp);
    fwrite(historyEta.data(), sizeof(float), count, fp); // Append the ETA array
    fclose(fp);
    
    // Construct the command string to pass the number as an argument
    std::string cmd = "python3 historyPlotter.py " + std::to_string(fileNumber) + " &";
    system(cmd.c_str());
}

int main(int argc, char **argv)
{
	STLStructCPU STLCPUIntake;
	readSTL( STLCPUIntake, STLPathIntake );
	STLStruct STLIntake( STLCPUIntake );
	checkSTLEdges( STLIntake );
	
	std::vector<STLStruct> STLs = { STLIntake };
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = (int)( widthGlobal / grids[0].Info.res );
	grids[0].Info.cellCountY = (int)( (yMaxGlobal - yMinGlobal) / grids[0].Info.res );
	grids[0].Info.cellCountZ = (int)( (zMaxGlobal - zMinGlobal) / grids[0].Info.res );
	grids[0].Info.ox = 0.5f * ( - grids[0].Info.res * grids[0].Info.cellCountX ) + 0.5f * grids[0].Info.res;
	grids[0].Info.oy = yMinGlobal;
	grids[0].Info.oz = zMinGlobal;
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
	
	std::vector<float> historyMassFlow( iterationCount, 0.f );
	std::vector<float> historyEta( iterationCount, 0.f );
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		if (iteration%20 == 0 && iteration != 0)
		{
			int iTemp, jTemp, kCutBelow;
			const float xTemp = 0.f; 
			const float yTemp = 0.f;
			const float zBelow = -100.f;
			getIJKCellIndexFromXYZ( iTemp, jTemp, kCutBelow, xTemp, yTemp, zBelow, grids[gridLevelCount-1].Info);
			FlowReportStruct FlowReportBelow;
			XYZBoundsStruct BoundsBelow;
			BoundsBelow.xmin = -30.f;
			BoundsBelow.xmax = 30.f;
			BoundsBelow.ymin = -60.f;
			BoundsBelow.ymax = -25.f;
			getFlowReportXY( grids, kCutBelow, BoundsBelow, FlowReportBelow );
			
			int kCutAbove;
			const float zAbove = 30.f;
			getIJKCellIndexFromXYZ( iTemp, jTemp, kCutAbove, xTemp, yTemp, zAbove, grids[gridLevelCount-1].Info);
			FlowReportStruct FlowReportAbove;
			XYZBoundsStruct BoundsAbove;
			BoundsAbove.xmin = -20.f;
			BoundsAbove.xmax = 20.f;
			BoundsAbove.ymin = -20.f;
			BoundsAbove.ymax = 20.f;
			getFlowReportXY( grids, kCutAbove, BoundsAbove, FlowReportAbove );
			
			convertToPhysicalVelocity( FlowReportBelow.uz, FlowReportAbove.uz, FlowReportAbove.uy, grids[gridLevelCount-1].Info );
			float pAvg = FlowReportAbove.rho;
			float pAvgLake = FlowReportBelow.rho;
			convertToPhysicalPressure( pAvg );
			convertToPhysicalPressure( pAvgLake );
			
			const float massFlow = (FlowReportAbove.areamm2 / 1000000.f) * FlowReportAbove.uz * rhoNominalPhys;

			float lakePower = 0.5f * massFlow * FlowReportBelow.uz * FlowReportBelow.uz;
			float intakePower = 0.5f * massFlow * FlowReportAbove.uz * FlowReportAbove.uz + massFlow * (pAvg - pAvgLake) / rhoNominalPhys;
			float eta = intakePower / lakePower;
			
			for ( int shifter = 0; shifter < 20; shifter++ )
			{
				historyMassFlow[iteration-shifter] = massFlow;
				historyEta[iteration-shifter] = eta;
			}
		}
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			int iCut, jCut, kCut;
			float xCut = 0.f; 
			float yCut = 0.f;
			float zCut = 0.f;
			
			xCut = 0.f; 
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZY( grids, iCut, iteration );
			system("python3 ../include/plotter/plotterGridID.py");
			
			xCut = 5.f; 
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZY( grids, iCut, iteration + 1 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			xCut = 10.f; 
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZY( grids, iCut, iteration + 2 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			zCut = -30.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotXY( grids, kCut, iteration + 10 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			zCut = 0.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotXY( grids, kCut, iteration + 11 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			zCut = 30.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotXY( grids, kCut, iteration + 12 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			yCut = -40.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZX( grids, jCut, iteration + 20 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			yCut = -30.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZX( grids, jCut, iteration + 21 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			yCut = -20.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZX( grids, jCut, iteration + 22 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			yCut = -10.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZX( grids, jCut, iteration + 23 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			yCut = 0.f;
			getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[gridLevelCount-1].Info);
			exportSectionCutPlotZX( grids, jCut, iteration + 24 );
			system("python3 ../include/plotter/plotterGridID.py");
			
			exportHistoryData( historyMassFlow, historyEta, iteration, caseID );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
