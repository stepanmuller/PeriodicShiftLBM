constexpr float resGlobal = 0.4f; 														// mm
constexpr int gridLevelCount = 3;
constexpr int wallRefinementSpan = 1;

constexpr int iterationCount = 150000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES
constexpr float SmagorinskyZoneLength = 10.f;

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number	
constexpr float rhoOutlet = 1.f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 16.94f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float RIn = 6.f;																// mm
constexpr float ROut = 12.f;															// mm
constexpr float C = 0.092784f;															// m2/s
const float boundaryLayerThickness = 0.2f;												// mm
constexpr int bladeCount = 5;

constexpr int iterationChunk = 2000;

#include "../include/types.h"

#include <sstream>

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
std::string STLPathBlade = "blade.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r2 = x * x + y * y;
	if ( r2 >= ROut * ROut )
	{
		Marker.bounceback = 1;
		return;
	}
	if ( r2 <= RIn * RIn )
	{
		Marker.bounceback = 1;
		return;
	}
	if ( kCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r = std::sqrt( x * x + y * y );
	const float vtPhys = (1.f / (r / 1000.f)) * C;
	const float vt = vtPhys * ( uzInlet / uzInletPhys );
	
	const float wallDistancePhys = std::max(0.f, std::min(r - RIn, ROut - r));
	const float delta = std::max( 0.f, std::min( 1.f, wallDistancePhys / boundaryLayerThickness ));
	const float velocityMultiplier = delta * delta * (3.0f - 2.0f * delta);
	
	ux = - vt * (y / r) * velocityMultiplier;
	uy = vt * (x / r) * velocityMultiplier;
	uz = uzInlet * velocityMultiplier;
	rho = rhoOutlet + Info.regulator;
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	const float zMax = Info.oz + Info.cellCountZ * Info.res;
	const float z = Info.oz + kCell * Info.res;
	if ( zMax - z < SmagorinskyZoneLength ) return 1.f; // 10mm from the end
	else return SmagorinskyConstantGlobal;
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

void fillOutletCellList( DIADGridStruct &Grid )
{
	auto kView = Grid.IJK.kArray.getView();
	auto fetch = [ = ] __cuda_callable__( const int cell )
		{
			if ( kView[ cell ] == Grid.Info.cellCountZ-2 ) return 1;
			else return 0;
		};
		auto reduction = [] __cuda_callable__( const int& a, const int& b )
		{
			return a + b;
		};	
	int outletCellCount = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, Grid.Info.cellCount, fetch, reduction, 0 );
	
	IntArrayTypeCPU kArrayCPU;
	kArrayCPU = Grid.IJK.kArray;
	IntArrayTypeCPU outletCellArrayCPU;
	outletCellArrayCPU.setSize( outletCellCount );
	int counter = 0;
	for ( int cell = 0; cell < Grid.Info.cellCount; cell++ ) 
	{
		if ( kArrayCPU[ cell ] == Grid.Info.cellCountZ-2 ) outletCellArrayCPU[ counter++ ] = cell;
	}
	Grid.outletCellArray = outletCellArrayCPU;
}

void applyNonReflectiveOutletBC( DIADGridStruct &Grid )
{
	const int outletCellCount = Grid.outletCellArray.getSize();
	
	auto outletCellView = outletCellArray.getView();
	
	auto fArrayView  = Grid.fArray.getView();
	
	auto iView = Grid.IJK.iArray.getConstView();
	auto jView = Grid.IJK.jArray.getConstView();
	auto kView = Grid.IJK.kArray.getConstView();
	
	bool esotwistFlipper = Grid.esotwistFlipper;
	auto iNbrView = Grid.EsotwistNbrArray.iNbrArray.getConstView();
	auto jNbrView = Grid.EsotwistNbrArray.jNbrArray.getConstView();
	auto kNbrView = Grid.EsotwistNbrArray.kNbrArray.getConstView();
	auto ijNbrView = Grid.EsotwistNbrArray.ijNbrArray.getConstView();
	auto ikNbrView = Grid.EsotwistNbrArray.ikNbrArray.getConstView();
	auto jkNbrView = Grid.EsotwistNbrArray.jkNbrArray.getConstView();
	auto ijkNbrView = Grid.EsotwistNbrArray.ijkNbrArray.getConstView();
	
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	InfoStruct Info = Grid.Info;
	
	auto cellLambda = [=] __cuda_callable__ ( const int index ) mutable
	{
		const int cell = 
		const int iCell = iView( cell );
		const int jCell = jView( cell );
		const int kCell = kView( cell );
		
		DIADEsotwistNbrStruct Nbr;
		Nbr.i = iNbrView( cell );
		Nbr.j = jNbrView( cell );
		Nbr.k = kNbrView( cell );
		Nbr.ij = ijNbrView( cell );
		Nbr.ik = ikNbrView( cell );
		Nbr.jk = jkNbrView( cell );
		Nbr.ijk = ijkNbrView( cell ); 
		
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( Marker.bounceback )
		{
			return; // bounceback gets implicitly applied by Esotwist
		}
		
		float f[27];
		int cellReadIndex[27];
		int fReadIndex[27];
		getEsotwistReadIndex( cell, cellReadIndex, fReadIndex, Nbr, esotwistFlipper, Info );
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fReadIndex[direction], cellReadIndex[direction]);
		
		float rho, ux, uy, uz;
		
		if ( Marker.fluid )
		{
			// do nothing, just skip the else block below
		}
		else
		{
			int outerNormalX, outerNormalY, outerNormalZ;
			getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info ); 
			if ( Marker.periodicX ) outerNormalX = 0;
			if ( Marker.periodicY ) outerNormalY = 0;
			if ( Marker.periodicZ ) outerNormalZ = 0;
			getGivenRhoUxUyUz( iCell, jCell, kCell, rho, ux, uy, uz, Info );
			if ( Marker.mirror )
			{
				applyMirror( outerNormalX, outerNormalY, outerNormalZ, f );
			}
			else if ( Marker.givenRho && !Marker.givenUxUyUz )
			{
				restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			else if ( !Marker.givenRho && Marker.givenUxUyUz )
			{
				restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			else if ( !Marker.givenRho && !Marker.givenUxUyUz )
			{
				restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
		}
		const float SmagorinskyConstant = getSmagorinskyConstant( iCell, jCell, kCell, Info );
		applyCollision( f, Info.nu, SmagorinskyConstant );
		
		int cellWriteIndex[27];
		int fWriteIndex[27];
		getEsotwistWriteIndex( cell, cellWriteIndex, fWriteIndex, Nbr, esotwistFlipper, Info );
		
		for ( int direction = 0; direction < 27; direction++ ) fArrayView( fWriteIndex[direction], cellWriteIndex[direction] ) = f[direction];
		
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, outletCellCount, cellLambda );
}

int main(int argc, char **argv)
{
	int caseID = 0;
    if (argc > 1) 
    {
		caseID = std::stoi(argv[1]);
	}

	STLStructCPU STLCPUBlade;
	readSTL( STLCPUBlade, STLPathBlade );
	STLStruct STLBlade( STLCPUBlade );
	checkSTLEdges( STLBlade );
	
	std::vector<STLStruct> STLs = { STLBlade };
	for ( int bladeIndex = 1; bladeIndex < bladeCount; bladeIndex++ )
	{	
		float radians = bladeIndex * 3.14159f * 2.f * (1.f / (float)bladeCount);
		STLStruct STLBladeRotated;
		STLBladeRotated = STLBlade;
		rotateSTLAlongZ( STLBladeRotated, radians );
		STLs.push_back( STLBladeRotated );
	}
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = (int)( (2.f*ROut) / grids[0].Info.res ) + 2;
	grids[0].Info.cellCountY = grids[0].Info.cellCountX;
	grids[0].Info.cellCountZ = (STLBlade.zmax - STLBlade.zmin + SmagorinskyZoneLength + 20.f) / grids[0].Info.res;
	grids[0].Info.ox = - (grids[0].Info.cellCountX / 2) * grids[0].Info.res;
	grids[0].Info.oy = grids[0].Info.ox;
	grids[0].Info.oz = STLBlade.zmin - 10.f;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	buildDIADGrids( grids, STLs, 0 );
	
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		fillOutletCellList( grids[level] );
	}
	
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
	float uOutletMovingAvg = uzInlet;
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		if (iteration%20 == 0 && iteration != 0)
		{
			FlowReportStruct FlowReport;
			XYZBoundsStruct Bounds;
			Bounds.xmin = -12.5f;
			Bounds.xmax = 12.5f;
			Bounds.ymin = -12.5f;
			Bounds.ymax = 12.5f;
			int kCut = 0;
			getFlowReportXY( grids, kCut, Bounds, FlowReport );
			float pAvg = FlowReport.rho;
			convertToPhysicalPressure( pAvg );	
			for ( int shifter = 0; shifter < 20; shifter++ )
			{
				historyVector[iteration-shifter] = pAvg;
			}
			
			// REGULATING OUTLET RHO TO STRONGLY DAMPEN ACUSTIC WAVES
			kCut = grids[gridLevelCount-1].Info.cellCountZ-1;
			getFlowReportXY( grids, kCut, Bounds, FlowReport );
			for ( int level = 0; level < gridLevelCount; level++ )
			{
				grids[level].Info.regulator = ( FlowReport.uz - uOutletMovingAvg ) / invSqrt3;
			}
			uOutletMovingAvg = 0.999 * uOutletMovingAvg + 0.001 * FlowReport.uz;
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
				exportSectionCutPlotToiletPaperZ( grids, r, iteration + 10*0 + counter );
				system("python3 ../include/plotter/plotterGridID.py");
				counter++;
			}	
			int kCut = 0;
			exportSectionCutPlotXY( grids, kCut, iteration + 10);
			system("python3 ../include/plotter/plotterGridID.py");
			kCut = grids[gridLevelCount-1].Info.cellCountZ-1;
			exportSectionCutPlotXY( grids, kCut, iteration + 11);
			system("python3 ../include/plotter/plotterGridID.py");
			lapTimer.reset();
			lapTimer.start();	
		}
	}
	return EXIT_SUCCESS;
}



/*

void writeCaseResult( const std::vector<float>& historyVector, const int caseID ) {
	const int averagingIntervalStart = (iterationCount / 10 * 7);
	float sum = 0.f;
	int counter = 0;
	for ( int i = averagingIntervalStart; i < iterationCount; i++ )
	{
		sum += historyVector[i];
		counter++;
	}
	const float result = sum / (float)counter;
    std::ofstream outFile("optimizerResults.txt", std::ios::app);
    outFile << caseID << "; " << result << "\n";
    outFile.close();
}
*/
