constexpr int caseID = 1;

constexpr float resGlobal = 0.84f; 														// mm
constexpr int gridLevelCount = 3;
constexpr int iterationCount = 1000000;
constexpr int iterationChunk = 10000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float hullAngle = 0.f;														// degrees
constexpr float uyInlet = (2.f * 3.14159f * hullAngle / 360.f) * uzInlet;				// this gives hull angle		
constexpr float rhoOutlet = 1.000f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 20.f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

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

#include "../STLFunctions.h"
std::string STLPathLake = "lake.STL";
std::string STLPathIntake = "intake.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	
	if ( Info.gridID == 0 ) // if zero, we are on the coarsest grid
	{
		if ( kCell == 0 || jCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.givenUxUyUz = 1;
		else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}	
	else // we are on a finer grid
	{
		if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.ghost = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.ghost = 1;
		//else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.ghost = 1;
		else if ( kCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}
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

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const InfoStruct &Info )
{
	rho = 1.f;
	ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"
#include "../gridRefinementFunctions.h"

void updateGrid( std::vector<GridStruct>& grids, int level ) 
{
    applyStreaming(grids[level]);
    applyLocalCellUpdate( grids[level] );
    if (level < gridLevelCount - 1) 
    {
        for (int i = 0; i < 2; i++) updateGrid(grids, level + 1);
        applyCoarseFineGridCommunication(grids[level], grids[level + 1]);
    }
}

void getFlowReport( const int kCell, GridStruct &Grid, const int &iStart, const int &jStart, const int &iEnd, const int &jEnd,
						float &uzAvgPhys, float &massFlowPhys, float &pPhys )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	const int iSpan = iEnd - iStart;
	const int jSpan = jEnd - jStart;
	
	const int start = 0;
	const int end = iSpan * jSpan;
	
	auto fetchCellCount = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0;
		else return 1;
	};
	auto reductionCellCount = [] __cuda_callable__( const int& a, const int& b )
	{
		return a + b;
	};
	
	auto fetchUz = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		return uz;
	};
	auto reductionUz = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	auto fetchRho = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		return rho;
	};
	auto reductionRho = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	const int cellCount = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchCellCount, reductionCellCount, 0 );
	float uzSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchUz, reductionUz, 0.f );
	float rhoSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchRho, reductionRho, 0.f );
		
	float uzAvg = uzSum / (float)cellCount;
	float rhoAvg = rhoSum / (float)cellCount;
	float ux, uy = 0;
	convertToPhysicalVelocity( ux, uy, uzAvg, Info );
	convertToPhysicalVelocity( ux, uy, uzSum, Info );
	float pAvg = rhoAvg;
	convertToPhysicalPressure( pAvg );
	
	float volumetricFlow = uzSum * (Info.res / 1000.f) * (Info.res / 1000.f);
	float massFlow = volumetricFlow * rhoAvg * rhoNominalPhys;
	
	uzAvgPhys = uzAvg;
	massFlowPhys = massFlow;
	pPhys = pAvg;
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
	STLStructCPU STLCPULake;
	readSTL( STLCPULake, STLPathLake );
	STLStruct STLLake( STLCPULake );
	checkSTLEdges( STLLake );
	STLStructCPU STLCPUIntake;
	readSTL( STLCPUIntake, STLPathIntake );
	STLStruct STLIntake( STLCPUIntake );
	checkSTLEdges( STLIntake );
	
	const float STLxminGlobal = std::min({ STLLake.xmin, STLIntake.xmin });
	const float STLxmaxGlobal = std::max({ STLLake.xmax, STLIntake.xmax });
	const float STLyminGlobal = std::min({ STLLake.ymin, STLIntake.ymin });
	const float STLymaxGlobal = std::max({ STLLake.ymax, STLIntake.ymax });
	const float STLzminGlobal = std::min({ STLLake.zmin, STLIntake.zmin });
	const float STLzmaxGlobal = std::max({ STLLake.zmax, STLIntake.zmax });
	
	std::vector<GridStruct> grids(gridLevelCount);
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	std::cout << "Sizing domain around the STL" << std::endl;
	grids[0].Info.cellCountX = (int)( (STLxmaxGlobal - STLxminGlobal) / grids[0].Info.res );
	grids[0].Info.cellCountY = (int)( (STLymaxGlobal - STLyminGlobal) / grids[0].Info.res );
	grids[0].Info.cellCountZ = (int)( (STLzmaxGlobal - STLzminGlobal) / grids[0].Info.res );
	grids[0].Info.ox = STLxminGlobal + 0.5f * ( (STLxmaxGlobal - STLxminGlobal) - grids[0].Info.res * grids[0].Info.cellCountX ) + 0.5f * grids[0].Info.res;
	grids[0].Info.oy = STLyminGlobal + 0.5f * ( (STLymaxGlobal - STLyminGlobal) - grids[0].Info.res * grids[0].Info.cellCountY ) + 0.5f * grids[0].Info.res;
	grids[0].Info.oz = STLzminGlobal + 0.5f * ( (STLzmaxGlobal - STLzminGlobal) - grids[0].Info.res * grids[0].Info.cellCountZ ) + 0.5f * grids[0].Info.res;
	grids[0].Info.cellCountY = grids[0].Info.cellCountY + 2; // adding wall layers on top
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	grids[0].fArray.setSizes( 27, grids[0].Info.cellCount );
	grids[0].shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( grids[0] );
	grids[0].bouncebackMarkerArray = BoolArrayType( grids[0].Info.cellCount, 0 );
	const bool insideMarkerValue = 0;
	for ( int scope = 0; scope < 1; scope++ )
	{
		BoolArrayType bouncebackLake = BoolArrayType( grids[0].Info.cellCount, 0 );
		BoolArrayType bouncebackIntake = BoolArrayType( grids[0].Info.cellCount, 0 );
		applyMarkersInsideSTL( bouncebackLake, STLLake, insideMarkerValue, grids[0].Info );
		applyMarkersInsideSTL( bouncebackIntake, STLIntake, insideMarkerValue, grids[0].Info );
		multiplyBoolArrays( bouncebackLake, bouncebackIntake, grids[0].bouncebackMarkerArray );
	}
	std::cout << "Cell count on grid " << 0 << ": " << grids[0].Info.cellCount << std::endl;
	
	for ( int level = 1; level < gridLevelCount; level++ )
	{
		grids[level].Info.gridID = grids[level-1].Info.gridID + 1;
		grids[level].Info.res = grids[level-1].Info.res * 0.5f;
		grids[level].Info.dtPhys = grids[level-1].Info.dtPhys * 0.5f;
		grids[level].Info.nu = (grids[level].Info.dtPhys * nuPhys) / ((grids[level].Info.res/1000.f) * (grids[level].Info.res/1000.f));
		
		float progress = (float)level / (float)(gridLevelCount-1);
		progress = std::pow( progress, 0.5f );
		
		const float xStart = (1 - progress) * grids[0].Info.ox + progress * (-20.f);
		const float xEnd = (1 - progress) * (-grids[0].Info.ox) + progress * 20.f;
		grids[level-1].Info.iSubgridStart = (int)((xStart - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.iSubgridEnd = (int)((xEnd - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f) + 1;
		
		const float yStart = (1 - progress) * grids[0].Info.oy + progress * (-45.f);
		grids[level-1].Info.jSubgridStart = (int)((yStart - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.jSubgridEnd = (int)((STLymaxGlobal - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f) + 1;
		
		grids[level-1].Info.kSubgridStart = 0;
		grids[level-1].Info.kSubgridEnd = grids[level-1].Info.cellCountZ;
		
		grids[level].Info.ox = grids[level-1].Info.ox + grids[level-1].Info.iSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oy = grids[level-1].Info.oy + grids[level-1].Info.jSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oz = grids[level-1].Info.oz + grids[level-1].Info.kSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		
		grids[level].Info.cellCountX = (grids[level-1].Info.iSubgridEnd - grids[level-1].Info.iSubgridStart) * 2;
		grids[level].Info.cellCountY = (grids[level-1].Info.jSubgridEnd - grids[level-1].Info.jSubgridStart) * 2;
		grids[level].Info.cellCountZ = (grids[level-1].Info.kSubgridEnd - grids[level-1].Info.kSubgridStart) * 2;
		grids[level].Info.cellCount = grids[level].Info.cellCountX * grids[level].Info.cellCountY * grids[level].Info.cellCountZ;
		
		grids[level].fArray.setSizes( 27, grids[level].Info.cellCount );	
		grids[level].shifter = IntArrayType( 27, 0 );
		fillEquilibriumFromFunction( grids[level] );
		
		grids[level].bouncebackMarkerArray = BoolArrayType( grids[level].Info.cellCount, 0 );
		BoolArrayType bouncebackLake = BoolArrayType( grids[level].Info.cellCount, 0 );
		BoolArrayType bouncebackIntake = BoolArrayType( grids[level].Info.cellCount, 0 );
		applyMarkersInsideSTL( bouncebackLake, STLLake, insideMarkerValue, grids[level].Info );
		applyMarkersInsideSTL( bouncebackIntake, STLIntake, insideMarkerValue, grids[level].Info );
		multiplyBoolArrays( bouncebackLake, bouncebackIntake, grids[level].bouncebackMarkerArray );
		
		std::cout << "Cell count on grid " << level << ": " << grids[level].Info.cellCount << std::endl;
	}
	
	int cellCountTotal = 0;
	long int cellUpdatesPerIteration = 0;
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		const int cellCountLevel = grids[level].Info.cellCount;
		cellCountTotal += cellCountLevel; 
		cellUpdatesPerIteration += cellCountLevel * std::pow(2, level);
		cellUpdatesPerIteration -= (grids[level].Info.iSubgridEnd - grids[level].Info.iSubgridStart - 2) 
								* (grids[level].Info.jSubgridEnd - grids[level].Info.jSubgridStart - 2) 
								* (grids[level].Info.kSubgridEnd - grids[level].Info.kSubgridStart - 2)
								* std::pow(2, level);
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
		const int iStart = (-26.f - grids[gridLevelCount-1].Info.ox) / grids[gridLevelCount-1].Info.res;
		const int iEnd = (26.f - grids[gridLevelCount-1].Info.ox) / grids[gridLevelCount-1].Info.res;
		const int jStart = (-26.f - grids[gridLevelCount-1].Info.oy) / grids[gridLevelCount-1].Info.res;
		const int jEnd = grids[gridLevelCount-1].Info.cellCountY-1;
		
		float uzAvg, massFlow, pAvg;
		int iTemp, jTemp, flowReportK;
		const float xTemp = 0.f; 
		const float yTemp = 0.f;
		float flowReportZ = 30.f;
		getIJKCellIndexFromXYZ( iTemp, jTemp, flowReportK, xTemp, yTemp, flowReportZ, grids[gridLevelCount-1].Info);
		getFlowReport( flowReportK, grids[gridLevelCount-1], iStart, jStart, iEnd, jEnd, uzAvg, massFlow, pAvg );
		
		float uzAvgInlet, massFlowInlet, pAvgInlet;
		flowReportZ = -100.f;
		getIJKCellIndexFromXYZ( iTemp, jTemp, flowReportK, xTemp, yTemp, flowReportZ, grids[gridLevelCount-1].Info);
		getFlowReport( flowReportK, grids[gridLevelCount-1], 0, 0, grids[gridLevelCount-1].Info.cellCountX, grids[gridLevelCount-1].Info.cellCountY, uzAvgInlet, massFlowInlet, pAvgInlet );

		float lakePower = 0.5f * massFlow * uzAvgInlet * uzAvgInlet;
		float intakePower = 0.5f * massFlow * uzAvg * uzAvg + massFlow * (pAvg - pAvgInlet) / rhoNominalPhys;
		float eta = intakePower / lakePower;
		
		historyMassFlow[iteration] = massFlow;
		historyEta[iteration] = eta;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
	
			for ( int level = gridLevelCount-2; level >= 0; level-- )
			{
				fillCoarseGridFromFine( grids[level], grids[level+1] );
			}
			
			for ( int level = 0; level < gridLevelCount; level++ )
			{
				int iCut, jCut, kCut;
				const float xCut = 0.f; 
				const float yCut = -30.f;
				const float zCut = 30.f;
				getIJKCellIndexFromXYZ( iCut, jCut, kCut, xCut, yCut, zCut, grids[level].Info);
				exportSectionCutPlotXY( grids[level], kCut, iteration + level );
				system("python3 ../plotter/plotter.py");
				exportSectionCutPlotZY( grids[level], iCut, iteration + level + 10 );
				system("python3 ../plotter/plotter.py");
				exportSectionCutPlotZX( grids[level], jCut, iteration + level + 20 );
				system("python3 ../plotter/plotter.py");
			}
			
			//export3DPlot( grids[gridLevelCount-1], iteration + 30 );
			//system("python3 plotter3D.py");
			
			exportHistoryData( historyMassFlow, historyEta, iteration, caseID );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
