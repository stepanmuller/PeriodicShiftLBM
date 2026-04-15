constexpr float reynoldsNumber = 100.f;

constexpr float sphereDiameterPhys = 2000.f;											// mm
constexpr float resGlobal = 80.f; 														// mm
constexpr float uxInlet = 0.07f; 														// also works as nominal LBM Mach number
constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float sphereRadiusPhys = 0.5 * sphereDiameterPhys;							// mm
constexpr float uxInletPhys = uxInlet; 													// m/s, physical velocity set to same as LBM velocity
constexpr float nuPhys = uxInletPhys * (sphereDiameterPhys / 1000.f) / reynoldsNumber;	// m2/s
constexpr float rhoNominalPhys = 1.225f;												// kg/m3 air
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float domainSizePhys = 11.f * sphereDiameterPhys;								// mm
constexpr float sphereXPhys = 0.f;														// mm
constexpr float sphereYPhys = 0.f;														// mm
constexpr float sphereZPhys = 0.f;														// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys / resGlobal));
const int cellCountY = cellCountX;
const int cellCountZ = cellCountX;

constexpr int iterationCount = 400000;
constexpr int iterationChunk = 20000;
constexpr int gridLevelCount = 6;
constexpr int wallRefinementSpan = 4;

#include "../../include/types.h"

#include "../../include/cellFunctions.h"
#include "../../include/applyStreaming.h"
#include "../../include/applyCollision.h"

#include "../../include/boundaryConditions/applyBounceback.h"
#include "../../include/boundaryConditions/applyMirror.h"
#include "../../include/boundaryConditions/restoreRho.h"
#include "../../include/boundaryConditions/restoreUxUyUz.h"
#include "../../include/boundaryConditions/restoreRhoUxUyUz.h"
#include "../../include/boundaryConditions/applyMBBC.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
   	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float zPhys = kCell * Info.res + Info.oz;
	const float r2 = (xPhys-sphereXPhys) * (xPhys - sphereXPhys) + (yPhys - sphereYPhys) * (yPhys - sphereYPhys) + (zPhys - sphereZPhys) * (zPhys - sphereZPhys);
	const float axialDistance = (yPhys - sphereYPhys) * (yPhys - sphereYPhys) + (zPhys - sphereZPhys) * (zPhys - sphereZPhys);
	
	// Enlarge the refinement area
	if ( Info.gridID == 0 )
	{
		const float widthHalf = 1.5f * sphereDiameterPhys;
		const float length = 4.f * sphereDiameterPhys;
		if (xPhys > 0 && xPhys < length	&& axialDistance < widthHalf * widthHalf ) Marker.refinement = 1;
		else if ( r2 < widthHalf * widthHalf ) Marker.refinement = 1;
	}
	else if ( Info.gridID == 1 )
	{
		const float widthHalf = 1.f * sphereDiameterPhys;
		const float length = 2.5f * sphereDiameterPhys;
		if (xPhys > 0 && xPhys < length	&& axialDistance < widthHalf * widthHalf ) Marker.refinement = 1;
		else if ( r2 < widthHalf * widthHalf ) Marker.refinement = 1;
	}
	else if ( Info.gridID == 2 )
	{
		const float widthHalf = 0.7f * sphereDiameterPhys;
		const float length = 1.0f * sphereDiameterPhys;
		if (xPhys > 0 && xPhys < length	&& axialDistance < widthHalf * widthHalf ) Marker.refinement = 1;
		else if ( r2 < widthHalf * widthHalf ) Marker.refinement = 1;
	}
	else if ( Info.gridID == 3 )
	{
		const float widthHalf = 0.6f * sphereDiameterPhys;
		const float length = 0.6f * sphereDiameterPhys;
		if (xPhys > 0 && xPhys < length	&& axialDistance < widthHalf * widthHalf ) Marker.refinement = 1;
		else if ( r2 < widthHalf * widthHalf ) Marker.refinement = 1;
	}
	
	if ( Info.gridID != 0 ) // if not zero, we are on a finer grid
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) Marker.bounceback = 1;
		else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.ghost = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.ghost = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.ghost = 1;
		else Marker.fluid = 1;
	}	
	else // we are on the coarse grid
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) Marker.bounceback = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.givenUxUyUz = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.givenUxUyUz = 1;
		else if ( iCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( iCell == Info.cellCountX-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}
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

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	//if ( Marker.bounceback ) ux = uxInlet;
	//else ux = 0.f;
	ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../../include/applyLocalCellUpdate.h"
#include "../../include/plotter/exportSectionCutPlot.h"
#include "../../include/fillEquilibrium.h"
#include "../../include/gridRefinementFunctions.h"

#include "../../include/DIADFunctions.h"

float getSphereDrag( DIADGridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getConstView();
	InfoStruct Info = Grid.Info;
	
	bool esotwistFlipper = Grid.esotwistFlipper;
	auto iNbrView = Grid.EsotwistNbrArray.iNbrArray.getConstView();
	auto jNbrView = Grid.EsotwistNbrArray.jNbrArray.getConstView();
	auto kNbrView = Grid.EsotwistNbrArray.kNbrArray.getConstView();
	auto ijNbrView = Grid.EsotwistNbrArray.ijNbrArray.getConstView();
	auto ikNbrView = Grid.EsotwistNbrArray.ikNbrArray.getConstView();
	auto jkNbrView = Grid.EsotwistNbrArray.jkNbrArray.getConstView();
	auto ijkNbrView = Grid.EsotwistNbrArray.ijkNbrArray.getConstView();
	
	auto iView = Grid.IJK.iArray.getConstView();
	auto jView = Grid.IJK.jArray.getConstView();
	auto kView = Grid.IJK.kArray.getConstView();
	
	bool useBouncebackArray = ( Grid.bouncebackMarkerArray.getSize() > 0 );
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	
	auto fetch = [ = ] __cuda_callable__( const int cell )
	{		
		const int iCell = iView( cell );
		const int jCell = jView( cell );
		const int kCell = kView( cell );
		
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( !Marker.bounceback ) return 0.f;
		
		DIADEsotwistNbrStruct Nbr;
		Nbr.i = iNbrView( cell );
		Nbr.j = jNbrView( cell );
		Nbr.k = kNbrView( cell );
		Nbr.ij = ijNbrView( cell );
		Nbr.ik = ikNbrView( cell );
		Nbr.jk = jkNbrView( cell );
		Nbr.ijk = ijkNbrView( cell );
		
		float f[27];
		int cellReadIndex[27];
		int fReadIndex[27];
		getEsotwistWriteIndex( cell, cellReadIndex, fReadIndex, Nbr, esotwistFlipper, Info ); 
		// Using the write index because plotting is happening after collision and we want to be consistent with that last write
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fReadIndex[direction], cellReadIndex[direction]);
		
		float rho, ux, uy, uz;
		getRhoUxUyUz( rho, ux, uy, uz, f );
		applyBounceback( f );
		float rhoPrev, uxPrev, uyPrev, uzPrev;
		getRhoUxUyUz( rhoPrev, uxPrev, uyPrev, uzPrev, f );
		const float gx = rho * (ux -  uxPrev);
		return gx;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float gxSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, Info.cellCount, fetch, reduction, 0.f );
	float gy, gz = 0;
	convertToPhysicalForce( gxSum, gy, gz, Info );
	return gxSum;
}

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

int main(int argc, char **argv)
{
	std::vector<STLStruct> STLs = { };
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = cellCountX;
	grids[0].Info.cellCountY = cellCountY;
	grids[0].Info.cellCountZ = cellCountZ;
	grids[0].Info.ox = - 2.f * sphereDiameterPhys;
	grids[0].Info.oy = - 5.5f * sphereDiameterPhys;
	grids[0].Info.oz = - 5.5f * sphereDiameterPhys;
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
	
	for ( int level = 0; level < gridLevelCount; level++ )
	{	
		grids[level].fArray.setSizes( 27, grids[level].Info.cellCount );
		fillEquilibriumFromFunction( grids[level] );
	}
	
	std::vector<float> historyVector( iterationCount, 0.f );
	
	std::cout << "Finest cells per sphere diameter: " << (int)std::round(sphereDiameterPhys / grids[gridLevelCount-1].Info.res) << std::endl;
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		const float drag = getSphereDrag( grids[gridLevelCount-1] );
		const float dragCoefficient = - (8 * drag) / (rhoNominalPhys * uxInletPhys * uxInletPhys * 3.14159f * (sphereDiameterPhys / 1000.f) * (sphereDiameterPhys / 1000.f));
		
		historyVector[iteration] = dragCoefficient;
		
		if (iteration%iterationChunk == 0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			const int kCut = grids[gridLevelCount-1].Info.cellCountZ / 2;
			exportSectionCutPlotXY( grids, kCut, iteration );
			system("python3 ../plotterGridID.py");
			
			exportHistoryData( historyVector, iteration );
			system("python3 ../historyPlotter.py &"); 
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
