// input physics
constexpr float res = 4.8f; 											// mm
constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number

constexpr float nuPhys = 1.5e-5;										// m2/s air
constexpr float rhoNominalPhys = 1.225f;								// kg/m3 water
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES

// calculated from input
constexpr float uzInletPhys = 90.f; 								// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * (res/1000); 		// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (res/1000) / dtPhys; 	// m/s
constexpr float nu = (dtPhys * nuPhys) / ((res/1000) * (res/1000));	// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau

constexpr int iterationCount = 300000;

#include "../includesTypes.h"
#include "exportSectionCutPlot.h"

std::string STLPath = "NACA0012_10deg_STL.STL";

void applyMarkers( MarkerStruct& Marker, CellCountStruct &cellCount )
{
	auto fluidArrayView = Marker.fluidArray.getView();
	auto bouncebackArrayView = Marker.bouncebackArray.getView();
	auto givenRhoArrayView = Marker.givenRhoArray.getView();
	auto givenUxUyUzArrayView = Marker.givenUxUyUzArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 3, int >& tripleIndex) mutable
	{
		const size_t i = tripleIndex.x();
		const size_t j = tripleIndex.y();
		const size_t k = tripleIndex.z();
		size_t cell = convertIndex(i, j, k, cellCount);
		if (bouncebackArrayView[cell] == 0)
		{
			if (k==0)  
			{
				givenUxUyUzArrayView[cell] = 1;
			}
			else if ( k==cellCount.nz-1 || i==0 || i==cellCount.nx-1 || j==0 || j==cellCount.ny-1 )
			{
				givenRhoArrayView[cell] = 1;
			}
			else
			{
				fluidArrayView[cell] = 1;
			}
		}
	};
	TNL::Containers::StaticArray< 3, size_t > start{ 0, 0, 0 };
	TNL::Containers::StaticArray< 3, size_t > end{ cellCount.nx, cellCount.ny, cellCount.nz };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

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
