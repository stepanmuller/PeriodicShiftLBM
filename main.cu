// input physics
constexpr float res = 0.5f; 											// mm
constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number

constexpr float nuPhys = 1.4607e-5;										// m2/s air
constexpr float rhoNominalPhys = 1.225f;								// kg/m3 air
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES

// box dimensions
constexpr size_t boxStartJ = 440;
constexpr size_t boxEndJ = 560;
constexpr size_t boxStartK = 340;
constexpr size_t boxEndK = 460;

// calculated from input
constexpr float uzInletPhys = 1.f; 									// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * res; 			// s
const float soundspeedPhys = (1.f / sqrt(3.f)) * res / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / (res * res);				// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau

constexpr int iterationCount = 100;

#include "includesTypes.h"

int main(int argc, char **argv)
{
	STLArbeiterStructCPU STLArbeiterCPU;
	
	readSTL( STLArbeiterCPU );
	
	std::cout << "Sizing domain around the STL" << std::endl;
	CellCountStruct cellCount;
	cellCount.res = res;
	cellCount.nx = static_cast<size_t>((STLArbeiterCPU.xmax - STLArbeiterCPU.xmin) / cellCount.res) + 1;
	cellCount.ny = static_cast<size_t>((STLArbeiterCPU.ymax - STLArbeiterCPU.ymin) / cellCount.res) + 1;
	cellCount.nz = static_cast<size_t>((STLArbeiterCPU.zmax - STLArbeiterCPU.zmin) / cellCount.res) + 1;
	
	cellCount.ox = STLArbeiterCPU.xmin + ( 0.5f * ( (STLArbeiterCPU.xmax - STLArbeiterCPU.xmin) - cellCount.res * (cellCount.nx-1) ) );
	cellCount.oy = STLArbeiterCPU.ymin + ( 0.5f * ( (STLArbeiterCPU.ymax - STLArbeiterCPU.ymin) - cellCount.res * (cellCount.ny-1) ) );
	cellCount.oz = STLArbeiterCPU.zmin + ( 0.5f * ( (STLArbeiterCPU.zmax - STLArbeiterCPU.zmin) - cellCount.res * (cellCount.nz-1) ) );
	
	cellCount.ny = cellCount.ny + 1; // adding one more "wall" layer on top
	
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
	
	#ifdef __CUDACC__
	std::cout << "Starting simulation" << std::endl;
	TNL::Timer timer;
	TNL::Timer lapTimer;
	timer.start();
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
			lapTimer.start();
		}
	}
	timer.stop();
	#endif
	auto totalTime = timer.getRealTime();
	std::cout << "This took " << totalTime << " s" << std::endl;
	float glups = (cellCount.n * iterationCount) / totalTime / 1000000000;
	std::cout << "Total average GLUPS: " << glups << std::endl;
	
	std::cout << "Saving result" << std::endl;
	
	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	// Write metadata first so Python knows the dimensions
	int dims[2] = {(int)cellCount.ny, (int)cellCount.nz};
	fwrite(dims, sizeof(int), 2, fp);
	
	DistributionStructCPU FCPU;
	FCPU.shifter = IndexArrayType( 27, 0 );
	FCPU.fArray.setSizes( 27, cellCount.n );
	FCPU.shifter = F.shifter;
	FCPU.fArray = F.fArray;
	
	for (size_t j = 0; j < cellCount.ny; j++)
	{
		for (size_t k = 0; k < cellCount.nz; k++)
		{
			size_t i = cellCount.nx / 2; 
			size_t cell = convertIndex(i, j, k, cellCount);
			size_t shiftedIndex[27];
			for (size_t i = 0; i < 27; i++) 
			{
				const size_t shift = FCPU.shifter[i];
				shiftedIndex[i] = cell + shift;
				if (shiftedIndex[i] >= cellCount.n) { shiftedIndex[i] -= cellCount.n; }
			}
			float f[27];
			float rho, ux, uy, uz;
			for (size_t i = 0; i < 27; i++)	f[i] = FCPU.fArray.getElement(i, shiftedIndex[i]);
			getRhoUxUyUz(rho, ux, uy, uz, f);
			float uMag = sqrt(uy * uy + uz * uz);
			fwrite(&uMag, sizeof(float), 1, fp);
		}
	}
	fclose(fp);
	system("python3 visualizer.py");
	return EXIT_SUCCESS;
}
