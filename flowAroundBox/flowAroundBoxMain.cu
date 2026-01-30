constexpr float uzInlet = 0.05f; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.0f; 							// set to zero to turn off LES
constexpr float nu = 1e-6;												// LBM nu
constexpr float tau = 3.f * nu + 0.5f;									// LBM tau

constexpr int cellCountX = 50;
constexpr int cellCountY = 1000;
constexpr int cellCountZ = 2000;

constexpr int boxStartY = (int)(cellCountY * 0.4f);
constexpr int boxStartZ = (int)(cellCountY * 0.3f);
constexpr int boxEndY = (int)(cellCountY * 0.6f);
constexpr int boxEndZ = (int)(cellCountY * 0.5f);

constexpr int iterationCount = 10000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"
#include "../fillDefaultEquilibrium.h"
#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

void applyLocalCellUpdate( FStruct &F, InfoStruct &Info )
{
	auto fArrayView  = F.fArray.getView();
	auto shifterView  = F.shifter.getConstView();

	auto cellLambda = [=] __cuda_callable__ (size_t cell) mutable
	{
		size_t cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if (shiftedIndex[i] >= cellCount) shiftedIndex[i] -= cellCount;
		}
		bool fluidMarker = 0;
		bool bouncebackMarker = 0;
		bool givenRhoMarker = 0;
		bool givenUxUyUzMarker = 0;
		size_t iCell, jCell, kCell;
		getCellIJK( cell, iCell, jCell, kCell, Info );
		
		if ( jCell >= boxStartY && jCell < boxEndY && kCell >= boxStartZ && kCell < boxEndZ ) bouncebackMarker = 1;
		else if ( iCell == 0 || iCell == Info.cellCountX-1 || jCell == 0 || jCell == Info.cellCountY-1 ) bouncebackMarker = 1;
		else if ( kCell == 0 ) givenUxUyUzMarker = 1;
		else if ( kCell == Info.cellCountZ-1 ) givenRhoMarker = 1;
		else fluidMarker = 1;
		
		float f[27];
		float rho, ux, uy, uz;
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView(i, shiftedIndex[i]);
		
		if ( bouncebackMarker )
		{
			applyBounceback(f);
		}
		else 
		{
			if ( fluidMarker )
			{
				getRhoUxUyUz( rho, ux, uy, uz, f );
			}
			else
			{
				int outerNormalX, outerNormalY, outerNormalZ;
				getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info );
				rho = 1.f;
				ux = 0.f;
				uy = 0.f;
				uz = uzInlet;
				if ( givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && givenUxUyUzMarker )
				{
					restoreRho( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				else if ( !givenRhoMarker && !givenUxUyUzMarker )
				{
					restoreRhoUxUyUz( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
				}
				applyMBBC( outerNormalX, outerNormalY, outerNormalZ, rho, ux, uy, uz, f );
			}
			applyCollision( rho, ux, uy, uz, f );
		}
		for ( size_t i = 0; i < 27; i++ ) fArrayView( i, shiftedIndex[i] ) = f[i];
	};
	size_t start = 0;
	size_t end = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
}

void exportSectionCutPlot( FStruct &F, InfoStruct &Info, const size_t &iCell, const int &plotNumber )
{
	std::cout << "Exporting section cut plot " << plotNumber << std::endl;
	auto fArrayView  = F.fArray.getConstView();
	auto shifterView  = F.shifter.getConstView();
	
	SectionCutStruct SectionCut;
	SectionCut.rhoArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uxArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uyArray.setSizes( Info.cellCountY, Info.cellCountZ );
	SectionCut.uzArray.setSizes( Info.cellCountY, Info.cellCountZ );
		
	auto rhoArrayView = SectionCut.rhoArray.getView();
	auto uxArrayView = SectionCut.uxArray.getView();
	auto uyArrayView = SectionCut.uyArray.getView();
	auto uzArrayView = SectionCut.uzArray.getView();

	auto cellLambda = [=] __cuda_callable__ (const LongIntPairType &doubleIndex) mutable
	{
		const size_t jCell = doubleIndex[0];
		const size_t kCell = doubleIndex[1];
		const size_t cell = kCell * (Info.cellCountX * Info.cellCountY) + jCell * Info.cellCountX + iCell;
		const size_t cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
		size_t shiftedIndex[27];
		for (size_t i = 0; i < 27; i++) 
		{
			const size_t shift = shifterView[i];
			shiftedIndex[i] = cell + shift;
			if ( shiftedIndex[i] >= cellCount ) shiftedIndex[i] -= cellCount;
		}		
		float f[27];
		for (size_t i = 0; i < 27; i++)	f[i] = fArrayView( i, shiftedIndex[i] );
				
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		rhoArrayView( jCell, kCell ) = rho;
		uxArrayView( jCell, kCell ) = ux;
		uyArrayView( jCell, kCell ) = uy;
		uzArrayView( jCell, kCell ) = uz;
	};
	LongIntPairType start{ 0, 0 };
	LongIntPairType end{ Info.cellCountY, Info.cellCountZ };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(start, end, cellLambda );
	
	SectionCutStructCPU SectionCutCPU;
	SectionCutCPU.rhoArray = SectionCut.rhoArray;
	SectionCutCPU.uxArray = SectionCut.uxArray;
	SectionCutCPU.uyArray = SectionCut.uyArray;
	SectionCutCPU.uzArray = SectionCut.uzArray;
	
	// Use /dev/shm/ for a pure RAM-based "file" on Linux
	FILE* fp = fopen("/dev/shm/sim_data.bin", "wb");
	// Write metadata first so Python knows the dimensions
	int header[4] = {plotNumber, Info.cellCountY, Info.cellCountZ, 4};
	fwrite(header, sizeof(int), 4, fp);
	
	for (int j = 0; j < Info.cellCountY; j++)
	{
		for (int k = 0; k < Info.cellCountZ; k++)
		{
			float rho = SectionCutCPU.rhoArray.getElement(j, k);
			float ux = SectionCutCPU.uxArray.getElement(j, k);
			float uy = SectionCutCPU.uyArray.getElement(j, k);
			float uz = SectionCutCPU.uzArray.getElement(j, k);
			float data[4] = {rho, ux, uy, uz};
			fwrite(data, sizeof(float), 4, fp);
		}
	}
	fclose(fp);
	system("python3 flowAroundBoxPlotter.py");
}

int main(int argc, char **argv)
{
	InfoStruct Info;
	Info.cellCountX = cellCountX;
	Info.cellCountY = cellCountY;
	Info.cellCountZ = cellCountZ;
	Info.iterationsFinished = 0;
	
	FStruct F;
	FloatArray2DType fArray;
	F.fArray.setSizes( 27, Info.cellCountX * Info.cellCountY * Info.cellCountZ );	
	F.shifter = LongIntArrayType( 27, 0 );
	
	fillDefaultEquilibrium( F, Info);
	
	std::cout << "Starting simulation" << std::endl;
	
	const int iCut = Info.cellCountX / 2;
	int plotNumber = 0;
	
	const int iterationChunk = 100;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<iterationCount; iteration++)
	{
		applyStreaming( F, Info );
		applyLocalCellUpdate( F, Info );
		Info.iterationsFinished++;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			size_t cellCount = Info.cellCountX * Info.cellCountY * Info.cellCountZ;
			float glups = (cellCount * (size_t)iterationChunk) / lapTime / 1000000000;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportSectionCutPlot( F, Info, iCut, plotNumber );
			plotNumber++;
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
