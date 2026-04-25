void fillOutletCellArray( DIADGridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	
	auto iView = Grid.IJK.iArray.getConstView();
	auto jView = Grid.IJK.jArray.getConstView();
	auto kView = Grid.IJK.kArray.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	
	BoolArrayType skipCellMarkerArray;
	skipCellMarkerArray.setSize( Info.cellCount );
	skipCellMarkerArray.setValue( 0 );
	auto skipCellMarkerView = skipCellMarkerArray.getView();
	
	auto cellLambda = [=] __cuda_callable__ ( const int cell ) mutable
	{
		const int iCell = iView( cell );
		const int jCell = jView( cell );
		const int kCell = kView( cell );
		int outerNormalX, outerNormalY, outerNormalZ;
		getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info ); 
		if ( abs(outerNormalX) + abs(outerNormalY) + abs(outerNormalZ) > 1 ) skipCellMarkerView( cell ) = 1;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, Info.cellCount, cellLambda );
	
	const int interfaceCellCountCF = Grid.coarseToFineReadArray.getSize();
	auto coarseToFineReadView = Grid.coarseToFineReadArray.getConstView();
	auto cellLambdaCF = [=] __cuda_callable__ ( const int index ) mutable
	{
		skipCellMarkerView( coarseToFineReadView( index ) ) = 1;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, interfaceCellCountCF, cellLambdaCF );
	
	const int interfaceCellCountFC = Grid.fineToCoarseWriteArray.getSize();
	auto fineToCoarseWriteView = Grid.fineToCoarseWriteArray.getConstView();
	auto cellLambdaFC = [=] __cuda_callable__ ( const int index ) mutable
	{
		skipCellMarkerView( fineToCoarseWriteView( index ) ) = 1;
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(0, interfaceCellCountFC, cellLambdaFC );
	
	auto fetch = [ = ] __cuda_callable__( const int cell )
	{
		const int iCell = iView( cell );
		const int jCell = jView( cell );
		const int kCell = kView( cell );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.givenRho && !skipCellMarkerView( cell ) ) return 1;
		else return 0;
	};
	auto reduction = [] __cuda_callable__( const int& a, const int& b )
	{
		return a + b;
	};	
	int outletCellCount = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, Grid.Info.cellCount, fetch, reduction, 0 );
	
	IntArrayTypeCPU iArrayCPU;
	IntArrayTypeCPU jArrayCPU;
	IntArrayTypeCPU kArrayCPU;
	BoolArrayTypeCPU bouncebackMarkerArrayCPU;
	BoolArrayTypeCPU skipCellMarkerArrayCPU;
	iArrayCPU = Grid.IJK.iArray;
	jArrayCPU = Grid.IJK.jArray;
	kArrayCPU = Grid.IJK.kArray;
	bouncebackMarkerArrayCPU = Grid.bouncebackMarkerArray;
	skipCellMarkerArrayCPU = skipCellMarkerArray;
	IntArrayTypeCPU outletCellArrayCPU;
	outletCellArrayCPU.setSize( outletCellCount );
	int counter = 0;
	for ( int cell = 0; cell < Grid.Info.cellCount; cell++ ) 
	{
		const int iCell = iArrayCPU[ cell ];
		const int jCell = jArrayCPU[ cell ];
		const int kCell = kArrayCPU[ cell ];
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayCPU[ cell ];
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.givenRho && !skipCellMarkerArrayCPU[ cell ] ) outletCellArrayCPU[ counter++ ] = cell;
	}
	Grid.outletCellArray = outletCellArrayCPU;
}

void getAverageOutletUSingle( DIADGridStruct &Grid, float &gridAreamm2, float &gridU )
{	
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
	
	const int outletCellCount = Grid.outletCellArray.getSize();
	gridAreamm2 = outletCellCount * Info.res * Info.res;
	auto outletCellView = Grid.outletCellArray.getConstView();
	
	auto fetch = [ = ] __cuda_callable__( const int index )
	{
		const int cell = outletCellView( index );
		
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
		
		float f[27];
		int cellWriteIndex[27];
		int fWriteIndex[27];
		getEsotwistWriteIndex( cell, cellWriteIndex, fWriteIndex, Nbr, esotwistFlipper, Info );
		for ( int direction = 0; direction < 27; direction++ )	f[direction] = fArrayView(fWriteIndex[direction], cellWriteIndex[direction]);
		
		float rho, ux, uy, uz;
		getRhoUxUyUz( rho, ux, uy, uz, f );
		
		int outerNormalX, outerNormalY, outerNormalZ;
		getOuterNormal( iCell, jCell, kCell, outerNormalX, outerNormalY, outerNormalZ, Info ); 
		
		if ( outerNormalX != 0 ) return (ux * (float)outerNormalX);
		else if ( outerNormalY != 0 ) return (uy * (float)outerNormalY);
		else return (uz * (float)outerNormalZ);
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float uSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( 0, outletCellCount, fetch, reduction, 0.f );
	gridU = uSum / outletCellCount;
}

float getAverageOutletU( std::vector<DIADGridStruct>& grids )
{
	float areamm2 = 0.f;
	float uTimesAreamm2 = 0.f;
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		float gridAreamm2 = 0.f;
		float gridU = 0.f;
		getAverageOutletUSingle( grids[level], gridAreamm2, gridU );
		areamm2 += gridAreamm2;
		uTimesAreamm2 += gridU * gridAreamm2;
	}
	float averageOutletU = uTimesAreamm2 / areamm2;
	return averageOutletU;
}

// REGULATING OUTLET RHO TO STRONGLY DAMPEN ACUSTIC WAVES
void regulateOutletPressure( std::vector<DIADGridStruct>& grids, float &uOutletMovingAvg )
{
	float uzAvg = getAverageOutletU( grids );
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		grids[level].Info.pRegulator = ( uzAvg - uOutletMovingAvg ) * 1.732050808f;
		grids[level].Info.iRegulator -= 0.02f * ( grids[0].Info.pRegulator + grids[0].Info.iRegulator );
	}
	uOutletMovingAvg = 0.999f * uOutletMovingAvg + 0.001f * uzAvg;
}
