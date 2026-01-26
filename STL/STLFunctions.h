void readSTL( STLArbeiterStructCPU& STLArbeiterCPU )
{
	std::ifstream file("STL/M35IntakeSTL.STL", std::ios::binary);
	if (!file.is_open()) throw std::runtime_error("Failed to open STL file");
	
	 // Skip header
    char header[80];
    file.read(header, 80);

    uint32_t triangleCount32;
	file.read(reinterpret_cast<char*>(&triangleCount32), sizeof(uint32_t));
	size_t triangleCount = static_cast<size_t>(triangleCount32);
	STLArbeiterCPU.triangleCount = triangleCount;
    
    std::cout<<"	triangleCount: "<< triangleCount << std::endl;

    STLArbeiterCPU.axArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.ayArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.azArray = FloatArrayTypeCPU(triangleCount);

    STLArbeiterCPU.bxArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.byArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.bzArray = FloatArrayTypeCPU(triangleCount);

    STLArbeiterCPU.cxArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.cyArray = FloatArrayTypeCPU(triangleCount);
    STLArbeiterCPU.czArray = FloatArrayTypeCPU(triangleCount);
    
    // Initialize minmax
    STLArbeiterCPU.xmin = std::numeric_limits<float>::max();
	STLArbeiterCPU.ymin = std::numeric_limits<float>::max();
	STLArbeiterCPU.zmin = std::numeric_limits<float>::max();

	STLArbeiterCPU.xmax = std::numeric_limits<float>::lowest();
	STLArbeiterCPU.ymax = std::numeric_limits<float>::lowest();
	STLArbeiterCPU.zmax = std::numeric_limits<float>::lowest();

    for (size_t triangle = 0; triangle < triangleCount; triangle++)
    {
        float nx, ny, nz;
        float ax, ay, az;
        float bx, by, bz;
        float cx, cy, cz;
        uint16_t attr;

        file.read(reinterpret_cast<char*>(&nx), 4);
        file.read(reinterpret_cast<char*>(&ny), 4);
        file.read(reinterpret_cast<char*>(&nz), 4);

        file.read(reinterpret_cast<char*>(&ax), 4);
        file.read(reinterpret_cast<char*>(&ay), 4);
        file.read(reinterpret_cast<char*>(&az), 4);

        file.read(reinterpret_cast<char*>(&bx), 4);
        file.read(reinterpret_cast<char*>(&by), 4);
        file.read(reinterpret_cast<char*>(&bz), 4);

        file.read(reinterpret_cast<char*>(&cx), 4);
        file.read(reinterpret_cast<char*>(&cy), 4);
        file.read(reinterpret_cast<char*>(&cz), 4);

        file.read(reinterpret_cast<char*>(&attr), 2);

        STLArbeiterCPU.axArray[triangle] = ax;
        STLArbeiterCPU.ayArray[triangle] = ay;
        STLArbeiterCPU.azArray[triangle] = az;

        STLArbeiterCPU.bxArray[triangle] = bx;
        STLArbeiterCPU.byArray[triangle] = by;
        STLArbeiterCPU.bzArray[triangle] = bz;

        STLArbeiterCPU.cxArray[triangle] = cx;
        STLArbeiterCPU.cyArray[triangle] = cy;
        STLArbeiterCPU.czArray[triangle] = cz;
        
         // Update bounding box (vertices only)
        STLArbeiterCPU.xmin = std::min(STLArbeiterCPU.xmin, std::min({ax, bx, cx}));
        STLArbeiterCPU.ymin = std::min(STLArbeiterCPU.ymin, std::min({ay, by, cy}));
        STLArbeiterCPU.zmin = std::min(STLArbeiterCPU.zmin, std::min({az, bz, cz}));

        STLArbeiterCPU.xmax = std::max(STLArbeiterCPU.xmax, std::max({ax, bx, cx}));
        STLArbeiterCPU.ymax = std::max(STLArbeiterCPU.ymax, std::max({ay, by, cy}));
        STLArbeiterCPU.zmax = std::max(STLArbeiterCPU.zmax, std::max({az, bz, cz}));
    }
    std::cout << "	xmin xmax: " << STLArbeiterCPU.xmin << " " << STLArbeiterCPU.xmax << "\n";
    std::cout << "	ymin ymax: " << STLArbeiterCPU.ymin << " " << STLArbeiterCPU.ymax << "\n";
    std::cout << "	zmin zmax: " << STLArbeiterCPU.zmin << " " << STLArbeiterCPU.zmax << "\n";
}

__host__ __device__ void rayHitDetector(const long long &ak, const long long &al, 
										const long long &bk, const long long &bl,
										const long long &ck, const long long &cl,
										bool &rayHit)
{
    // Calculate Edge Functions
    long long w0 = bk * cl - bl * ck;
    long long w1 = ck * al - cl * ak;
    long long w2 = ak * bl - al * bk;

    long long area = w0 + w1 + w2;
    if (area == 0) {
        rayHit = false; 
        return;
    }

    bool flipped = (area < 0);
    if (flipped) { w0 = -w0; w1 = -w1; w2 = -w2; }
    
    // EDGE 0 CHECK
    if (w0 < 0) { rayHit = false; return; }
    if (w0 == 0) {
        long long dx = flipped ? (ak - bk) : (bk - ak);
        long long dy = flipped ? (al - bl) : (bl - al);
        if (!((dy < 0) || (dy == 0 && dx > 0))) { rayHit = false; return; }
    }

    // EDGE 1 CHECK
    if (w1 < 0) { rayHit = false; return; }
    if (w1 == 0) {
        long long dx = flipped ? (bk - ck) : (ck - bk);
        long long dy = flipped ? (bl - cl) : (cl - bl);
        if (!((dy < 0) || (dy == 0 && dx > 0))) { rayHit = false; return; }
    }

    // EDGE 2 CHECK
    if (w2 < 0) { rayHit = false; return; }
    if (w2 == 0) {
        long long dx = flipped ? (ck - ak) : (ak - ck);
        long long dy = flipped ? (cl - al) : (al - cl);
        if (!((dy < 0) || (dy == 0 && dx > 0))) { rayHit = false; return; }
    }

    // If we reached this point, all checks passed!
    rayHit = true;
}

void applyMarkersFromSTL( MarkerStruct &Marker, STLArbeiterStruct &STLArbeiter, CellCountStruct &cellCount )
{
	auto fluidMarkerArrayView = Marker.fluidArray.getView();
	auto bouncebackMarkerArrayView = Marker.bouncebackArray.getView();
	auto givenRhoMarkerArrayView = Marker.givenRhoArray.getView();
	auto givenUxUyUzMarkerArrayView = Marker.givenUxUyUzArray.getView();

	auto axArrayView = STLArbeiter.axArray.getConstView();
	auto ayArrayView = STLArbeiter.ayArray.getConstView();
	auto azArrayView = STLArbeiter.azArray.getConstView();
	auto bxArrayView = STLArbeiter.bxArray.getConstView();
	auto byArrayView = STLArbeiter.byArray.getConstView();
	auto bzArrayView = STLArbeiter.bzArray.getConstView();
	auto cxArrayView = STLArbeiter.cxArray.getConstView();
	auto cyArrayView = STLArbeiter.cyArray.getConstView();
	auto czArrayView = STLArbeiter.czArray.getConstView();
	
	CounterArray2DType intersectionCounterArray;
	intersectionCounterArray.setSizes(cellCount.nx, cellCount.ny);
	intersectionCounterArray.setValue(0);
	auto intersectionCounterArrayView = intersectionCounterArray.getView();

    auto counterLambda = [ = ] __cuda_callable__( size_t triangleIndex ) mutable
    {
		const float a0x = axArrayView[ triangleIndex ];
		const float a0y = ayArrayView[ triangleIndex ];
		const float b0x = bxArrayView[ triangleIndex ];
		const float b0y = byArrayView[ triangleIndex ];
		const float c0x = cxArrayView[ triangleIndex ];
		const float c0y = cyArrayView[ triangleIndex ];
		// transform STL floats to integer grid that is 100x finer than the LBM grid to prevent float errors
		// transform into coordinate system of the LBM grid
		// make the STL coords odd, rays will be even, this prevents hitting a vortex
		const float scale = 50.0f / cellCount.res;
		const long long ok = round( cellCount.ox * scale );
		const long long ol = round( cellCount.oy * scale );
		const long long a0k = (long long)(round( a0x * scale ) - ok) * 2 + 1;
		const long long a0l = (long long)(round( a0y * scale ) - ol) * 2 + 1;
		const long long b0k = (long long)(round( b0x * scale ) - ok) * 2 + 1;
		const long long b0l = (long long)(round( b0y * scale ) - ol) * 2 + 1;
		const long long c0k = (long long)(round( c0x * scale ) - ok) * 2 + 1;
		const long long c0l = (long long)(round( c0y * scale ) - ol) * 2 + 1;
		
		const long long kmin = std::min({ a0k, b0k, c0k });
		const long long kmax = std::max({ a0k, b0k, c0k });
		const long long lmin = std::min({ a0l, b0l, c0l });
		const long long lmax = std::max({ a0l, b0l, c0l });
		
		const size_t imin = min( (size_t)max(0LL, (kmin / 100) ), cellCount.nx-1 );
		const size_t imax = min( (size_t)max(0LL, (kmax / 100) + 1), cellCount.nx-1 );
		const size_t jmin = min( (size_t)max(0LL, (lmin / 100) ), cellCount.ny-1 );
		const size_t jmax = min( (size_t)max(0LL, (lmax / 100) + 1), cellCount.ny-1 );
		
		for (size_t j = jmin; j <= jmax; j++) 
		{
			for (size_t i = imin; i <= imax; i++) 
			{
				const long long rayK = i * 100;
				const long long rayL = j * 100;
				// finally, transform the triangle into coordinate system where ray is [0, 0]
				const long long ak = a0k - rayK;
				const long long al = a0l - rayL;
				const long long bk = b0k - rayK;
				const long long bl = b0l - rayL;
				const long long ck = c0k - rayK;
				const long long cl = c0l - rayL;

				bool rayHit = true;
				
				rayHitDetector(ak, al, bk, bl, ck, cl, rayHit);

				if (rayHit) 
				{
					TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(intersectionCounterArrayView(i, j), 1);
				}
			}
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STLArbeiter.triangleCount, counterLambda );	
	
	auto fetch = [ = ] __cuda_callable__( const size_t singleIndex )
	{
		const size_t i = singleIndex % cellCount.nx;
		const size_t j = singleIndex / cellCount.ny;
		return intersectionCounterArrayView(i, j);
	};
	auto reduction = [] __cuda_callable__( const int& a, const int& b )
	{
		return TNL::max( a, b );
	};
	size_t start = 0;
	size_t end = cellCount.nx * cellCount.ny;
	int intersectionCountMax = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.0 );
	std::cout << "intersectionCountMax: " << intersectionCountMax << std::endl; 
	
	//test
	for (size_t j = 0; j < cellCount.ny; j++) 
	{
		for (size_t i = 0; i < cellCount.nx; i++) 
		{
			int number = intersectionCounterArray.getElement(i, j);
			if (number % 2 != 0) 
			{
				std::cout << "Intersection count for i: " << i << ", j: " << j << ", : " << number << std::endl;
			}
		}
	}	
	
	// another test
	for (size_t triangleIndex = 0; triangleIndex < STLArbeiter.triangleCount; triangleIndex++ )
    {
		const float a0x = STLArbeiter.axArray.getElement(triangleIndex);
		const float a0y = STLArbeiter.ayArray.getElement(triangleIndex);
		const float b0x = STLArbeiter.bxArray.getElement(triangleIndex);
		const float b0y = STLArbeiter.byArray.getElement(triangleIndex);
		const float c0x = STLArbeiter.cxArray.getElement(triangleIndex);
		const float c0y = STLArbeiter.cyArray.getElement(triangleIndex);
		// transform STL floats to integer grid that is 100x finer than the LBM grid to prevent float errors
		// transform into coordinate system of the LBM grid
		// make the STL coords odd, rays will be even, this prevents hitting a vortex
		const float scale = 50.0f / cellCount.res;
		const long long ok = round( cellCount.ox * scale );
		const long long ol = round( cellCount.oy * scale );
		const long long a0k = (long long)(round( a0x * scale ) - ok) * 2 + 1;
		const long long a0l = (long long)(round( a0y * scale ) - ol) * 2 + 1;
		const long long b0k = (long long)(round( b0x * scale ) - ok) * 2 + 1;
		const long long b0l = (long long)(round( b0y * scale ) - ol) * 2 + 1;
		const long long c0k = (long long)(round( c0x * scale ) - ok) * 2 + 1;
		const long long c0l = (long long)(round( c0y * scale ) - ol) * 2 + 1;
		
		const long long kmin = std::min({ a0k, b0k, c0k });
		const long long kmax = std::max({ a0k, b0k, c0k });
		const long long lmin = std::min({ a0l, b0l, c0l });
		const long long lmax = std::max({ a0l, b0l, c0l });
		
		const size_t imin = min( (size_t)max(0LL, (kmin / 100) ), cellCount.nx-1 );
		const size_t imax = min( (size_t)max(0LL, (kmax / 100) + 1), cellCount.nx-1 );
		const size_t jmin = min( (size_t)max(0LL, (lmin / 100) ), cellCount.ny-1 );
		const size_t jmax = min( (size_t)max(0LL, (lmax / 100) + 1), cellCount.ny-1 );
		
		for (size_t j = jmin; j <= jmax; j++) 
		{
			for (size_t i = imin; i <= imax; i++) 
			{
				const long long rayK = i * 100;
				const long long rayL = j * 100;
				// finally, transform the triangle into coordinate system where ray is [0, 0]
				const long long ak = a0k - rayK;
				const long long al = a0l - rayL;
				const long long bk = b0k - rayK;
				const long long bl = b0l - rayL;
				const long long ck = c0k - rayK;
				const long long cl = c0l - rayL;

				bool rayHit = true;
				
				rayHitDetector(ak, al, bk, bl, ck, cl, rayHit);

				if (rayHit) 
				{
					if (i==57 && j==121)
					{
						std::cout << "i: " << i << ", j: " << j << std::endl;
						std::cout << "triangleIndex: " << triangleIndex << std::endl;
						std::cout << "ak: " << ak << ", al: " << al << std::endl;
						std::cout << "bk: " << bk << ", bl: " << bl << std::endl;
						std::cout << "ck: " << ck << ", cl: " << cl << std::endl;
					}
				}
			}
		}
	}		     
}

void checkSTLEdges( STLArbeiterStruct &STLArbeiter )
{
	auto axArrayView = STLArbeiter.axArray.getConstView();
	auto ayArrayView = STLArbeiter.ayArray.getConstView();
	auto azArrayView = STLArbeiter.azArray.getConstView();
	auto bxArrayView = STLArbeiter.bxArray.getConstView();
	auto byArrayView = STLArbeiter.byArray.getConstView();
	auto bzArrayView = STLArbeiter.bzArray.getConstView();
	auto cxArrayView = STLArbeiter.cxArray.getConstView();
	auto cyArrayView = STLArbeiter.cyArray.getConstView();
	auto czArrayView = STLArbeiter.czArray.getConstView();
	
	IntArrayType ABcounterArray = IntArrayType(STLArbeiter.triangleCount, 0);
	IntArrayType BCcounterArray = IntArrayType(STLArbeiter.triangleCount, 0);
	IntArrayType CAcounterArray = IntArrayType(STLArbeiter.triangleCount, 0);
	auto ABcounterArrayView = ABcounterArray.getView();
	auto BCcounterArrayView = BCcounterArray.getView();
	auto CAcounterArrayView = CAcounterArray.getView();

    auto counterLambda = [ = ] __cuda_callable__( size_t triangleIndex ) mutable
    {		
		const float ax1 = axArrayView[ triangleIndex ];
		const float ay1 = ayArrayView[ triangleIndex ];
		const float bx1 = bxArrayView[ triangleIndex ];
		const float by1 = byArrayView[ triangleIndex ];
		const float cx1 = cxArrayView[ triangleIndex ];
		const float cy1 = cyArrayView[ triangleIndex ];
		
		for (size_t comparedTriangleIndex = 0; comparedTriangleIndex < STLArbeiter.triangleCount; comparedTriangleIndex++) 
		{
			const float ax2 = axArrayView[ comparedTriangleIndex ];
			const float ay2 = ayArrayView[ comparedTriangleIndex ];
			const float bx2 = bxArrayView[ comparedTriangleIndex ];
			const float by2 = byArrayView[ comparedTriangleIndex ];
			const float cx2 = cxArrayView[ comparedTriangleIndex ];
			const float cy2 = cyArrayView[ comparedTriangleIndex ];
			
			// ---------------------------------------------------------
			// CHECK EDGE AB (Triangle 1)
			// ---------------------------------------------------------
			// vs AB2
			if (ax1 == ax2 && ay1 == ay2 && bx1 == bx2 && by1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// vs BA2
			if (ax1 == bx2 && ay1 == by2 && bx1 == ax2 && by1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// vs BC2
			if (ax1 == bx2 && ay1 == by2 && bx1 == cx2 && by1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// vs CB2
			if (ax1 == cx2 && ay1 == cy2 && bx1 == bx2 && by1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// vs CA2
			if (ax1 == cx2 && ay1 == cy2 && bx1 == ax2 && by1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// vs AC2
			if (ax1 == ax2 && ay1 == ay2 && bx1 == cx2 && by1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(ABcounterArrayView(triangleIndex), 1);
			// ---------------------------------------------------------
			// CHECK EDGE BC (Triangle 1)
			// ---------------------------------------------------------
			// vs AB2
			if (bx1 == ax2 && by1 == ay2 && cx1 == bx2 && cy1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// vs BA2
			if (bx1 == bx2 && by1 == by2 && cx1 == ax2 && cy1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// vs BC2
			if (bx1 == bx2 && by1 == by2 && cx1 == cx2 && cy1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// vs CB2
			if (bx1 == cx2 && by1 == cy2 && cx1 == bx2 && cy1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// vs CA2
			if (bx1 == cx2 && by1 == cy2 && cx1 == ax2 && cy1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// vs AC2
			if (bx1 == ax2 && by1 == ay2 && cx1 == cx2 && cy1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(BCcounterArrayView(triangleIndex), 1);
			// ---------------------------------------------------------
			// CHECK EDGE CA (Triangle 1)
			// ---------------------------------------------------------
			// vs AB2
			if (cx1 == ax2 && cy1 == ay2 && ax1 == bx2 && ay1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
			// vs BA2
			if (cx1 == bx2 && cy1 == by2 && ax1 == ax2 && ay1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
			// vs BC2
			if (cx1 == bx2 && cy1 == by2 && ax1 == cx2 && ay1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
			// vs CB2
			if (cx1 == cx2 && cy1 == cy2 && ax1 == bx2 && ay1 == by2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
			// vs CA2
			if (cx1 == cx2 && cy1 == cy2 && ax1 == ax2 && ay1 == ay2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
			// vs AC2
			if (cx1 == ax2 && cy1 == ay2 && ax1 == cx2 && ay1 == cy2) TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(CAcounterArrayView(triangleIndex), 1);
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STLArbeiter.triangleCount, counterLambda );	

	for (size_t triangleIndex = 0; triangleIndex < STLArbeiter.triangleCount; triangleIndex++ )
    {
		int ABcount = ABcounterArray.getElement(triangleIndex);
		int BCcount = BCcounterArray.getElement(triangleIndex);
		int CAcount = CAcounterArray.getElement(triangleIndex);
		if (ABcount != 2 || BCcount != 2 || CAcount != 2)
		{
			std::cout << "	Edge problem on triangle " << triangleIndex << ", ABcount: " << ABcount << ", BCcount: " << BCcount << ", CAcount: " << CAcount << std::endl;
		}
	}
		     
}
