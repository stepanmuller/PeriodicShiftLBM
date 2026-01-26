void readSTL( STLArbeiterStructCPU& STLArbeiterCPU )
{
	std::cout << "Reading STL" << std::endl;
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

void checkSTLEdges( STLArbeiterStruct &STLArbeiter )
{
	std::cout << "Starting STL check for shared edges" << std::endl;
	auto axArrayView = STLArbeiter.axArray.getConstView();
	auto ayArrayView = STLArbeiter.ayArray.getConstView();
	auto azArrayView = STLArbeiter.azArray.getConstView();
	auto bxArrayView = STLArbeiter.bxArray.getConstView();
	auto byArrayView = STLArbeiter.byArray.getConstView();
	auto bzArrayView = STLArbeiter.bzArray.getConstView();
	auto cxArrayView = STLArbeiter.cxArray.getConstView();
	auto cyArrayView = STLArbeiter.cyArray.getConstView();
	auto czArrayView = STLArbeiter.czArray.getConstView();
	
	CounterArray2DType edgeCounterArray;
	edgeCounterArray.setSizes(STLArbeiter.triangleCount, 3);
	edgeCounterArray.setValue(0);
	auto edgeCounterArrayView = edgeCounterArray.getView();

    auto counterLambda = [ = ] __cuda_callable__( size_t triangle1Index ) mutable
    {		
		float triangle1[9];
		triangle1[0] = axArrayView[ triangle1Index ];
		triangle1[1] = ayArrayView[ triangle1Index ];
		triangle1[2] = azArrayView[ triangle1Index ];
		triangle1[3] = bxArrayView[ triangle1Index ];
		triangle1[4] = byArrayView[ triangle1Index ];
		triangle1[5] = bzArrayView[ triangle1Index ];
		triangle1[6] = cxArrayView[ triangle1Index ];
		triangle1[7] = cyArrayView[ triangle1Index ];
		triangle1[8] = czArrayView[ triangle1Index ];
		
		for (size_t triangle2Index = 0; triangle2Index < STLArbeiter.triangleCount; triangle2Index++) 
		{
			float triangle2[9];
			triangle2[0] = axArrayView[ triangle2Index ];
			triangle2[1] = ayArrayView[ triangle2Index ];
			triangle2[2] = azArrayView[ triangle2Index ];
			triangle2[3] = bxArrayView[ triangle2Index ];
			triangle2[4] = byArrayView[ triangle2Index ];
			triangle2[5] = bzArrayView[ triangle2Index ];
			triangle2[6] = cxArrayView[ triangle2Index ];
			triangle2[7] = cyArrayView[ triangle2Index ];
			triangle2[8] = czArrayView[ triangle2Index ];
			
			for (int edge1 = 0; edge1 < 3; edge1++)
			{
				for (int edge2 = 0; edge2 < 3; edge2++)
				{
					float ax1 = triangle1[3*edge1];
					float ay1 = triangle1[3*edge1+1];
					float az1 = triangle1[3*edge1+2];
					float bx1 = triangle1[3*((edge1+1)%3)];
					float by1 = triangle1[3*((edge1+1)%3)+1];
					float bz1 = triangle1[3*((edge1+1)%3)+2];
					
					float ax2 = triangle2[3*edge2];
					float ay2 = triangle2[3*edge2+1];
					float az2 = triangle2[3*edge2+2];
					float bx2 = triangle2[3*((edge2+1)%3)];
					float by2 = triangle2[3*((edge2+1)%3)+1];
					float bz2 = triangle2[3*((edge2+1)%3)+2];
					// testing same orientation of the edge
					if (ax1 == ax2 && ay1 == ay2 && az1 == az2 && bx1 == bx2 && by1 == by2 && bz1 == bz2) 
					{
						TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(edgeCounterArrayView(triangle1Index, edge1), 1);
					}
					// testing reverse orientation of the edge
					if (ax1 == bx2 && ay1 == by2 && az1 == bz2 && bx1 == ax2 && by1 == ay2 && bz1 == az2) 
					{
						TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(edgeCounterArrayView(triangle1Index, edge1), 1);
					}
				}
			}
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STLArbeiter.triangleCount, counterLambda );
	int errorCounter = 0;
	for (size_t triangleIndex = 0; triangleIndex < STLArbeiter.triangleCount; triangleIndex++ )
    {
		int ABcount = edgeCounterArray.getElement(triangleIndex, 0);
		int BCcount = edgeCounterArray.getElement(triangleIndex, 1);
		int CAcount = edgeCounterArray.getElement(triangleIndex, 2);
		if (ABcount != 2 || BCcount != 2 || CAcount != 2)
		{
			errorCounter++;
			std::cout << "	Edge problem on triangle " << triangleIndex << ", ABcount: " << ABcount << ", BCcount: " << BCcount << ", CAcount: " << CAcount << std::endl;
		}
	}    
	std::cout<< "	Total shared edge problems: " << errorCounter << std::endl; 
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

__host__ __device__ void rayHitCoordinate(	const float &ax, const float &ay, const float &az,
											const float &bx, const float &by, const float &bz,
											const float &cx, const float &cy, const float &cz,
											const float &rayX, const float &rayY, float &rayZ)
{
    float v1x = bx - ax;
    float v1y = by - ay;
    float v1z = bz - az;

    float v2x = cx - ax;
    float v2y = cy - ay;
    float v2z = cz - az;

    float nx = v1y * v2z - v1z * v2y;
    float ny = v1z * v2x - v1x * v2z;
    float nz = v1x * v2y - v1y * v2x;

    if (nz != 0.0f) rayZ = az - (nx * (rayX - ax) + ny * (rayY - ay)) / nz;
    else rayZ = std::numeric_limits<float>::max();
}

void applyMarkersFromSTL( MarkerStruct &Marker, STLArbeiterStruct &STLArbeiter, CellCountStruct &cellCount )
{
	std::cout << "Applying markers from STL" << std::endl;
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
	std::cout << "	intersectionCountMax: " << intersectionCountMax << std::endl; 
	
	// check for odd intersections
	int oddIntersectionCounter = 0;
	for (size_t j = 0; j < cellCount.ny; j++) 
	{
		for (size_t i = 0; i < cellCount.nx; i++) 
		{
			int number = intersectionCounterArray.getElement(i, j);
			if (number % 2 != 0) 
			{
				std::cout << "	Odd intersection count for i: " << i << ", j: " << j << ", : " << number << std::endl;
				oddIntersectionCounter++;
			}
		}
	}		
	std::cout<< "	Total rays with odd intersection count: " << oddIntersectionCounter << std::endl;      
	
	// extracting the intersection indexes
	IndexArray3DType intersectionIndexArray;
	intersectionIndexArray.setSizes(intersectionCountMax, cellCount.nx, cellCount.ny);
	intersectionIndexArray.setValue(cellCount.nz+1);
	auto intersectionIndexArrayView = intersectionIndexArray.getView();
	intersectionCounterArray.setValue(0);
	
	auto rayHitIndexLambda = [ = ] __cuda_callable__( size_t triangleIndex ) mutable
    {
		const float a0x = axArrayView[ triangleIndex ];
		const float a0y = ayArrayView[ triangleIndex ];
		const float a0z = azArrayView[ triangleIndex ];
		const float b0x = bxArrayView[ triangleIndex ];
		const float b0y = byArrayView[ triangleIndex ];
		const float b0z = bzArrayView[ triangleIndex ];
		const float c0x = cxArrayView[ triangleIndex ];
		const float c0y = cyArrayView[ triangleIndex ];
		const float c0z = czArrayView[ triangleIndex ];
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
					const int writePosition = TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(intersectionCounterArrayView(i, j), 1);
					// here I need to find the Z coordinate of the hit and the nearest larger k index in Z
					const float rayX = i * cellCount.res + cellCount.ox;
					const float rayY = j * cellCount.res + cellCount.oy;
					float rayZ = std::numeric_limits<float>::max();
					rayHitCoordinate(a0x, a0y, a0z, b0x, b0y, b0z, c0x, c0y, c0z, rayX, rayY, rayZ);
					size_t k = (size_t)max(0.0f, ceilf((rayZ - cellCount.oz) / cellCount.res));
					k = min(k, cellCount.nz+1);
					intersectionIndexArrayView(writePosition, i, j) = k;
				}
			}
		}
	};
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STLArbeiter.triangleCount, rayHitIndexLambda );
	
	auto rayLambda = [=] __cuda_callable__ (const TNL::Containers::StaticArray< 2, int >& doubleIndex) mutable
	{
		const size_t i = doubleIndex.x();
		const size_t j = doubleIndex.y();
		
		for (int layer = 1; layer < intersectionCountMax; layer++) 
		{
			size_t key = intersectionIndexArrayView(layer, i, j);
			int slider = layer - 1;
			while (slider >= 0 && intersectionIndexArrayView(slider, i, j) > key) 
			{
				intersectionIndexArrayView(slider + 1, i, j) = intersectionIndexArrayView(slider, i, j);
				slider = slider - 1;
			}
			intersectionIndexArrayView(slider + 1, i, j) = key;
		}
		bool markerValue = 1;
		for (int interval = 0; interval <= intersectionCountMax; interval++)
		{
			size_t start = 0;
			size_t end = 0;
			if (interval == 0) end = intersectionIndexArrayView(0, i, j);
			else if (interval == intersectionCountMax) 
			{
				start = intersectionIndexArrayView(intersectionCountMax-1, i, j);
				end = cellCount.nz;
			}
			else
			{
				start = intersectionIndexArrayView(interval-1, i, j);
				end = intersectionIndexArrayView(interval, i, j);
			}
			size_t runner = start;
			while (runner < end)
			{
				size_t cell = convertIndex(i, j, runner, cellCount);
				bouncebackMarkerArrayView[cell] = markerValue;
				runner++;
			}
			markerValue = !markerValue;
		}		
	};
	TNL::Containers::StaticArray< 2, size_t > startList{ 0, 0 };
	TNL::Containers::StaticArray< 2, size_t > endList{ cellCount.nx, cellCount.ny };
	TNL::Algorithms::parallelFor<TNL::Devices::Cuda>(startList, endList, rayLambda );	
}
