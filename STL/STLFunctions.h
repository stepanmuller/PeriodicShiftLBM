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
		float ax = axArrayView[ triangleIndex ];
		float ay = ayArrayView[ triangleIndex ];
		float az = azArrayView[ triangleIndex ];
		float bx = bxArrayView[ triangleIndex ];
		float by = byArrayView[ triangleIndex ];
		float bz = bzArrayView[ triangleIndex ];
		float cx = cxArrayView[ triangleIndex ];
		float cy = cyArrayView[ triangleIndex ];
		float cz = czArrayView[ triangleIndex ];
		
		float xmin = fminf(ax, fminf(bx, cx));
		float ymin = fminf(ay, fminf(by, cy));
		float xmax = fmaxf(ax, fmaxf(bx, cx)); 
		float ymax = fmaxf(ay, fmaxf(by, cy));
		size_t imin = (size_t)max(0LL, (long long)((xmin - cellCount.ox) / cellCount.res));
		size_t jmin = (size_t)max(0LL, (long long)((ymin - cellCount.oy) / cellCount.res));
		size_t imax = (size_t)min((long long)cellCount.nx - 1, (long long)((xmax - cellCount.ox) / cellCount.res) + 1);
		size_t jmax = (size_t)min((long long)cellCount.ny - 1, (long long)((ymax - cellCount.oy) / cellCount.res) + 1);
		
		for (size_t j = jmin; j <= jmax; j++) 
		{
			for (size_t i = imin; i <= imax; i++) 
			{
				float rayX = i * cellCount.res + cellCount.ox;
				float rayY = j * cellCount.res + cellCount.oy;

				// 2. 2D Edge Functions (Point-in-Triangle test in the XY plane)
				auto edgeFunction = [](float px, float py, float x0, float y0, float x1, float y1) {
					return (px - x0) * (y1 - y0) - (py - y0) * (x1 - x0);
				};

				float w0 = edgeFunction(rayX, rayY, bx, by, cx, cy);
				float w1 = edgeFunction(rayX, rayY, cx, cy, ax, ay);
				float w2 = edgeFunction(rayX, rayY, ax, ay, bx, by);

				// 3. Determine if the ray center is inside the triangle projection
				bool intersectionHappens = (w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0);

				if ( intersectionHappens ) 
				{
					TNL::Algorithms::AtomicOperations<TNL::Devices::Cuda>::add(intersectionCounterArrayView(i, j), 1);
				}
			}
		}
    };
    TNL::Algorithms::parallelFor<TNL::Devices::Cuda>( 0, STLArbeiter.triangleCount, counterLambda );	    
}
