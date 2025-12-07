#pragma once

using ArrayType = TNL::Containers::Array< float, TNL::Devices::Cuda, size_t >;
using IndexArrayType = TNL::Containers::Array< size_t, TNL::Devices::Cuda, size_t >;
using normalArrayType = TNL::Containers::Array< uint8_t, TNL::Devices::Cuda, size_t >;

struct DistributionFunctionStruct { IndexArrayType shifter; ArrayType fArray[27]; };

struct RhoUGStruct { 	ArrayType rhoArray; 
						ArrayType uxArray; ArrayType uyArray; ArrayType uzArray; 
						ArrayType gxArray; ArrayType gyArray; ArrayType gzArray; };

struct CellGroupStruct { size_t groupSize = 0; size_t temp = 0; IndexArrayType indexArray; IndexArrayType sourceIndexArray; normalArrayType normalArray; };
