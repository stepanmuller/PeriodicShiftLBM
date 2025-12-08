#pragma once

#include "config.h"

__host__ __device__ void getOuterNormal(const size_t& cell, short& outerNormalX, short& outerNormalY, short& outerNormalZ)
{
    const size_t xy = cellCountX * cellCountY;
    const size_t k = cell / xy;
    size_t remainder = cell % xy;
    const size_t j = remainder / cellCountX;
    const size_t i = remainder % cellCountX;
    outerNormalX = 0;
    outerNormalY = 0;
    outerNormalZ = 0;
    if 			( i == 0 ) 				outerNormalX = -1;
    else if 	( i == cellCountX - 1 ) outerNormalX = 1;
    if 			( j == 0 ) 				outerNormalY = -1;
    else if 	( j == cellCountY - 1) 	outerNormalY = 1;
    if 			( k == 0 ) 				outerNormalZ = -1;
    else if 	( k == cellCountZ-1 ) 	outerNormalZ = 1;
}
