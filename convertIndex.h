#pragma once

#include "config.h"

__host__ __device__ size_t convertIndex(size_t i, size_t j, size_t k)
{
	size_t cell = k * (cellCountX * cellCountY) + j * cellCountX + i;
	return cell;
}
