#include "config.h"

size_t convertIndex(size_t i, size_t j, size_t k)
{
	size_t cell = k * (cellCountX * cellCountY) + j * cellCountX + i;
	return cell;
}
