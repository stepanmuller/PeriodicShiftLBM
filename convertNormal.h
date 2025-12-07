#pragma once

__host__ __device__ void encodeNormal(  uint8_t& normalCode, const int nx, const int ny, const int nz )
{
    uint8_t ex, ey, ez;

    // encode nx
    ex = 0u;
    if (nx == -1) ex = 1u;
    else if (nx == 1) ex = 2u;

    // encode ny
    ey = 0u;
    if (ny == -1) ey = 1u;
    else if (ny == 1) ey = 2u;

    // encode nz
    ez = 0u;
    if (nz == -1) ez = 1u;
    else if (nz == 1) ez = 2u;

    normalCode = ex | (ey << 2) | (ez << 4);
}

__host__ __device__ void decodeNormal(const uint8_t normalCode, int& nx, int& ny, int& nz)
{
    uint8_t ex, ey, ez;

    ex = (normalCode >> 0) & 3u;
    ey = (normalCode >> 2) & 3u;
    ez = (normalCode >> 4) & 3u;

    // decode nx
    nx = 0;
    if (ex == 1u) nx = -1;
    else if (ex == 2u) nx = 1;

    // decode ny
    ny = 0;
    if (ey == 1u) ny = -1;
    else if (ey == 2u) ny = 1;

    // decode nz
    nz = 0;
    if (ez == 1u) nz = -1;
    else if (ez == 2u) nz = 1;
}
