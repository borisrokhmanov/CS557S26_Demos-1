#pragma once

#include "Parameters.h"

//compute laplacian of u, scale by -1, add y, write result into z
void GroupLine2(const float (&u)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
                float (&z)[XDIM][YDIM][ZDIM], float (&x)[XDIM][YDIM][ZDIM]);

//compute laplacian of p, return inner product of p and z
float GroupLine6(const float (&p)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM]);