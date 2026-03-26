#include "Reductions.h"

#include <algorithm>

//#define DO_NOT_USE_MKL
#ifndef DO_NOT_USE_MKL
#include <mkl.h>
#endif

float Norm(const float (&x)[XDIM][YDIM][ZDIM])
{
#ifdef DO_NOT_USE_MKL
    float result = 0.;

#pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        result = std::max(result, std::abs(x[i][j][k]));

    return result;
#else
    int result = cblas_isamax(
        XDIM * YDIM * ZDIM, //# of elements in x
        &x[0][0][0],        //input vector x
        1                   //step size of 1
    );
    return std::abs((&x[0][0][0])[result]); //isamax returns index of max value
#endif
}

float InnerProduct(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM])
{
#ifdef DO_NOT_USE_MKL
    double result = 0.;

#pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        result += (double) x[i][j][k] * (double) y[i][j][k];

    return (float) result;
#else
    float result = cblas_sdot(
        XDIM * YDIM * ZDIM, //# of elements in x and y
        &x[0][0][0],        //input vector x
        1,                  //step size 1 for x
        &y[0][0][0],        //input vector y
        1                   //step size 1 for y
    );
    return result;
#endif
}
