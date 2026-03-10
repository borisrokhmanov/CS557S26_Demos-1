#include "Part2KernelGroup.h"

void GroupLine2(const float (&a)[XDIM][YDIM][ZDIM], const float (&b)[XDIM][YDIM][ZDIM],
                float (&c)[XDIM][YDIM][ZDIM], float (&d)[XDIM][YDIM][ZDIM])
{
    #pragma omp parallel for collapse(3)
        for (int i = 1; i < XDIM-1; i++){
            for (int j = 1; j < YDIM-1; j++){
                for (int k = 1; k < ZDIM-1; k++){
                    c[i][j][k] =
                    (-6 * a[i][j][k]
                    + a[i+1][j][k]
                    + a[i-1][j][k]
                    + a[i][j+1][k]
                    + a[i][j-1][k]
                    + a[i][j][k+1]
                    + a[i][j][k-1]);
                    d[i][j][k] = -1 * c[i][j][k] + b[i][j][k];
                }
            }
        }
}

float GroupLine6(const float (&p)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM])
{
    double result = 0.;

#pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM-1; i++){
        for (int j = 1; j < YDIM-1; j++){
            for (int k = 1; k < ZDIM-1; k++){
                z[i][j][k] = -6 * p[i][j][k]
                    + p[i+1][j][k]
                    + p[i-1][j][k]
                    + p[i][j+1][k]
                    + p[i][j-1][k]
                    + p[i][j][k+1]
                    + p[i][j][k-1];
                result += (double) p[i][j][k] * (double) z[i][j][k];
            }
        }
    }
    
    return (float) result;
}