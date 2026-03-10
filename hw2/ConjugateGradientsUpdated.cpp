#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"

#include "GroupLine2.h"

#include <iostream>

Timer ComputeLaplacianTimer;
Timer InnerProductTimer;
Timer NormTimer;
Timer SaxpyTimer;
Timer CopyTimer;

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
    // Algorithm : Line 2
    ComputeLaplacianTimer.Start();
    ComputeLaplacian(x, z);
    ComputeLaplacianTimer.Pause();
    SaxpyTimer.Start();
    Saxpy(z, f, r, -1);
    SaxpyTimer.Pause();
    NormTimer.Start();
    float nu = Norm(r);
    NormTimer.Pause();

    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    CopyTimer.Start();
    Copy(r, p);
    CopyTimer.Pause();
    InnerProductTimer.Start();
    float rho=InnerProduct(p, r);
    InnerProductTimer.Pause();
        
    // Beginning of loop from Line 5
    for(int k=0;;k++)
    {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        ComputeLaplacianTimer.Restart();
        ComputeLaplacian(p, z);
        ComputeLaplacianTimer.Pause();
        InnerProductTimer.Restart();
        float sigma=InnerProduct(p, z);
        InnerProductTimer.Pause();

        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        SaxpyTimer.Restart();
        Saxpy(z, r, r, -alpha);
        SaxpyTimer.Pause();
        NormTimer.Restart();
        nu=Norm(r);
        NormTimer.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            SaxpyTimer.Restart();
            Saxpy(p, x, x, alpha);
            SaxpyTimer.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            //disabled for now
            //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            ComputeLaplacianTimer.Print("ComputeLaplacian kernel: ");
            InnerProductTimer.Print("Inner Product kernel: ");
            NormTimer.Print("Norm kernel: ");
            SaxpyTimer.Print("Saxpy kernel: ");
            CopyTimer.Print("Copy kernel: ");
            return;
        }
            
        // Algorithm : Line 13
        CopyTimer.Restart();
        Copy(r, z);
        CopyTimer.Pause();
        InnerProductTimer.Restart();
        float rho_new = InnerProduct(z, r);
        InnerProductTimer.Pause();

        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        SaxpyTimer.Restart();
        Saxpy(p, x, x, alpha);
        Saxpy(p, r, p, beta);
        SaxpyTimer.Pause();

        //disabled for now
        //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
}
