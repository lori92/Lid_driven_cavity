#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include <sstream>
#include <iomanip>
 
#include <chrono>

double* solve(double* A, double* B, double* C, double* D, int n) {

    double* upperDiag = C;
    double*  Diag      = B;
    double*  lowerDiag = A;

    double*  b         = D;

    double* x = new double[n];

    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    */
    n--; // since we start from x0 (not x1)
    upperDiag[0] /= Diag[0];
    b[0] /= Diag[0];

    for (int i = 1; i < n; i++) {
        upperDiag[i] /= Diag[i] - lowerDiag[i] * upperDiag[i-1];
        b[i] = (b[i] - lowerDiag[i]*b[i-1]) / (Diag[i] - lowerDiag[i]*upperDiag[i-1]);
    }

    b[n] = (b[n] - lowerDiag[n]*b[n-1]) / (Diag[n] - lowerDiag[n]*upperDiag[n-1]);

    for (int i = n; i-- > 0;) {
        b[i] -= upperDiag[i]*b[i+1];
    }
    return b;
}

int main()
{
   const int n = 3;
   double A[n][n]={0};
   double b[n]={1};
   double* x = new double[n];

  
   double upperD[n];
   double mainD[n] ;
   double lowD[n] ;

    for (int k=0; k<n; k++)
   {
    upperD[k] = 1.;
    mainD[k] = 1;
    lowD[k] = 1;
   }


   for (int k=0; k<n; k++)
   {
    std::cout << lowD[k]<<"\n";

   }
      x = solve( lowD, mainD, upperD, b, n);

      
   for (int k=0; k<n; k++)
   {
    std::cout << lowD[k]<<"\n";
    std::cout << x[k]<<"\n";

   }



    return 0;
}