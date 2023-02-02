#include <iostream>
#include <cmath>
#include "mesh.cpp"
using namespace std;

/// number of cells (without ghost nodes) ///
#define Nx 5
#define Ny 5

int main()
{

 const double Lx = 1;
 const double Ly = 1;

 cell **mesh = buildMesh(Nx, Ny, Lx, Ly);
 cell **Ux   = buildMesh(Nx, Ny, Lx, Ly);

 cell **Uy   = buildMesh(Nx, Ny, Lx, Ly);
// cell **p    = buildMesh(Nx, Ny, Lx, Ly);

 //mesh[1][Nx+1].print_cellCentre();

 initFieldVal(Ux, Nx, Ny, 12.5);
 initFieldVal(Uy, Nx, Ny,  2.5);

 setBCu(Ux, Nx, Ny);
 setBCv(Uy, Nx, Ny);

/*
 for (int i = 0; i <= Nx+1; i++)
 {
   for(int j = 0; j<= Ny+1; j++)
   {
    Ux[i][j].print_cellVal();
   }
   cout<<"\n";
 } 
*/
  for (int j = 1; j <= Ny; j++)
 {
   for(int i = 0; i<= Nx; i++)
   {
    Uy[i][j].print_cellVal();
   }
   cout<<"\n";
 } 



 //free(Uy);
 //free(p);
double fUU_f, fUU_b;
double fUV_f, fUV_b;

double conv_u;

 for (int i = 1; i <= Nx; i++)
 {
   for(int j = 1; j<= Ny; j++)
   {
     fUU_f = 0.25 * (Ux[i][j].cellVal() + Ux[i][j].cellVal()) *(Ux[i+1][j].cellVal() + Ux[i+1][j].cellVal());
     fUU_b = 0.25 * (Ux[i][j].cellVal() + Ux[i][j].cellVal()) *(Ux[i-1][j].cellVal() + Ux[i-1][j].cellVal());

     fUV_f = 0.5*(Ux[i ][j+1].cellVal() + Ux[i][j].cellVal() ) * \
                                  0.5* (  Uy[i][ j ].cellVal() + Uy[i+1][ j ].cellVal() );
     fUV_b = 0.5*(Ux[i ][j ].cellVal() + Ux[i-1][j-1].cellVal() ) * \
                                  0.5* (  Uy[i][ j ].cellVal() + Uy[i+1][j-1].cellVal() );

     conv_u = (fUU_f - fUU_b) / (mesh[i][j].dx())  \
             +(fUV_f - fUV_b) / (mesh[i][j].dy()); 

   }
   cout<<"\n";
 } 




 free(mesh);
 free(Ux);
 free(Uy);
 return 0;	
}
