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

// cell **Uy   = buildMesh_v(Nx, Ny, Lx, Ly);
// cell **p    = buildMesh(Nx, Ny, Lx, Ly);

 mesh[1][Nx+1].print_cellCentre();

 initMeshVal(Ux, Nx, Ny, 12.5);
 Ux  [1][Nx+1].print_cellVal();

 free(mesh);
 free(Ux);
 //free(Uy);
 //free(p);
 return 0;	
}
