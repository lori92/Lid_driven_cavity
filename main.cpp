#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

void set_bc_u(double** u, int dim_x, int dim_y )
{
  for (int i=0; i<dim_x; i++) 
  { u[i][0] = 0; }
};


int main()
{

const double Lx = 1;
const double Ly = 1;

const int  Nx = 10;
const int  Ny = 10;

double *x = new double[Nx+2];
double *y = new double[Ny+2];

double** u = new double*[Nx+2];
for (int i = 0; i <= Nx+1; i++)
{
 u[i] = new double[Ny+2];
 for (int j = 0; j <= Ny+2; j++) 
 { u[i][j] = 1.0; }
}

set_bc_u(u,Nx+2,Ny+2);

/////// 
for (int i = 0; i<=Nx+1; i++) 
{ x[i] = 0 + (Lx/Nx)*i; }

cout << "x - position:" <<x[9]<<"\n";

for (int i = 0; i<Nx+2; i++) {
 for (int j = 0; j<Ny+2; j++)  {
  cout <<u[i][j]<<" ";
 }
 cout << "\n";
}

free(x);
free(y);
return 0;
}