#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "global_variables.h"

#include "mesh.cpp"
#include "poisson_solver.cpp"
#include "finiteVolume.cpp"
#include "predictor_step.cpp"

#include "write_output.C"

using namespace std;

/// number of internal cells (without ghost nodes) ///
#define Nx 200
#define Ny 200

int main()
{

 const double Lx = 1;
 const double Ly = 1;

 FVcell **mesh      = buildMesh(Nx, Ny, Lx, Ly);

 volumeField p      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField pmg      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField p_sor    = buildFVfield(Nx, Ny, Lx, Ly);

 volumeField Ux     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField div    = buildFVfield(Nx, Ny, Lx, Ly);

 initFieldVal(p , 0);
 initFieldVal(pmg, 0);
 initFieldVal(p_sor , 0);
 initFieldVal(div , 0);
 initFieldVal(Ux, 0);
 initFieldVal(Uy, 0);

// write first initialized field 
 std::stringstream ss;
 ss <<"test_poisson";
 std::string namefile = ss.str();

 cout<<"\n";

 //////// solve Poisson Equation for Pressure /////////
 //if (n_iter_ssor >100 &&  n_iter>4000)
 //{  // multigrid(&p, div, 1e-14, rho, dt, 1.3, Lx, Ly, 2);
//   multigrid(&p, div, 1e-14, 1, 1, 1.3, Lx, Ly, 1e-3,1e-3,2);
  int vcycle = 0;
  // multigrid_rec(&pmg, div, 1e-14,  1.3, Lx, Ly, 1e-3,1e-3, vcycle);

  //}
  //else {
  poisson(&p_sor, div, 1, 1, 1e-14, 1.3, n_iter_ssor);
   //cout << "convergence achieved on pressure in "<< n_iter_ssor <<" iterations\n";//}

 

     // interpolate on the coarser grid //
 volumeField  Ux_coarser = buildFVfield(Nx/2, Ny/2, Lx, Ly);
 volumeField  Uy_coarser = buildFVfield(Nx/2, Ny/2, Lx, Ly);
 volumeField  p_coarser  = buildFVfield(Nx/2, Ny/2, Lx, Ly);

 interpolateFieldVal(Ux , &Ux_coarser); 
 interpolateFieldVal(Uy , &Uy_coarser); 
 interpolateFieldVal(p  , &p_coarser); 
 write_output(Ux_coarser, Uy_coarser, p_coarser, "0_coarse_grid");


      // prolongate on the coarser grid //
 volumeField  Ux_prl = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  Uy_prl = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  p_prl  = buildFVfield(Nx, Ny, Lx, Ly);

 prolongateFieldVal(Ux_coarser , &Ux_prl); 
 prolongateFieldVal(Uy_coarser , &Uy_prl); 
 prolongateFieldVal(p_coarser , &p_prl); 
 write_output(Ux_prl, Uy_prl, p_prl, "0_fine_grid");


 write_output(Ux, Uy, p, "0_1cycle");
 write_output(Ux, Uy, pmg, "0_multigrid");
 write_output(Ux, Uy, p_sor, "0_sor");

 free(mesh);
 free(p.mesh);

 



 return 0;	
}

