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
#define Nx 300
#define Ny 300

int main()
{

 const double Lx = 1;
 const double Ly = 1;

 FVcell **mesh      = buildMesh(Nx, Ny, Lx, Ly);

 volumeField Ux     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Ux_star= buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy_star= buildFVfield(Nx, Ny, Lx, Ly);
 volumeField div    = buildFVfield(Nx, Ny, Lx, Ly);

 volumeField p      = buildFVfield(Nx, Ny, Lx, Ly);



 FILE *OutFile_meshx = fopen("x.dat","w++");
 
   for(int i=1; i<=Nx; i++){
       double x = (Uy.mesh[i][0].xFace(1) );
      fprintf(OutFile_meshx,"\n %f ", x); }  

fclose(OutFile_meshx);

 FILE *OutFile_meshy = fopen("y.dat","w++");
 
   for(int j=1; j<=Ny; j++){
       double y = (Uy.mesh[0][j].yFace(1) );
      fprintf(OutFile_meshy,"\n %f ", y); }  

fclose(OutFile_meshy);

 // initialize initial fiels
 initFieldVal(Ux, 0);
 initFieldVal(Uy, 0);
 initFieldVal(p , 0);

 initFieldVal(div , 0);

 volumeField  conv_u_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  conv_v_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  diff_u_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  diff_v_expl      = buildFVfield(Nx, Ny, Lx, Ly);


 double mu;
 mu = 0.001;
 double dt;
 double rho;
 rho = 1;
   double rho_f, rho_b;


 dt = 1e-4;
 double t0   = 0.;
 double tend = 100;
 double t;

 list<volumeField> Ux_time_series;
 list<volumeField> Uy_time_series;
 list<volumeField>  p_time_series;
  
 t = t0;

 setBCuFV(&Ux);
 setBCvFV(&Uy);

 double dt_write = 0.002;
 double t_last   = 0;
 int n_file = 0;
 n_iter_ssor = 0;
 int n_iter=0;

// write first initialized field 
 std::stringstream ss;
 ss << std::setw(10) << std::setfill('0') << n_file;
 std::string namefile = ss.str();

 while (t<=tend)
 {
  t = t + dt;
 
  cout << " -- time = "<<t<<" ---\n";
  cout << "\n";
 
  convection_u_expl (&conv_u_expl, Ux, Uy, rho);
  convection_v_expl (&conv_v_expl, Ux, Uy, rho);
 
  diffusion_u_expl (&diff_u_expl, Ux, Uy, mu);
  diffusion_v_expl (&diff_v_expl, Ux, Uy, mu); 

   predictor_u (&Ux_star, Ux, Uy,
               conv_u_expl, 
               diff_u_expl, rho, dt, mu, 1e-10, n_iter_u); 
  cout << "Ux iterations: "<< n_iter_u<<" iterations\n";
 
  predictor_v (&Uy_star, Ux, Uy,
               conv_v_expl, 
               diff_v_expl, rho, dt, mu, 1e-10, n_iter_v);
  cout << "Uy iterations: "<< n_iter_v<<" iterations\n"; 
  

  //// predictor step for velocity Ux and Uy
 // Ux_star =   Ux+(conv_u_expl*(-1./rho) + diff_u_expl*(1./rho) )*dt ;
  //Uy_star =   Uy+(conv_v_expl*(-1./rho) + diff_v_expl*(1./rho) )*dt ;
  cout<<"\n";

  setBCuFV(&Ux_star);
  setBCvFV(&Uy_star);

  ///////// compute divergence of predictor steps field /////////
  divergence(&div, Ux_star, Uy_star);
  cout << "mass imbalance "<< l1_norm(div)<<"\n";
  
  //////// solve Poisson Equation for Pressure /////////
  //if (n_iter_ssor >100 &&  n_iter>4000)
  //{
  // multigrid(&p, div, 1e-14, rho, dt, 1.3, Lx, Ly, 2);
  //}
  //else {
     poisson(&p, div, rho, 0, 1e-14, 0.7, n_iter_ssor);
  //  cout << "convergence achieved on pressure in "<< n_iter_ssor <<" iterations\n";//}

  /////// corrector step: projet the predicted velocity field in a div free space ////////
  correct (&Ux, &Uy, Ux_star, Uy_star, p, rho, dt); 
  setBCuFV(&Ux);
  setBCvFV(&Uy);

  divergence(&div, Ux, Uy);
  cout << "mass imbalance after correction "<< l1_norm(div)<<"\n";

  if (t - t_last > dt_write)
  {
    n_file ++;
    std::stringstream ss;
    ss << std::setw(10) << std::setfill('0') << n_file;
    std::string namefile = ss.str();

     write_output(Ux, Uy, p, namefile);
     t_last = t;
  
 }
  //Ux_time_series.push_back(Ux);
  //Uy_time_series.push_back(Uy);
  n_iter++;
 };

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


 free(mesh);
 free(Ux.mesh);
 free(Uy.mesh);
 free(p.mesh);

 



 return 0;	
}

