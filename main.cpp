#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "global_variables.h"

#include "mesh.cpp"
#include "poisson.cpp"
#include "finiteVolume.cpp"
#include "predictor_step.cpp"

#include "write_output.C"

using namespace std;

/// number of internal cells (without ghost nodes) ///
#define Nx 76
#define Ny 76
#define Nz 76

int main()
{

 const double Lx = 1;
 const double Ly = 1;
 const double Lz = 1;

 FVcell ***mesh      = buildMesh(Nx, Ny, Nz, Lx, Ly, Lz);

 volumeField Ux     = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField Uy     = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField Uz     = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);

 volumeField Ux_star= buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField Uy_star= buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField Uz_star= buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField div    = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);

 volumeField p      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);



 FILE *OutFile_meshx = fopen("x.dat","w++");
 
   for(int i=1; i<=Nx; i++){
       double x = (Uy.mesh[i][0][0].xFace(1) );
      fprintf(OutFile_meshx,"\n %f ", x); }  

fclose(OutFile_meshx);

 FILE *OutFile_meshy = fopen("y.dat","w++");
 
   for(int j=1; j<=Ny; j++){
       double y = (Uy.mesh[0][j][0].yFace(1) );
      fprintf(OutFile_meshy,"\n %f ", y); }  

fclose(OutFile_meshy);

 FILE *OutFile_meshz = fopen("z.dat","w++");
 
   for(int k=1; k<=Nz; k++){
       double z = (Uz.mesh[0][0][k].zFace(1) );
      fprintf(OutFile_meshz,"\n %f ", z); }  

fclose(OutFile_meshz);

 // initialize initial fiels
 initFieldVal(Ux, 0);
 initFieldVal(Uy, 0);
 initFieldVal(Uz, 0);
 initFieldVal(p , 0);

 initFieldVal(div , 0);

 volumeField  conv_u_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  conv_v_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  conv_w_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);

 volumeField  diff_u_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  diff_v_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  diff_w_expl      = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);


 double mu;
 mu = 0.001;
 double dt;
 double rho;
 rho = 1;
 double rho_f, rho_b;


 dt = 0.5e-3;
 double t0   = 0.;
 double tend = 100;
 double t;

 list<volumeField> Ux_time_series;
 list<volumeField> Uy_time_series;
 list<volumeField> Uz_time_series;

 list<volumeField>  p_time_series;
  
 t = t0;

 setBCuFV(&Ux);
 setBCvFV(&Uy);
 setBCvFV(&Uz);

 double dt_write = 0.1;
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
 
  convection_u_expl (&conv_u_expl, Ux, Uy, Uz, rho);
  convection_v_expl (&conv_v_expl, Ux, Uy, Uz, rho);
  convection_w_expl (&conv_w_expl, Ux, Uy, Uz, rho);
 
  diffusion_u_expl (&diff_u_expl, Ux, Uy, Uz, mu);
  diffusion_v_expl (&diff_v_expl, Ux, Uy, Uz, mu); 
  diffusion_w_expl (&diff_w_expl, Ux, Uy, Uz, mu); 

   predictor_u (&Ux_star, Ux, Uy, Uz, 
               conv_u_expl, 
               diff_u_expl, 
               rho, dt, mu, 1e-10, n_iter_u); 
  cout << "Ux iterations: "<< n_iter_u<<" iterations\n";
 
  predictor_v (&Uy_star, Ux, Uy, Uz, 
               conv_v_expl, 
               diff_v_expl, rho, dt, mu, 1e-10, n_iter_v);
  cout << "Uy iterations: "<< n_iter_v<<" iterations\n"; 

  predictor_w (&Uz_star, Ux, Uy, Uz, 
               conv_w_expl, 
               diff_w_expl, rho, dt, mu, 1e-10, n_iter_w);
  cout << "Uz iterations: "<< n_iter_w<<" iterations\n"; 
  

  cout<<"\n";

  setBCuFV(&Ux_star);
  setBCvFV(&Uy_star);
  setBCvFV(&Uz_star);

  ///////// compute divergence of predictor steps field /////////
  divergence(&div, Ux_star, Uy_star, Uz_star);
  cout << "mass imbalance "<< l1_norm(div)<<"\n";
  
  //////// solve Poisson Equation for Pressure /////////
  //if (n_iter_ssor >100 &&  n_iter>4000)
  //{
 //   multigrid(&p, div, 1e-14, rho, dt, 1.3, Lx, Ly, 2);
  //}
  //else {
     poisson(&p, div, rho, dt, 1e-14, 0.7, n_iter_ssor);
    cout << "convergence achieved on pressure in "<< n_iter_ssor <<" iterations\n";//}

  /////// corrector step: projet the predicted velocity field in a div free space ////////
  correct (&Ux, &Uy, &Uz, Ux_star, Uy_star, Uz_star, p, rho, dt); 
  setBCuFV(&Ux);
  setBCvFV(&Uy);
  setBCvFV(&Uz);


  divergence(&div, Ux, Uy, Uz);
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
 volumeField  Ux_coarser = buildFVfield(Nx/2, Ny/2, Nz/2, Lx, Ly, Lz);
 volumeField  Uy_coarser = buildFVfield(Nx/2, Ny/2, Nz/2, Lx, Ly, Lz);
 volumeField  Uz_coarser = buildFVfield(Nx/2, Ny/2, Nz/2, Lx, Ly, Lz);

 volumeField  p_coarser  = buildFVfield(Nx/2, Ny/2, Nz/2, Lx, Ly, Lz);

 interpolateFieldVal(Ux , &Ux_coarser); 
 interpolateFieldVal(Uy , &Uy_coarser); 
 interpolateFieldVal(Uz , &Uz_coarser); 

 interpolateFieldVal(p  , &p_coarser); 
 write_output(Ux_coarser, Uy_coarser, Uz_coarser, p_coarser, "0_coarse_grid");


      // prolongate on the coarser grid //
 volumeField  Ux_prl = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  Uy_prl = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);
 volumeField  Uz_prl = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);

 volumeField  p_prl  = buildFVfield(Nx, Ny, Nz, Lx, Ly, Lz);

 prolongateFieldVal(Ux_coarser , &Ux_prl); 
 prolongateFieldVal(Uy_coarser , &Uy_prl); 
 prolongateFieldVal(Uz_coarser , &Uz_prl); 

 prolongateFieldVal(p_coarser , &p_prl); 
 write_output(Ux_prl, Uy_prl, Uz_prl, p_prl, "0_fine_grid");


 free(mesh);
 free(Ux.mesh);
 free(Uy.mesh);
 free(Uz.mesh);

 free(p.mesh);

 



 return 0;	
}

