#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "global_variables.h"

#include "mesh.h"
#include "poisson_solver.h"
#include "finiteVolume.h"
#include "numerics.h"
#include "predictor_step_new.h"
#include "write_output.h"
#include <chrono>
using namespace std;

/// number of internal cells (without ghost nodes) ///
#define Nx 64
#define Ny 64

int main()
{

 int restart;
 //
 
 FILE* fid;
 fid = fopen("input.in","r");
 if (fid)
 {
  fscanf(fid, "%d\n",  &restart);
  fscanf(fid, "%lf\n", &t0);
  fscanf(fid, "%lf\n", &dt);
  fscanf(fid, "%lf\n", &tend);
  fscanf(fid, "%lf\n", &dt_write);

  fscanf(fid, "%lf\n", &mu);
  fscanf(fid, "%lf\n", &rho);
  fscanf(fid, "%lf\n", &Lx);
  fscanf(fid, "%lf\n", &Ly);  
 }
 else
 {
   cout << "***** input.in not found *****\n";
   return 4;
 };

 fclose(fid);


 volumeField Ux     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Ux_star= buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy_star= buildFVfield(Nx, Ny, Lx, Ly);
 volumeField div    = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField p      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField b      = buildFVfield(Nx, Ny, Lx, Ly);

 // get the dsicretization coefficients:
 // lap(p)_(i,j) = aE*p_(i,j) + aW*p_(i-1,j) +aP*p_(i,j) +aN*p_(i,j) +aS*p_(i,j) 
 
 double** aP; //[Nx][Ny];
 double** aS; //[Nx][Ny];
 double** aN; //[Nx][Ny];
 double** aW; //[Nx][Ny];
 double** aE; //[Nx][Ny];

 aP = new double* [Nx+2];
 aS = new double* [Nx+2];
 aN = new double* [Nx+2];
 aW = new double* [Nx+2];
 aE = new double* [Nx+2];

 for (int i = 0; i<=Nx+1; i++)
 {
  aP[i] = new double [Ny+2];
  aS[i] = new double [Ny+2];
  aN[i] = new double [Ny+2];
  aW[i] = new double [Ny+2];
  aE[i] = new double [Ny+2];
 }
 std:: cout << "assemebling discretization coefficients .."<<"\n";
 assemble_coeff(p, aP, aN, aS, aW, aE);
 std:: cout << "done"<<"\n";

 // get the discretized poisson operator L: Lp = b :
 int N = Nx*Ny;
 int dim = N+1;

 double** L; 
 double* L_reshaped;
 double* L_reshaped_tmp;

 L          = new double* [N+1];
 for (int n = 0; n<N+1; n++)
 {
  L[n] = new double [N+1];
 }

 L_reshaped     = new double  [(N+1)*(N+1)];
 L_reshaped_tmp = new double  [(N+1)*(N+1)];
 
 double L_array [(Nx*Ny+1)*(Nx*Ny+1)];
 double L_array_tmp [(Nx*Ny+1)*(Nx*Ny+1)];
 double b_array [Nx*Ny+1];
 

  // stuff for direct solver 

 int nrhs = 1;
 int lda = N+1;
 int ipiv[N+1];
 int ldb = N+1;
 int info = 0;
 int factorization=0;

 std:: cout << "assembling discretized poisson operator L: Lp = b .."<<"\n";
 assemble_matrix(L, Nx,Ny,aP,aS,aN,aW,aE,p);
 reshape(N+1,  &L_array[0], L);
  std:: cout << "done"<<"\n";

 /*
 std:: cout << "compute LU factorization of L: .."<<"\n";
 dgetrf_(&dim, &dim, &L_array[0], &lda, ipiv, &info);
 factorization=1;

 for (int i=0;i<N+1;i++)
 {
   for (int j=0;j<N+1;j++){
  cout<<L[i][j]<<" ";

   }
     cout<<"\n";

 }


  for (int i=0;i<(N+1)*(N+1);i++)
 {
     cout<<L_array[i]<<" ";
     cout<<"\n";
 }
 */
 // allocate the memory for the forcing term (div(u)) b :
 double* b_reshaped;
 b_reshaped          = new double [N+1];
 

 cout << "------- Lid Cavity Flow 2D ---------"<<"\n";
 cout << "\n";
 cout << "simulation parameters:"<<"\n";
 cout << "\n";

 switch (restart)
 {
   case -2:
   {
    std:: cout << "discrete Laplacian operator successfully written... stop here"<<"\n";
    return 0;
   }
   case -1:
   {
    std:: cout << "check on Poisson Solver with assigned forcing term.."<<"\n";
    initFieldVal(b, 0);
    b_reshaped = forcing_term(b,aP,aS,aN,aW,aE);

    write_output(Ux, Uy, b, t0, "0_CHECK_FORCING_TERM");

    int dim = N+1;

    dgesv_(&dim, &nrhs, &L_array[0], &lda, ipiv, &b_array[0], &ldb, &info);

    getVolumeField(&p,  &b_array[0]);

    write_output(Ux, Uy, p, t0, "0_poisson_DIRECT_SOLVER");

    initFieldVal(p, 0);

    return 0;
   }

   case 0:
   {
    read_data( &Ux,  &Uy,   &p,  t0);
    std:: cout << "starting from previous solution at time = "<<t0<<" s"<<"\n";
    break;
   }

   case 1:
   {
    std:: cout << "interpolating previous solution at time = "<<t0<<" s on the current grid"<<"\n";
    // interpolate solution from a coarser grid 
    volumeField  Ux_coarse = buildFVfield(Nx/2, Ny/2, Lx, Ly);
    volumeField  Uy_coarse = buildFVfield(Nx/2, Ny/2, Lx, Ly);
    volumeField  p_coarse  = buildFVfield(Nx/2, Ny/2, Lx, Ly);

    read_data( &Ux_coarse,  &Uy_coarse,   &p_coarse,  t0);

    prolongateFieldVal(Ux_coarse , &Ux); 
    prolongateFieldVal(Uy_coarse , &Uy); 
    prolongateFieldVal(p_coarse  , &p);
    cout<< "interpolation performed...\n";
    break;
   }
   case 2:
   {
    t0 = 0;
    std:: cout << "starting new computation from time = "<<t0<<" s"<<"\n";
    initFieldVal(Ux, 0);
    initFieldVal(Uy, 0);
    initFieldVal(p , 0);
        break;
   }
 }

 initFieldVal(div , 0);
 initFieldVal(b, 0);

 cout << "-----------------------------------"<<"\n";
 cout << "tend "<<tend <<"\n";
 cout << "dt "<<dt <<"\n";
 cout << "mu "<<mu <<"\n";
 cout << "rho "<<rho <<"\n";
 cout << "Lx "<<Lx <<"\n";
 cout << "Ly "<<Ly <<"\n";
 cout << "-----------------------------------"<<"\n";
 cout << "\n";

 volumeField  conv_u_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  conv_v_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  diff_u_expl      = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  diff_v_expl      = buildFVfield(Nx, Ny, Lx, Ly);

 list<volumeField> Ux_time_series;
 list<volumeField> Uy_time_series;
 list<volumeField>  p_time_series;
  
 t = t0;

 setBCuFV(&Ux);
 setBCvFV(&Uy);

 double t_last   = 0;
 int n_file = 0;
 n_iter_ssor = 0;
 int n_iter=0;

// write first initialized field 
 std::stringstream ss;
 ss << std::setw(10) << std::setfill('0') << n_file;
 std::string namefile = ss.str();

 double rho_f, rho_b;
 time_t start, end;

 time(&start);

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
               diff_u_expl, rho, dt, mu, 1e-12, n_iter_u); 
  cout << "Ux iterations: "<< n_iter_u<<" iterations\n";

  predictor_v (&Uy_star, Ux, Uy,
               conv_v_expl, 
               diff_v_expl, rho, dt, mu, 1e-12, n_iter_v);
  cout << "Uy iterations: "<< n_iter_v<<" iterations\n"; 
 

  //// predictor step for velocity Ux and Uy
  //Ux_star =   Ux+(conv_u_expl*(-1./rho) + diff_u_expl*(1./rho) )*dt ;
  //Uy_star =   Uy+(conv_v_expl*(-1./rho) + diff_v_expl*(1./rho) )*dt ;
  cout<<"\n";

  setBCuFV(&Ux_star);
  setBCvFV(&Uy_star);

  ///////// compute divergence of predictor steps field /////////
  divergence(&div, Ux_star, Uy_star);
  cout << "mass imbalance "<< l1_norm(div)<<"\n";
  
  //////// solve Poisson Equation for Pressure /////////
  solvePressure_test(  &L_array[0], &b_array[0], &p, div, dt, rho, factorization, info, ipiv);
  //solvePressure( L, &L_array[0],  &p, div, dt, rho, aP, aS, aN, aW, aE);
 
  //L_array = L_array_tmp;

  /////// corrector step: projet the predicted velocity field in a div free space ////////
  correct (&Ux, &Uy, Ux_star, Uy_star, p, rho, dt); 
  setBCuFV(&Ux);
  setBCvFV(&Uy);

  divergence(&div, Ux, Uy);
  cout << "mass imbalance after correction "<< l1_norm(div)<<"\n";

  if (t - t_last >= dt_write)
  {
    n_file ++;
    std::stringstream ss;
    ss << std::setw(10) << std::setfill('0') << n_file;
    std::string namefile = ss.str();

     write_output(Ux, Uy, p, t, namefile);
     write_data (Ux, Uy, p, t);
     write_centre_line(Ux, Uy, t);
     t_last = t;  
 }
  //Ux_time_series.push_back(Ux);
  //Uy_time_series.push_back(Uy);
  n_iter++;

 write_k(Ux, Uy, t); 

 };

     // interpolate on the coarser grid //
 volumeField  Ux_coarser = buildFVfield(Nx/2, Ny/2, Lx, Ly);
 volumeField  Uy_coarser = buildFVfield(Nx/2, Ny/2, Lx, Ly);
 volumeField  p_coarser  = buildFVfield(Nx/2, Ny/2, Lx, Ly);

 interpolateFieldVal(Ux , &Ux_coarser); 
 interpolateFieldVal(Uy , &Uy_coarser); 
 interpolateFieldVal(p  , &p_coarser); 
 write_output(Ux_coarser, Uy_coarser, p_coarser, tend, "0_coarse_grid");


      // prolongate on the coarser grid //
 volumeField  Ux_prl = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  Uy_prl = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField  p_prl  = buildFVfield(Nx, Ny, Lx, Ly);

 prolongateFieldVal(Ux_coarser , &Ux_prl); 
 prolongateFieldVal(Uy_coarser , &Uy_prl); 
 prolongateFieldVal(p_coarser , &p_prl); 
 write_output(Ux_prl, Uy_prl, p_prl, tend, "0_fine_grid");


 free(Ux.mesh);
 free(Uy.mesh);
 free(p.mesh);

 delete (b_reshaped);
 delete (L);

 
  time(&end);
  double time_taken = double(end - start);
  cout << "Execution time is : " << fixed << time_taken << setprecision(5);
  cout << " sec " << endl;  


 return 0;	
}
