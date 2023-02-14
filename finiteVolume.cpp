#include <iostream>
#include "mesh.h"
using namespace std;

///////////////// finite volume discretization of explicit convective terms //////////
void convection_u_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, const double rho)
{
  double fUU_f, fUU_b;
  double fUV_f, fUV_b;
  double rho_f, rho_b;
  
  int N_x = Ux.dimx;
  int N_y = Ux.dimy;
  int N_z = Ux.dimz;

  //// convective term d(uu)/dx +  d(uv)/dy -- u-control volume	 
  for (int i = 1; i <= N_x; i++)
  {
    for(int j = 1; j<= N_y; j++)
    {  
     rho_f = rho;
     rho_b = rho;   
     fUU_f = 0.25 * rho_f * (Ux.mesh[i+1][j].cellVal() + Ux.mesh[ i ][j].cellVal()) * 
                            (Ux.mesh[i+1][j].cellVal() + Ux.mesh[ i ][j].cellVal());

     fUU_b = 0.25 * rho_b * (Ux.mesh[ i ][j].cellVal() + Ux.mesh[i-1][j].cellVal()) * 
                            (Ux.mesh[ i ][j].cellVal() + Ux.mesh[i-1][j].cellVal());

     rho_f = rho;
     rho_b = rho; 
     fUV_f = 0.5*rho_f*(Ux.mesh[i ][j+1].cellVal() + Ux.mesh[i][j].cellVal() ) * 
                  0.5* (Uy.mesh[i][ j ].cellVal() + Uy.mesh[i+1][ j ].cellVal() );

     fUV_b = 0.5*rho_b*(Ux.mesh[i ][j ].cellVal() + Ux.mesh[i ][j-1].cellVal() ) * 
                  0.5* (Uy.mesh[i][ j-1].cellVal() + Uy.mesh[i+1][j-1].cellVal() );

     conv->mesh[i][j].setVal(  (fUU_f - fUU_b) /  (Ux.mesh[i+1][j].xCentroid() - Ux.mesh[i][j].xCentroid())  
                               +(fUV_f - fUV_b) / (Ux.mesh[i][j].dy()) ); 
    };
  };

  return;
};



void convection_v_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const double rho)
{

  double fVU_f, fVU_b;
  double fVV_f, fVV_b;
  double rho_f, rho_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;

  //// convective term d(vu)/dx +  d(vv)/dy -- v-control volume     
  for (int i = 1; i <= N_x; i++)
  {
    for(int j = 1; j<= N_y; j++)
    {  
     rho_f = rho;
     rho_b = rho;
     fVU_f = 0.25 * rho_f * (Ux.mesh[ i ][j+1].cellVal() + Ux.mesh[ i ][j].cellVal()) * 
                            (Uy.mesh[i+1][ j ].cellVal() + Uy.mesh[ i ][j].cellVal()); 
     fVU_b = 0.25 * rho_b * (Ux.mesh[i-1][j+1].cellVal() + Ux.mesh[i-1][j].cellVal()) * 
                            (Uy.mesh[ i ][ j ].cellVal() + Uy.mesh[i-1][j].cellVal() );       
     
     rho_f = rho;
     rho_b = rho;
     fVV_f = 0.25 * rho_f * (Uy.mesh[i][j+1].cellVal() + Uy.mesh[i][ j ].cellVal()) *
                            (Uy.mesh[i][j+1].cellVal() + Uy.mesh[i][ j ].cellVal()); 
     fVV_b = 0.25 * rho_b * (Uy.mesh[i][ j ].cellVal() + Uy.mesh[i][j-1].cellVal()) *
                            (Uy.mesh[i][ j ].cellVal() + Uy.mesh[i][j-1].cellVal()); 

     conv->mesh[i][j].setVal(  (fVU_f - fVU_b) / (Uy.mesh[i][j].dx())  
                                +(fVV_f - fVV_b) / (Uy.mesh[i][j].dy()) ); 
    };
  };
  return;
};

///////////////// finite volume discretization of explicit diffusive terms //////////

void diffusion_u_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, const double mu)
{
  double ddU_dxx, dU_dx_f, dU_dx_b;
  double ddU_dyy, dU_dy_f, dU_dy_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;

  for (int i = 1; i <= N_x; i++)
  {
    for(int j = 1; j<= N_y; j++)
    {    
     // diffusive term dd(u)/dxx +  dd(u)/dyy -- u-control volume     
     dU_dx_f = mu * (  (Ux.mesh[i+1][j].cellVal()   - Ux.mesh[ i ][j].cellVal()) /
                       (Ux.mesh[i+1][j].dx()  ) );
     dU_dx_b = mu * (  (Ux.mesh[ i ][j].cellVal() - Ux.mesh[i-1][j].cellVal()) / 
                       (Ux.mesh[ i ][j].dx( ) ) );
     ddU_dxx = (dU_dx_f - dU_dx_b) / (Ux.mesh[i+1][j].xCentroid() - Ux.mesh[i][j].xCentroid());

     dU_dy_f = mu * (  (Ux.mesh[i][j+1].cellVal()     - Ux.mesh[i][ j ].cellVal()) / 
                       (Ux.mesh[ i ][j+1].yCentroid() - Ux.mesh[ i ][ j ].yCentroid()) );
     dU_dy_b = mu * (  (Ux.mesh[i][j  ].cellVal()     - Ux.mesh[i][j-1].cellVal()) / 
                       (Ux.mesh[ i ][ j ].yCentroid() - Ux.mesh[ i ][j-1].yCentroid()) );
     ddU_dyy = (dU_dy_f - dU_dy_b) / Ux.mesh[i][j].dy();       

     diff->mesh[i][j].setVal( ddU_dxx + ddU_dyy );
    };
  };
  return;
};


void diffusion_v_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, double mu)
{
  double ddV_dxx, dV_dx_f, dV_dx_b;
  double ddV_dyy, dV_dy_f, dV_dy_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;

  //// advective term dd(v)/dxx +  dd(v)/dyy -- v-control volume
  for (int i = 1; i <= N_x; i++)
  {
    for(int j = 1; j<= N_y; j++)
    {  
     dV_dx_f = mu * (  (Uy.mesh[i+1][j].cellVal() - Uy.mesh[ i ][j].cellVal()) / 
                       (Uy.mesh[i+1][j].xCentroid() - Uy.mesh[ i ][j].xCentroid()));
     dV_dx_b = mu * (  (Uy.mesh[ i ][j].cellVal() - Uy.mesh[i-1][j].cellVal()) / (Uy.mesh[ i ][j].xCentroid() - Uy.mesh[i-1][j].xCentroid()));
     ddV_dxx = (dV_dx_f - dV_dx_b) / (Uy.mesh[i][j].dx());

     dV_dy_f = mu * (  (Uy.mesh[i][j+1].cellVal() - Uy.mesh[i][ j ].cellVal()) / (Uy.mesh[i][j+1].dy()));
     dV_dy_b = mu * (  (Uy.mesh[i][ j ].cellVal() - Uy.mesh[i][j-1].cellVal()) / (Uy.mesh[i][ j ].dy()));
     ddV_dyy = (dV_dy_f - dV_dy_b) / (Uy.mesh[i][j+1].yCentroid() - Uy.mesh[i][j].yCentroid());

     diff->mesh[i][j].setVal( ddV_dxx + ddV_dyy ); 
    };
  };
  return;
};