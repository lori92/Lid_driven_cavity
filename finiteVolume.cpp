#include <iostream>
#include "mesh.h"
using namespace std;

///////////////// finite volume discretization of explicit convective terms //////////
void convection_u_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, const double rho)
{
  double fUU_f, fUU_b;
  double fUV_f, fUV_b;
  double fUW_f, fUW_b;

  double rho_f, rho_b;
  
  int N_x = Ux.dimx;
  int N_y = Ux.dimy;
  int N_z = Ux.dimz;

  //// convective term d(uu)/dx +  d(uv)/dy +  d(uw)/dz -- u-control volume	 
  for (int i = 1; i <= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {  
    for(int k = 1; k<= N_z; k++)
    {  
     rho_f = rho;
     rho_b = rho;   
     fUU_f = 0.25 * rho_f * (Ux.mesh[i+1][j][k].cellVal() + Ux.mesh[ i ][j][k].cellVal()) * 
                            (Ux.mesh[i+1][j][k].cellVal() + Ux.mesh[ i ][j][k].cellVal());

     fUU_b = 0.25 * rho_b * (Ux.mesh[ i ][j][k].cellVal() + Ux.mesh[i-1][j][k].cellVal()) * 
                            (Ux.mesh[ i ][j][k].cellVal() + Ux.mesh[i-1][j][k].cellVal());

     rho_f = rho;
     rho_b = rho; 
     fUV_f = 0.5*rho_f*(Ux.mesh[i ][j+1][k].cellVal() + Ux.mesh[i][j][k].cellVal() ) * 
                  0.5* (Uy.mesh[i][ j ][k].cellVal() + Uy.mesh[i+1][ j ][k].cellVal() );

     fUV_b = 0.5*rho_b*(Ux.mesh[i ][j ][k].cellVal() + Ux.mesh[i ][j-1][k].cellVal() ) * 
                  0.5* (Uy.mesh[i][ j-1][k].cellVal() + Uy.mesh[i+1][j-1][k].cellVal() );


     rho_f = rho;
     rho_b = rho; 
     fUW_f = 0.5*rho_f*(Ux.mesh[i][j][k+1].cellVal() + Ux.mesh[i][j][k].cellVal() ) * 
                  0.5* (Uz.mesh[i+1][j][k].cellVal() + Uz.mesh[i][j][k].cellVal() );

     fUW_b = 0.5*rho_b*(Ux.mesh[i][j][k].cellVal()     + Ux.mesh[i][j][k-1].cellVal() ) * 
                  0.5* (Uz.mesh[i+1][j][k-1].cellVal() + Uz.mesh[i][j][k-1].cellVal() );

     conv->mesh[i][j][k].setVal(  (fUU_f - fUU_b) /  (Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[i][j][k].xCentroid())  
                                 +(fUV_f - fUV_b) /  (Ux.mesh[i][j][k].dy())
                                 +(fUW_f - fUW_b) /  (Ux.mesh[i][j][k].dz()) ); 
    };
   };
  };
  return;
};



void convection_v_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, const double rho)
{

  double fVU_f, fVU_b;
  double fVV_f, fVV_b;
  double fVW_f, fVW_b;

  double rho_f, rho_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;
  int N_z = Uy.dimz;

  //// convective term d(vu)/dx +  d(vv)/dy   +  d(vw)/dz -- v-control volume     
  for (int i = 1; i <= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {  
    for(int k = 1; k<= N_z; k++)
    { 
     rho_f = rho;
     rho_b = rho;
     fVU_f = 0.25 * rho_f * (Ux.mesh[ i ][j+1][k].cellVal() + Ux.mesh[ i ][j][k].cellVal()) * 
                            (Uy.mesh[i+1][ j ][k].cellVal() + Uy.mesh[ i ][j][k].cellVal()); 
     fVU_b = 0.25 * rho_b * (Ux.mesh[i-1][j+1][k].cellVal() + Ux.mesh[i-1][j][k].cellVal()) * 
                            (Uy.mesh[ i ][ j ][k].cellVal() + Uy.mesh[i-1][j][k].cellVal() );       
     
     rho_f = rho;
     rho_b = rho;
     fVV_f = 0.25 * rho_f * (Uy.mesh[i][j+1][k].cellVal() + Uy.mesh[i][ j ][k].cellVal()) *
                            (Uy.mesh[i][j+1][k].cellVal() + Uy.mesh[i][ j ][k].cellVal()); 
     fVV_b = 0.25 * rho_b * (Uy.mesh[i][ j ][k].cellVal() + Uy.mesh[i][j-1][k].cellVal()) *
                            (Uy.mesh[i][ j ][k].cellVal() + Uy.mesh[i][j-1][k].cellVal()); 

     rho_f = rho;
     rho_b = rho; 
     fVW_f = 0.5*rho_f*(Uy.mesh[i][j][k+1].cellVal() + Uy.mesh[i][j][k].cellVal() ) * 
                  0.5* (Uz.mesh[i][j+1][k].cellVal() + Uz.mesh[i][j][k].cellVal() );

     fVW_b = 0.5*rho_b*(Uy.mesh[i][j][k].cellVal()     + Uy.mesh[i][j][k-1].cellVal() ) * 
                  0.5* (Uz.mesh[i][j+1][k-1].cellVal() + Uz.mesh[i][j][k-1].cellVal() );


     conv->mesh[i][j][k].setVal(   (fVU_f - fVU_b) / (Uy.mesh[i][j][k].dx())  
                                 + (fVV_f - fVV_b) / (Uy.mesh[i][j][k].dy()) 
                                 + (fVW_f - fVW_b) / (Uy.mesh[i][j][k].dz()) ); 
      };
    };
  };
  return;
};




void convection_w_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, const double rho)
{

  double fWU_f, fWU_b;
  double fWV_f, fWV_b;
  double fWW_f, fWW_b;

  double rho_f, rho_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;
  int N_z = Uy.dimz;

  //// convective term d(wu)/dx +  d(wv)/dy   +  d(ww)/dz -- w-control volume     
  for (int i = 1; i <= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {  
    for(int k = 1; k<= N_z; k++)
    { 
     rho_f = rho;
     rho_b = rho;
     fWU_f = 0.25 * rho_f * (Uz.mesh[ i ][j][k+1].cellVal() + Uz.mesh[ i ][j][ k ].cellVal()) * 
                            (Ux.mesh[i][ j ][k+1].cellVal() + Ux.mesh[i-1][j][k+1].cellVal()); 
     fWU_b = 0.25 * rho_f * (Uz.mesh[ i ][j][k].cellVal() + Uz.mesh[ i ][j][k-1].cellVal()) * 
                            (Ux.mesh[i][ j ][k].cellVal() + Ux.mesh[i-1][j][ k ].cellVal());     
     
     rho_f = rho;
     rho_b = rho;
     fWV_f = 0.25 * rho_f * (Uz.mesh[i][j+1][k].cellVal() + Uz.mesh[i][ j ][k].cellVal()) *
                            (Uy.mesh[i][j][k+1].cellVal() + Uy.mesh[i][ j ][k].cellVal()); 
     fWV_b = 0.25 * rho_b * (Uz.mesh[i][ j ][k].cellVal() + Uz.mesh[i][j-1][k].cellVal()) *
                            (Uy.mesh[i][j-1][k+1].cellVal() + Uy.mesh[i][j-1][k].cellVal()); 

     rho_f = rho;
     rho_b = rho; 
     fWW_f = 0.25* rho_f* (Uz.mesh[i][j][k+1].cellVal() + Uz.mesh[i][j][k].cellVal() ) * 
                          (Uz.mesh[i][j][k+1].cellVal() + Uz.mesh[i][j][k].cellVal() );

     fWW_b = 0.25* rho_b* (Uz.mesh[i][j][ k ].cellVal() + Uz.mesh[i][j][k-1].cellVal() ) * 
                          (Uz.mesh[i][j][ k ].cellVal() + Uz.mesh[i][j][k-1].cellVal() );


     conv->mesh[i][j][k].setVal(   (fVU_f - fVU_b) / (Uz.mesh[i][j][k].dx())  
                                 + (fVV_f - fVV_b) / (Uz.mesh[i][j][k].dy()) 
                                 + (fWW_f - fWW_b) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid()) ); 
      };
    };
  };
  return;
};

///////////////// finite volume discretization of explicit diffusive terms //////////

void diffusion_u_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, const double mu)
{
  double ddU_dxx, dU_dx_f, dU_dx_b;
  double ddU_dyy, dU_dy_f, dU_dy_b;
  double ddU_dzz, dU_dz_f, dU_dz_b;

  int N_x = Uy.dimx;
  int N_y = Uy.dimy;
  int N_z = Uy.dimz;

  for (int i = 1; i<= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {
    for(int k = 1; k<= N_z; k++)
    {         
     // diffusive term dd(u)/dxx +  dd(u)/dyy + dd(u)/dzz-- u-control volume     
     dU_dx_f = mu * (  (Ux.mesh[i+1][j][k].cellVal()   - Ux.mesh[ i ][j][k].cellVal()) /
                       (Ux.mesh[i+1][j][k].dx()  ) );
     dU_dx_b = mu * (  (Ux.mesh[ i ][j][k].cellVal() - Ux.mesh[i-1][j][k].cellVal()) / 
                       (Ux.mesh[ i ][j][k].dx( ) ) );
     ddU_dxx = (dU_dx_f - dU_dx_b) / (Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[i][j][k].xCentroid());

     dU_dy_f = mu * (  (Ux.mesh[i][j+1][k].cellVal()     - Ux.mesh[i][ j ][k].cellVal()) / 
                       (Ux.mesh[ i ][j+1][k].yCentroid() - Ux.mesh[ i ][ j ][k].yCentroid()) );
     dU_dy_b = mu * (  (Ux.mesh[i][j  ][k].cellVal()     - Ux.mesh[i][j-1][k].cellVal()) / 
                       (Ux.mesh[ i ][ j ][k].yCentroid() - Ux.mesh[ i ][j-1][k].yCentroid()) );
     ddU_dyy = (dU_dy_f - dU_dy_b) / Ux.mesh[i][j][k].dy();       

     dU_dz_f = mu * (  (Ux.mesh[i][j][k+1].cellVal()   - Ux.mesh[i][j][k].cellVal()  ) / 
                       (Ux.mesh[i][j][k+1].dz()) );
     dU_dz_b = mu * (  (Ux.mesh[i][j][k].cellVal()     - Ux.mesh[i][j][k-1].cellVal()) / 
                       (Ux.mesh[i][j][k].dz())   );                       
     ddU_dzz = (dU_dz_f - dU_dz_b) / (Ux.mesh[i][j][k+1].zCentroid() - Ux.mesh[i][j][k].zCentroid());       

     diff->mesh[i][j][k].setVal( ddU_dxx + ddU_dyy + ddU_dzz);
    };
   };
  };
  return;
};


void diffusion_v_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, double mu)
{
  double ddV_dxx, dV_dx_f, dV_dx_b;
  double ddV_dyy, dV_dy_f, dV_dy_b;
  double ddV_dzz, dV_dz_f, dV_dz_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;
  int N_z = Uz.dimz;

  //// advective term dd(v)/dxx +  dd(v)/dyy  +  dd(w)/dzz -- v-control volume
  for (int i = 1; i <= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {  
    for(int k = 1; k<= N_z; k++)
    {  
     dV_dx_f = mu * (  (Uy.mesh[i+1][j][k].cellVal() - Uy.mesh[ i ][j][k].cellVal()) / 
                       (Uy.mesh[i+1][j][k].xCentroid() - Uy.mesh[ i ][j][k].xCentroid()));
     dV_dx_b = mu * (  (Uy.mesh[ i ][j][k].cellVal() - Uy.mesh[i-1][j][k].cellVal()) / 
                       (Uy.mesh[ i ][j][k].xCentroid() - Uy.mesh[i-1][j][k].xCentroid()));
     ddV_dxx = (dV_dx_f - dV_dx_b) / (Uy.mesh[i][j][k].dx());

     dV_dy_f = mu * (  (Uy.mesh[i][j+1][k].cellVal() - Uy.mesh[i][ j ][k].cellVal()) / (Uy.mesh[i][j+1][k].dy()));
     dV_dy_b = mu * (  (Uy.mesh[i][ j ][k].cellVal() - Uy.mesh[i][j-1][k].cellVal()) / (Uy.mesh[i][ j ][k].dy()));
     ddV_dyy = (dV_dy_f - dV_dy_b) / (Uy.mesh[i][j+1][k].yCentroid() - Uy.mesh[i][j][k].yCentroid());

     dV_dz_f = mu * (  (Uy.mesh[i][j][k+1].cellVal() - Uy.mesh[i][ j ][k].cellVal()) / (Uy.mesh[i][j][k+1].zCentroid() - Uy.mesh[i][j][k].zCentroid() ) );
     dV_dz_b = mu * (  (Uy.mesh[i][j][ k ].cellVal() - Uy.mesh[i][j][k-1].cellVal()) / (Uy.mesh[i][j][ k ].zCentroid() - Uy.mesh[i][j][k-1].zCentroid() ) );
     ddV_dzz = (dV_dz_f - dV_dz_b) / (Uy.mesh[i][j][k].dz());

     diff->mesh[i][j][k].setVal( ddV_dxx + ddV_dyy + ddW_dzz); 
     };
    };
  };
  return;
};


void diffusion_w_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, const volumeField& Uz, double mu)
{
  double ddW_dxx, dW_dx_f, dW_dx_b;
  double ddW_dyy, dW_dy_f, dW_dy_b;
  double ddW_dzz, dW_dz_f, dW_dz_b;
  
  int N_x = Uy.dimx;
  int N_y = Uy.dimy;
  int N_z = Uz.dimz;

  //// advective term dd(v)/dxx +  dd(v)/dyy  +  dd(w)/dzz -- w-control volume
  for (int i = 1; i <= N_x; i++)
  {
   for(int j = 1; j<= N_y; j++)
   {  
    for(int k = 1; k<= N_z; k++)
    {  
     dW_dx_f = mu * (  (Uz.mesh[i+1][j][k].cellVal() - Uz.mesh[ i ][j][k].cellVal()) / 
                       (Uz.mesh[i+1][j][k].xCentroid() - Uz.mesh[ i ][j][k].xCentroid()));
     dW_dx_b = mu * (  (Uz.mesh[ i ][j][k].cellVal() - Uz.mesh[i-1][j][k].cellVal()) / 
                       (Uz.mesh[ i ][j][k].xCentroid() - Uz.mesh[i-1][j][k].xCentroid()));
     ddW_dxx = (dW_dx_f - dW_dx_b) / (Uz.mesh[i][j][k].dx());

     dW_dy_f = mu * (  (Uz.mesh[i][j+1][k].cellVal() - Uz.mesh[i][ j ][k].cellVal()) / 
                       (Uz.mesh[i][j+1][k].yCentroid() - Uz.mesh[i][j][k].yCentroid())  );
     dW_dy_b = mu * (  (Uz.mesh[i][ j ][k].cellVal() - Uz.mesh[i][j-1][k].cellVal()) /         
                       (Uz.mesh[i][ j ][k].yCentroid() - Uz.mesh[i][j][k].zCentroid())  );
     ddW_dyy = (dW_dy_f - dW_dy_b) / (Uz.mesh[i][j+1][k].dy());

     dW_dz_f = mu * (  (Uz.mesh[i][j][k+1].cellVal() - Uz.mesh[i][ j ][k].cellVal()) / (Uz.mesh[i][j][k+1].dz() ) );
     dW_dz_b = mu * (  (Uz.mesh[i][j][ k ].cellVal() - Uz.mesh[i][j][k-1].cellVal()) / (Uz.mesh[i][j][ k ].dz() ) );
     ddW_dzz = (dW_dz_f - dW_dz_b) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid());

     diff->mesh[i][j][k].setVal( ddV_dxx + ddV_dyy + ddW_dzz); 
     };
    };
  };
  return;
};