#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include "mesh.cpp"
using namespace std;

/// number of cells (without ghost nodes) ///
#define Nx 5
#define Ny 5



int main()
{

 const double Lx = 1;
 const double Ly = 1;

 FVcell **mesh = buildMesh(Nx, Ny, Lx, Ly);
 FVcell **Ux   = buildMesh(Nx, Ny, Lx, Ly);
 FVcell **Uy   = buildMesh(Nx, Ny, Lx, Ly);
 FVcell **Ux_star   = buildMesh(Nx, Ny, Lx, Ly);
 FVcell **Uy_star   = buildMesh(Nx, Ny, Lx, Ly);

 volumeField p    = buildFVfield(Nx, Ny, Lx, Ly);

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

double fVU_f, fVU_b;
double fVV_f, fVV_b;

double ddU_dxx, dU_dx_f, dU_dx_b;
double ddU_dyy, dU_dy_f, dU_dy_b;

double ddV_dxx, dV_dx_f, dV_dx_b;
double ddV_dyy, dV_dy_f, dV_dy_b;

double conv_u, conv_v;
double diff_u, diff_v;

double mu;
double dt;
double rho;

dt = 0.001;
double t0   = 0.;
double tend = 0.001;
double t;

list<FVcell**> Ux_time_series; 	
list<FVcell**> Uy_time_series; 	
list<FVcell**> p_time_series; 	

t = t0;

while (t<=tend)
{
 t = t + dt;

 for (int i = 1; i <= Nx; i++)
 {
   for(int j = 1; j<= Ny; j++)
   {

     //// convective term d(uu)/dx +  d(uv)/dy -- u-control volume	   
     fUU_f = 0.25 * (Ux[i][j].cellVal() + Ux[i][j].cellVal()) *(Ux[i+1][j].cellVal() + Ux[i+1][j].cellVal());
     fUU_b = 0.25 * (Ux[i][j].cellVal() + Ux[i][j].cellVal()) *(Ux[i-1][j].cellVal() + Ux[i-1][j].cellVal());

     fUV_f = 0.5*rho*(Ux[i ][j+1].cellVal() + Ux[i][j].cellVal() ) * \
                                  0.5* (  Uy[i][ j ].cellVal() + Uy[i+1][ j ].cellVal() );
     fUV_b = 0.5*rho*(Ux[i ][j ].cellVal() + Ux[i-1][j-1].cellVal() ) * \
                                  0.5* (  Uy[i][ j ].cellVal() + Uy[i+1][j-1].cellVal() );

     conv_u = (fUU_f - fUU_b) / (mesh[i][j].dx())  \
             +(fUV_f - fUV_b) / (mesh[i][j].dy()); 

     //// convective term d(vu)/dx +  d(vv)/dy -- v-control volume     
     fVU_f = 0.5 * rho*(Ux[i][j+1].cellVal() + Ux[i] [ j].cellVal()) * 0.5 * (Uy[i+1][ j ].cellVal() + Uy[i][ j ].cellVal() ); 
     fVU_b = 0.5 * rho*(Ux[i][ j ].cellVal() + Ux[i][j-1].cellVal()) * 0.5 * (Uy[i+1][j-1].cellVal() + Uy[i][j-1].cellVal() );       


     conv_v = (fVU_f - fVU_b) / (mesh[i][j].dx())  \
             +(fVV_f - fVV_b) / (mesh[i][j].dy());
    
     //// advective term dd(u)/dxx +  dd(u)/dyy -- u-control volume     
     dU_dx_f = mu * (  (Ux[i+1][j].cellVal() - Ux[ i ][j].cellVal()) / (mesh[i+1][j].xFace(1) - mesh[i+1][j].xFace(0)));
     dU_dx_b = mu * (  (Ux[ i ][j].cellVal() - Ux[i-1][j].cellVal()) / (mesh[ i ][j].xFace(1) - mesh[ i ][j].xFace(0)));
     ddU_dxx = (dU_dx_f - dU_dx_b) / (mesh[i+1][j].xCentroid() - mesh[i][j].xCentroid());


     dU_dy_f = mu * ( 0.5 * (Ux[i][j+1].cellVal() - Ux[i][ j ].cellVal()) / (mesh[i+1][j+1].yCentroid() - mesh[i+1][ j ].yCentroid()) );
     dU_dy_b = mu * ( 0.5 * (Ux[i][j  ].cellVal() - Ux[i][j-1].cellVal()) / (mesh[ i ][ j ].yCentroid() - mesh[ i ][j-1].yCentroid()) );
     ddU_dyy = (dU_dy_f - dU_dy_b) / mesh[i][j].dy();       

     diff_u = ddU_dxx + ddU_dyy;


     //// advective term dd(v)/dxx +  dd(v)/dyy -- v-control volume
     dV_dx_f = mu * (  (Uy[i+1][j].cellVal() - Uy[ i ][j].cellVal()) / (mesh[i+1][j].xCentroid() - mesh[ i ][j].xCentroid()));
     dV_dx_b = mu * (  (Uy[ i ][j].cellVal() - Uy[i-1][j].cellVal()) / (mesh[ i ][j].xCentroid() - mesh[i-1][j].xCentroid()));
     ddV_dxx = (dV_dx_f - dV_dx_b) / (mesh[i][j].dx());

     dV_dy_f = mu * (  (Uy[i][j+1].cellVal() - Uy[i][ j ].cellVal()) / (mesh[i][j+1].dy()));
     dV_dy_b = mu * (  (Uy[i][ j ].cellVal() - Uy[i][j-1].cellVal()) / (mesh[i][ j ].dy()));
     ddV_dyy = (dV_dy_f - dV_dy_b) / (mesh[i][j+1].yCentroid() - mesh[i][j].yCentroid());

     diff_v = ddV_dxx + ddV_dyy;

     Ux_star[i][j].setVal( Ux[i][j].cellVal() - dt * (conv_u + diff_u) );
     Uy_star[i][j].setVal( Uy[i][j].cellVal() - dt * (conv_v + diff_v) );
   }
   cout<<"\n";
 } 
 setBCu(Ux_star, Nx, Ny);
 setBCv(Uy_star, Nx, Ny);

Ux_time_series.push_back(Ux_star);

};

/*
 for (int i = Nx+1; i <= Nx+1; i++)
 {
   for(int j = 0; j<= Ny+1; j++)
   {
    Uy_star[i][j].print_cellVal();
   }
   cout<<"\n";
 } 
*/


 free(mesh);
 free(Ux);
 free(Uy);
 return 0;	
}










///////////////////////////////////////////////////////////////////////////////////////////
double l1_norm(FVcell** const field, int N_x, int N_y) 
{
    double accum = 0.;

    for (int i = 1; i <= N_x; ++i) 
    {
         for (int j= 1; j <= N_y; ++i) 
         {
          accum += abs( field[i][j].cellVal() ) ;
         };
    };
    return sqrt(accum);
};
/////////////////////////////////////////////////////////////////////////////////////////////
void poisson (volumeField* p, FVcell** Ux_pred, FVcell** Uy_pred, double tol, double omega)
{
    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y;
    N_x = p->dimx;
    N_y = p->dimy;

    FVcell **p_ptr;

    // initialize mesh pointer to input pointer to pressure
    p_ptr = p->mesh; 

    int niter = 0;

    while (ratio > tol)
    {
        p_ptr = p->mesh; 
        l1norm_old =  l1_norm (p_ptr, N_x, N_y);
        niter+= 1;

        for (int i = 1; i <= N_x; ++i) {
            for (int j= 1; j <= N_y; ++i) {
                //p_ptr[i][j].setVal( (1-omega)*(p[i][j].cellVal() + omega * (-a2_x_b[i]*p[i-1][j].cellVal() - a2_x_f[i]*p[i+1][j].cellVal()
                //           -a2_y_b[j]*pn[i,j-1] - a2_y_f[j]*pn[i,j+1]  ) / (a2_x_c[i] + a2_y_c[j]) );
            };
        };
        setBCp(p_ptr, Nx, Ny);
        l1norm =  l1_norm (p_ptr, N_x, N_y);

        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old + 1e-10);
    };
return;
};


