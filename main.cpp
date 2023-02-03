#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include "mesh.cpp"
using namespace std;

/// number of cells (without ghost nodes) ///
#define Nx 50
#define Ny 50


///////////////////////////////////////////////////////////////////////////////////////////
double l1_norm(volumeField const field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;

    double accum = 0.;

    for (int i = 1; i <= N_x; ++i)
    {
         for (int j= 1; j <= N_y; ++j)
         {
          accum += abs( field.mesh[i][j].cellVal() ) ;
         };
    };
    return sqrt(accum);
};
/////////////////////////////////////////////////////////////////////////////////////////////
void poisson (volumeField* p, volumeField const Ux_pred, volumeField const Uy_pred, double tol, double omega)
{
    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y;
    N_x = p->dimx;
    N_y = p->dimy;

    // east and west coefficients for discretized Poisson equation
    double aE, aW, aP;
    double aN, aS;

    // initialize mesh pointer to input pointer to pressure
    volumeField p_star;
    p_star = *p;
 
    setBCpFV(&p_star);
    int niter = 0;

//    cout << "start poisson iterations....\n";

    while (ratio > tol)
    {
        p_star.mesh = p->mesh;
  //      cout << "update old pressure....\n";

        l1norm_old =  l1_norm (p_star);
        cout << "computing l1- norm....\n";

        niter+= 1;

        for (int i = 1; i <= N_x; ++i) {
            for (int j= 1; j <= N_y; ++j) {

                aP = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() )
                     - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );

                aW = - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() );
                aE = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() );

                aS = - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );
                aN = - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() );

                p->mesh[i][j].setVal( (1.-omega)*p_star.mesh[i][j].cellVal() + omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal())   );


            };
        };
        setBCpFV(p);
        l1norm =  l1_norm (*p);
        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old + 1e-10);
    };
    cout << "convergence achieved in "<< niter<<" iterations\n";
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{

 const double Lx = 50;
 const double Ly = 50;

 FVcell **mesh      = buildMesh(Nx, Ny, Lx, Ly);

 volumeField Ux     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy     = buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Ux_star= buildFVfield(Nx, Ny, Lx, Ly);
 volumeField Uy_star= buildFVfield(Nx, Ny, Lx, Ly);

 volumeField p      = buildFVfield(Nx, Ny, Lx, Ly);

 // initialize initial fiels
 initFieldVal(Ux, 0);
 initFieldVal(Uy, 0);
 initFieldVal(p , 0);

 // impose boundary conditions on initial fields
 setBCuFV(&Ux);
 setBCvFV(&Uy);
 setBCpFV( &p);

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
 double tend = 0.;
 double t;

 list<volumeField> Ux_time_series;
 list<volumeField> Uy_time_series;
 list<volumeField>  p_time_series;
  
 t = t0;

 while (t<=tend)
 {
 t = t + dt;

 for (int i = 1; i <= Nx; i++)
 {
   for(int j = 1; j<= Ny; j++)
   {

     //// convective term d(uu)/dx +  d(uv)/dy -- u-control volume	   
     fUU_f = 0.25 * (Ux.mesh[i][j].cellVal() + Ux.mesh[i][j].cellVal()) *(Ux.mesh[i+1][j].cellVal() + Ux.mesh[i+1][j].cellVal());
     fUU_b = 0.25 * (Ux.mesh[i][j].cellVal() + Ux.mesh[i][j].cellVal()) *(Ux.mesh[i-1][j].cellVal() + Ux.mesh[i-1][j].cellVal());

     fUV_f = 0.5*rho*(Ux.mesh[i ][j+1].cellVal() + Ux.mesh[i][j].cellVal() ) * 
                                  0.5* (  Uy.mesh[i][ j ].cellVal() + Uy.mesh[i+1][ j ].cellVal() );
     fUV_b = 0.5*rho*(Ux.mesh[i ][j ].cellVal() + Ux.mesh[i-1][j-1].cellVal() ) * 
                                  0.5* (  Uy.mesh[i][ j ].cellVal() + Uy.mesh[i+1][j-1].cellVal() );

     conv_u = (fUU_f - fUU_b) / (mesh[i][j].dx())  
             +(fUV_f - fUV_b) / (mesh[i][j].dy()); 

     //// convective term d(vu)/dx +  d(vv)/dy -- v-control volume     
     fVU_f = 0.5 * rho*(Ux.mesh[i][j+1].cellVal() + Ux.mesh[i] [ j].cellVal()) * 0.5 * (Uy.mesh[i+1][ j ].cellVal() + Uy.mesh[i][ j ].cellVal() ); 
     fVU_b = 0.5 * rho*(Ux.mesh[i][ j ].cellVal() + Ux.mesh[i][j-1].cellVal()) * 0.5 * (Uy.mesh[i+1][j-1].cellVal() + Uy.mesh[i][j-1].cellVal() );       


     conv_v = (fVU_f - fVU_b) / (mesh[i][j].dx())  
             +(fVV_f - fVV_b) / (mesh[i][j].dy());
    
     //// advective term dd(u)/dxx +  dd(u)/dyy -- u-control volume     
     dU_dx_f = mu * (  (Ux.mesh[i+1][j].cellVal() - Ux.mesh[ i ][j].cellVal()) / (mesh[i+1][j].xFace(1) - mesh[i+1][j].xFace(0)));
     dU_dx_b = mu * (  (Ux.mesh[ i ][j].cellVal() - Ux.mesh[i-1][j].cellVal()) / (mesh[ i ][j].xFace(1) - mesh[ i ][j].xFace(0)));
     ddU_dxx = (dU_dx_f - dU_dx_b) / (mesh[i+1][j].xCentroid() - mesh[i][j].xCentroid());


     dU_dy_f = mu * ( 0.5 * (Ux.mesh[i][j+1].cellVal() - Ux.mesh[i][ j ].cellVal()) / (mesh[i+1][j+1].yCentroid() - mesh[i+1][ j ].yCentroid()) );
     dU_dy_b = mu * ( 0.5 * (Ux.mesh[i][j  ].cellVal() - Ux.mesh[i][j-1].cellVal()) / (mesh[ i ][ j ].yCentroid() - mesh[ i ][j-1].yCentroid()) );
     ddU_dyy = (dU_dy_f - dU_dy_b) / mesh[i][j].dy();       

     diff_u = ddU_dxx + ddU_dyy;


     //// advective term dd(v)/dxx +  dd(v)/dyy -- v-control volume
     dV_dx_f = mu * (  (Uy.mesh[i+1][j].cellVal() - Uy.mesh[ i ][j].cellVal()) / (mesh[i+1][j].xCentroid() - mesh[ i ][j].xCentroid()));
     dV_dx_b = mu * (  (Uy.mesh[ i ][j].cellVal() - Uy.mesh[i-1][j].cellVal()) / (mesh[ i ][j].xCentroid() - mesh[i-1][j].xCentroid()));
     ddV_dxx = (dV_dx_f - dV_dx_b) / (mesh[i][j].dx());

     dV_dy_f = mu * (  (Uy.mesh[i][j+1].cellVal() - Uy.mesh[i][ j ].cellVal()) / (mesh[i][j+1].dy()));
     dV_dy_b = mu * (  (Uy.mesh[i][ j ].cellVal() - Uy.mesh[i][j-1].cellVal()) / (mesh[i][ j ].dy()));
     ddV_dyy = (dV_dy_f - dV_dy_b) / (mesh[i][j+1].yCentroid() - mesh[i][j].yCentroid());

     diff_v = ddV_dxx + ddV_dyy;

     //// predictor step for velocity Ux and Uy
     Ux_star.mesh[i][j].setVal( Ux.mesh[i][j].cellVal() - dt * (conv_u + diff_u) );
     Uy_star.mesh[i][j].setVal( Uy.mesh[i][j].cellVal() - dt * (conv_v + diff_v) );
   }
   cout<<"\n";
  } 

  setBCuFV(&Ux_star);
  setBCvFV(&Uy_star);

  poisson(&p, Ux_star, Uy_star, 1e-9, 0.5);

 Ux_time_series.push_back(Ux_star);
 Uy_time_series.push_back(Ux_star);

 };


 for (int j = Ny; j >=1; j--)
 {
   for(int i = 1; i<= Nx; i++)
   {
    p.mesh[i][j].print_cellVal();
   }
   cout<<"\n";
 }

 
FILE *OutFile = fopen("p.dat","w++");

 for (int j = Ny; j >=1; j--) {
   fprintf(OutFile,"\n ");
   for(int i=1; i<=Nx; i++){
      double val = p.mesh[i][j].cellVal();
      fprintf(OutFile,"\t %f ", val); } }

fclose(OutFile);




 free(mesh);
 free(Ux.mesh);
 free(Uy.mesh);
 free(p.mesh);

 



 return 0;	
}

