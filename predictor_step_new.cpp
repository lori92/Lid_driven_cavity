#include <iostream>
#include <cmath>
#include "mesh.h"

void predictor_u(volumeField* Ux_pred, 
                  const volumeField& Ux, 
                  const volumeField& Uy, 
                   const volumeField& conv_expl, 
                  const volumeField& diff_expl, 
                  const double rho, const  double dt, 
                  const double mu,  const double tol, 
                  int& niter )
{

    /// Ux      velocity at old time step
    /// Ux_pred velocity at new time step

    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y, N_z;
    N_x = Ux_pred->dimx;
    N_y = Ux_pred->dimy;
 
    // east and west coefficients for discretized Poisson equation

    niter = 0;

    // initialize mesh pointer to input pointer to predicted velocity 
    volumeField Ux_pred_star;
    Ux_pred_star = *Ux_pred;
 
    setBCuFV(&Ux_pred_star);

    double aW, aE, aS, aN, aP;
 
    double xF_i, xC_i, xC_i1, dx;
    double yF_j, yC_j, yC_j1, dy;

    double lx, ly_f, ly_b;   
    double Fn;

    double mu_f, mu_b;

    while (ratio > tol)
    {
        Ux_pred_star.mesh = Ux_pred->mesh;

        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
          for (int j= 1; j <= N_y; j++) {

               mu_f = mu;
               mu_b = mu;

               xF_i  = Ux.mesh[i][j] .xFace(1);
               xC_i  = Ux.mesh[i][j] .xCentroid();
               xC_i1 = Ux.mesh[i+1][j] .xCentroid();
               lx    = (xF_i - xC_i)/(xC_i1 - xC_i);
 

               yF_j  = Uy.mesh[i][j] .yFace(1);
               yC_j  = Uy.mesh[i][j] .yCentroid();
               yC_j1 = Uy.mesh[i][j+1] .yCentroid();
               ly_f =  (yF_j - yC_j)/(yC_j1 - yC_j);

 

               yF_j  = Uy.mesh[i][j] .yFace(0);
               yC_j  = Uy.mesh[i][j-1] .yCentroid();
               yC_j1 = Uy.mesh[i][j] .yCentroid();
               ly_b =  (yF_j - yC_j)/(yC_j1 - yC_j);

 

               aE =   0.25 *rho * (Ux.mesh[i+1][j] .cellVal()   + 
                                   Ux.mesh[ i ][j] .cellVal()) /
                                  (Ux.mesh[i+1][j] .xCentroid() -
                                   Ux.mesh[ i ][j] .xCentroid());

               aE = aE -  
                     mu_f /((Ux.mesh[i+1][j] .xCentroid() - 
                            Ux.mesh[ i ][j] .xCentroid()) * 
                            Ux.mesh[i+1][j] .dx());

               aW = - 0.25 *rho * (Ux.mesh[i-1][j] .cellVal() + 
                                   Ux.mesh[ i ][j] .cellVal()) /
                                  (Ux.mesh[i+1][j] .xCentroid() - 
                                   Ux.mesh[ i ][j] .xCentroid());

               aW = aW -  
                     mu_b /((Ux.mesh[i+1][j] .xCentroid() - 
                            Ux.mesh[ i ][j] .xCentroid())*
                            Ux.mesh[i][j] .dx());


               aS = -  rho* (1.-ly_b) * (lx*   Uy.mesh[i+1][j-1] .cellVal() + 
                                       (1.-lx)*Uy.mesh[i][ j-1 ] .cellVal()) / 
                                       Uy.mesh[i][j] .dy();

               aS = aS - 
                     mu_b /((Ux.mesh[i][j] .yCentroid() -
                            Ux.mesh[i][j-1] .yCentroid())*Ux.mesh[i][j] .dy());

               aN =   rho * ly_f * (lx*Uy.mesh[i+1][j] .cellVal() +
                               (1.-lx)*Uy.mesh[i][j] .cellVal()) / 
                               Uy.mesh[i][j] .dy();

               aN = aN - 
                     mu_f /((Ux.mesh[i][j+1] .yCentroid() - 
                            Ux.mesh[i][j] .yCentroid())*Ux.mesh[i][j] .dy());
 

       
               aP =  0.25 *rho * (Ux.mesh[i+1][j] .cellVal()   + Ux.mesh[ i ][j] .cellVal())   /
                                 (Ux.mesh[i+1][j] .xCentroid() - Ux.mesh[ i ][j] .xCentroid())  
                   - 0.25 *rho * (Ux.mesh[ i ][j] .cellVal()   + Ux.mesh[i-1][j] .cellVal  ()) /
                                 (Ux.mesh[i+1][j] .xCentroid() - Ux.mesh[ i ][j] .xCentroid()) 
                 + rho *(1.-ly_f) * (lx*Uy.mesh[i+1][ j ] .cellVal() + (1.-lx)*Uy.mesh[i][ j ] .cellVal()) / Uy.mesh[i][j] .dy() 
                 - rho* (   ly_b) * (lx*Uy.mesh[i+1][j-1] .cellVal() + (1.-lx)*Uy.mesh[i][j-1] .cellVal()) / Uy.mesh[i][j] .dy() ; 
  
                aP = aP +  (
                     mu_b /((Ux.mesh[i+1][j] .xCentroid() - Ux.mesh[i][ j ] .xCentroid())*Ux.mesh[ i ][j] .dx()) +
                     mu_f /((Ux.mesh[i+1][j] .xCentroid() - Ux.mesh[i][ j ] .xCentroid())*Ux.mesh[i+1][j] .dx()) +
                     mu_b /((Ux.mesh[i][ j ] .yCentroid() - Ux.mesh[i][j-1] .yCentroid())*Ux.mesh[ i ][j] .dy()) +
                     mu_f /((Ux.mesh[i][j+1] .yCentroid() - Ux.mesh[i][ j ] .yCentroid())*Ux.mesh[ i ][j] .dy()) );

 
                Fn = rho*Ux.mesh[i][j] .cellVal() - 0.5 * dt * conv_expl.mesh[i][j] .cellVal() 
                                                  + 0.5 * dt * diff_expl.mesh[i][j] .cellVal();

/*
                if (i==1)
                {
                  aP  = aP; // + aW;
                  aW  = 0; 
                };

                if (i==N_x)
                {
                  aP  = aP; // +aE;
                  aE  = 0; 
                };

                if (j==1)
                {
                  aP  = aP - aS;
                  aS  = 0; 
                };

                if (j==N_y)
                {
                  aP  = aP - aN;
                  Fn = Fn + 0.2 * aN * dt *0.5;
                  aN  = 0; 
                };
*/
                Ux_pred->mesh[i][j] .setVal( 
                                Fn / (rho+0.5*dt*aP) 
                              - (0.5*dt/(rho+0.5*dt*aP)) *
                                        ( aE*Ux_pred_star.mesh[i+1][ j ] .cellVal() +
                                          aW*Ux_pred_star.mesh[i-1][ j ] .cellVal() + 
                                          aN*Ux_pred_star.mesh[ i ][j+1] .cellVal() + 
                                          aS*Ux_pred_star.mesh[ i ][j-1] .cellVal() ));       

 
                               
            };
        };
        setBCuFV(Ux_pred);
 
        l1norm =  l1_norm (*Ux_pred);
 
        ratio = fabs( l1norm -  l1norm_old ) / fabs( l1norm_old +1e-30);
 
        l1norm_old =  l1norm; 
    };
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
void  predictor_v(volumeField* Uy_pred, 
                  const volumeField& Ux, 
                  const volumeField& Uy, 
                  const volumeField& conv_v_expl, 
                  const volumeField& diff_v_expl, 
                  const double rho, const  double dt, const  double mu, const double tol, int &niter )
{

    /// Uy      velocity at old time step
    /// Uy_pred velocity at new time step

    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y;
    N_x = Uy_pred->dimx;
    N_y = Uy_pred->dimy;

    // east and west coefficients for discretized Poisson equation
    double aW, aE, aS, aN, aP;

    double mu_f, mu_b;

    // initialize mesh pointer to input pointer to predicted velocity 
    volumeField Uy_pred_star;
    Uy_pred_star = *Uy_pred;
 
    setBCvFV(&Uy_pred_star);
    niter = 0;

    double xF_i, xC_i, xC_i1;
    double yF_j, yC_j, yC_j1;
    double lx_b, lx_f, ly;

    double Fn;

    while (ratio > tol)
    {
        Uy_pred_star.mesh = Uy_pred->mesh;

        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
          for (int j= 1; j <= N_y; j++) {

               xF_i  = Uy.mesh[ i ][j] .xFace(1);
               xC_i  = Uy.mesh[ i ][j] .xCentroid();
               xC_i1 = Uy.mesh[i+1][j] .xCentroid();
               lx_f  = (xF_i - xC_i)/(xC_i1 - xC_i);

               xF_i  = Uy.mesh[ i ][j] .xFace(0);
               xC_i  = Uy.mesh[i-1][j] .xCentroid();
               xC_i1 = Uy.mesh[ i ][j] .xCentroid();
               lx_b = (xF_i - xC_i)/(xC_i1 - xC_i);                   

               yF_j  = Uy.mesh[i][ j ] .yFace(1);
               yC_j  = Uy.mesh[i][ j ] .yCentroid();
               yC_j1 = Uy.mesh[i][j+1] .yCentroid();
               ly =  (yF_j - yC_j)/(yC_j1 - yC_j);

              mu_f = mu;
              mu_b = mu;
              
               aE =   0.5*rho *    lx_f   * (ly*Ux.mesh[ i ][j+1] .cellVal() + (1-ly)*Ux.mesh [i ][j] .cellVal() ) / Uy.mesh[i][j] .dx();

               aE = aE - mu_f /((Ux.mesh[i+1][j] .xCentroid() - 
                                 Ux.mesh[ i ][j] .xCentroid()) * 
                                  Ux.mesh[i][j] .dx());

               aW = - 0.5*rho * (1.-lx_b) * (ly*Ux.mesh[i-1][j+1] .cellVal() + (1-ly)*Ux.mesh[i-1][j] .cellVal() ) / Uy.mesh[i][j] .dx();
 
               aW = aW - mu_f /((Ux.mesh[i][j] .xCentroid() - 
                                 Ux.mesh[i-1][j] .xCentroid()) * 
                                  Ux.mesh[i][j] .dx());

               aN =  0.125* rho * (Uy.mesh[i][j+1] .cellVal() + Uy.mesh[i][ j ] .cellVal()) / (Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j] .yCentroid());

               aN = aN - mu_f /((Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j].yCentroid()) * Uy.mesh[i][j+1].dy());

               aS = -0.125* rho * (Uy.mesh[i][ j ] .cellVal() + Uy.mesh[i][j-1] .cellVal()) / (Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j] .yCentroid());

               aS = aS - mu_f /((Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j].yCentroid()) * Uy.mesh[i][j].dy());

               aP =   0.5*rho *    lx_f   * (ly*Ux.mesh[ i ][j+1] .cellVal() + (1.-ly)*Ux.mesh [i ][j] .cellVal() ) / Uy.mesh[i][j] .dx() 
                    - 0.5*rho * (1.-lx_b) * (ly*Ux.mesh[i-1][j+1] .cellVal() + (1.-ly)*Ux.mesh[i-1][j] .cellVal() ) / Uy.mesh[i][j] .dx()
                    + 0.125* rho * (Uy.mesh[i][j+1] .cellVal() + Uy.mesh[i][ j ] .cellVal()) / (Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j] .yCentroid())
                    - 0.125* rho * (Uy.mesh[i][ j ] .cellVal() + Uy.mesh[i][j-1] .cellVal()) / (Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[i][j] .yCentroid());

              aP  = aP + (
                    mu_f /((Ux.mesh[i+1][j] .xCentroid() - Ux.mesh[ i ][j] .xCentroid()) * Ux.mesh[i][j] .dx()) +
                    mu_f /((Ux.mesh[ i ][j] .xCentroid() - Ux.mesh[i-1][j] .xCentroid()) * Ux.mesh[i][j] .dx()) +
                    mu_f /((Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[ i ][j].yCentroid()) * Uy.mesh[i][j+1].dy()) +
                    mu_f /((Uy.mesh[i][j+1] .yCentroid() - Uy.mesh[ i ][j].yCentroid()) * Uy.mesh[i][ j ].dy()) );

                Fn = rho * Uy.mesh[i][j] .cellVal() - 0.5 * dt*conv_v_expl.mesh[i][j] .cellVal() + 0.5*dt*diff_v_expl.mesh[i][j] .cellVal();

                if (i==1)
                {
                  aP  = aP - aW;
                  aW  = 0; 
                };

                if (i==N_x)
                {
                  aP  = aP - aE;
                  aE  = 0; 
                };

                if (j==1)
                {
                  aP  = aP; // - aS;
                  aS  = 0; 
                };

                if (j==N_y)
                {
                  aP  = aP; // - aN;
                  aN  = 0; 
                };

                Uy_pred->mesh[i][j] .setVal( 
                                     Fn/(rho + 0.5*dt*aP) +
                              - (0.5*dt/(rho + 0.5*dt*aP)) *
                                        ( aE*Uy_pred_star.mesh[i+1][ j ] .cellVal() +
                                          aW*Uy_pred_star.mesh[i-1][ j ] .cellVal() + 
                                          aN*Uy_pred_star.mesh[ i ][j+1] .cellVal() + 
                                          aS*Uy_pred_star.mesh[ i ][j-1] .cellVal()  ) ) ;
                                        

                
            };
        };
        setBCvFV(Uy_pred);
        l1norm =  l1_norm (*Uy_pred);
        ratio = fabs( l1norm -  l1norm_old ) / fabs( l1norm_old +1e-30);
        l1norm_old =  l1norm; 
    };
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////