#include <iostream>
#include "mesh.h"
using namespace std;

void predictor_u (volumeField* Ux_pred, 
                  const volumeField& Ux, 
                  const volumeField& Uy, 
                  const volumeField& Uz,
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
    N_z = Ux_pred->dimz;

    // east and west coefficients for discretized Poisson equation

    niter = 0;

    // initialize mesh pointer to input pointer to predicted velocity 
    volumeField Ux_pred_star;
    Ux_pred_star = *Ux_pred;
 
    setBCuFV(&Ux_pred_star);

    double aW, aN, aS, aE, aF, aB, aP;
 
    double xF_i, xC_i, xC_i1, dx;
    double yF_j, yC_j, yC_j1, dy;
    double zF_j, zC_j, zC_j1, dz;

    double lx, ly_f, ly_b, lz_f, lz_b;   
    double Fn;

    double mu_f, mu_b;

    while (ratio > tol)
    {
        Ux_pred_star.mesh = Ux_pred->mesh;

        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
          for (int j= 1; j <= N_y; j++) {
            for (int k= 1; k <= N_z; k++) {

               mu_f = mu;
               mu_b = mu;

               xF_i  = Ux.mesh[i][j][k].xFace(1);
               xC_i  = Ux.mesh[i][j][k].xCentroid();
               xC_i1 = Ux.mesh[i+1][j][k].xCentroid();
               lx    = (xF_i - xC_i)/(xC_i1 - xC_i);

               yF_j  = Uy.mesh[i][j][k].yFace(1);
               yC_j  = Uy.mesh[i][j][k].yCentroid();
               yC_j1 = Uy.mesh[i][j+1][k].yCentroid();
               ly_f =  (yF_j - yC_j)/(yC_j1 - yC_j);

               yF_j  = Uy.mesh[i][j][k].yFace(0);
               yC_j  = Uy.mesh[i][j-1][k].yCentroid();
               yC_j1 = Uy.mesh[i][j][k].yCentroid();
               ly_b =  (yF_j - yC_j)/(yC_j1 - yC_j);

               zF_k  = Uz.mesh[i][j][k].zFace(1);
               zC_k  = Uz.mesh[i][j][k].zCentroid();
               zC_k1 = Uz.mesh[i][j][k+1].zCentroid();
               lz_f =  (zF_j - zC_j)/(zC_j1 - zC_j);

               zF_k  = Uz.mesh[i][j][ k ].zFace(0);
               zC_k  = Uz.mesh[i][j][k-1].zCentroid();
               zC_k1 = Uz.mesh[i][j][ k ].zCentroid();
               lz_b =  (zF_k - zC_k)/(zC_k1 - zC_k);               

               aE = - 0.25 *rho * (Ux.mesh[ i ][j][k].cellVal()   + Ux.mesh[i-1][j][k].cellVal  ()) /
                                  (Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[ i ][j][k].xCentroid());

               aE = aE -  
                     mu_b /(Ux.mesh[i+1][j].xCentroid() - Ux.mesh[ i ][j].xCentroid())*Ux.mesh[i][j].dx();

               aW =   0.25 *rho * (Ux.mesh[i+1][j].cellVal() + Ux.mesh[ i ][j].cellVal()) /
                                  (Ux.mesh[i+1][j].xCentroid() - Ux.mesh[ i ][j].xCentroid());

               aW = aW -  
                     mu_f /(Ux.mesh[i+1][j].xCentroid() - Ux.mesh[i][j].xCentroid())*Ux.mesh[i+1][j].dx();


               aS = - rho* (1-ly_b) * (lx*Uy.mesh[i-1][j+1][k].cellVal() + (1.-lx)*Uy.mesh[i-1][j][k].cellVal()) / Uy.mesh[i][j].dy();

               aS = aS - 
                     mu_b /(Ux.mesh[i][j].yCentroid() - Ux.mesh[i][j-1].yCentroid())*Ux.mesh[i][j].dy();

               aN =  rho * ly_f * (lx*Uy.mesh[i][j+1][k].cellVal() + (1.-lx)*Uy.mesh[i][j][k].cellVal()) / Uy.mesh[i][j][k].dy();

               aN = aN - 
                     mu_f /(Ux.mesh[i][j+1][k].yCentroid() - Ux.mesh[i][j][k].yCentroid())*Ux.mesh[i][j][k].dy();


               aF =   rho* lz_f * (lx*Uz.mesh[i+1][j][k].cellVal() + (1.-lx)*Uy.mesh[i][j][k].cellVal()) / Uz.mesh[i][j][k].dz();

               aB = - rho* (1.-lz_b) * (lx*Uz.mesh[i+1][j][k-1].cellVal() + (1.-lx)*Uz.mesh[i][j][k-1].cellVal()) / Uz.mesh[i][j][k].dz();

               aP =  0.25 *rho * (Ux.mesh[i+1][j][k].cellVal()   + Ux.mesh[ i ][j][k].cellVal())   /
                                 (Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[ i ][j][k].xCentroid())  
                   - 0.25 *rho * (Ux.mesh[ i ][j][k].cellVal()   + Ux.mesh[i-1][j][k].cellVal  ()) /
                                 (Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[ i ][j][k].xCentroid()) 
                 + rho *   ly_f   * (lx*Uy.mesh[ i ][j+1][k].cellVal() + (1.-lx)*Uy.mesh[ i ][j][k].cellVal()) / Uy.mesh[i][j][k].dy() 
                 - rho* (1.-ly_b) * (lx*Uy.mesh[i-1][j+1][k].cellVal() + (1.-lx)*Uy.mesh[i-1][j][k].cellVal()) / Uy.mesh[i][j][k].dy()   
                 + rho *  lz_f    * (lx*Uz.mesh[i+1][j] [k ].cellVal() + (1.-lx)*Uy.mesh[i][j][ k ].cellVal()) / Uz.mesh[i][j][k].dz() 
                 - rho* (1.-lz_b) * (lx*Uz.mesh[i+1][j][k-1].cellVal() + (1.-lx)*Uz.mesh[i][j][k-1].cellVal()) / Uz.mesh[i][j][k].dz();

                aP = aP + 
                     mu_b /(Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[i][ j ][k].xCentroid())*Ux.mesh[ i ][j][k].dx() +
                     mu_f /(Ux.mesh[i+1][j][k].xCentroid() - Ux.mesh[i][ j ][k].xCentroid())*Ux.mesh[i+1][j][k].dx() +
                     mu_b /(Ux.mesh[i][ j ][k].yCentroid() - Ux.mesh[i][j-1][k].yCentroid())*Ux.mesh[ i ][j][k].dy() +
                     mu_f /(Ux.mesh[i][j+1][k].yCentroid() - Ux.mesh[i][ j ][k].yCentroid())*Ux.mesh[ i ][j][k].dy();

               
                Fn = rho*Ux.mesh[i][j][k].cellVal() - 0.5 * dt* conv_expl.mesh[i][j][k].cellVal() + 0.5 * dt*diff_expl.mesh[i][j][k].cellVal();

                Ux_pred->mesh[i][j][k].setVal( 
                                Fn / (rho+0.5*dt*aP) +
                              - (0.5*dt/(rho+0.5*dt*aP)) *
                                        ( aE*Ux_pred_star.mesh[i-1][ j ][k].cellVal() +
                                          aW*Ux_pred_star.mesh[i+1][ j ][k].cellVal() + 
                                          aN*Ux_pred_star.mesh[ i ][j+1][k].cellVal() + 
                                          aS*Ux_pred_star.mesh[ i ][j-1][k].cellVal() + 
                                          aF*Ux_pred_star.mesh[ i ][j][k+1].cellVal() + 
                                          aB*Ux_pred_star.mesh[ i ][j][k-1].cellVal() ));                                      


            };
        };
        setBCuFV(Ux_pred);
        l1norm =  l1_norm (*Ux_pred);
        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old );
        l1norm_old =  l1norm; //l1_norm (Ux_pred);

    };
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
void  predictor_v (volumeField* Uy_pred, 
                  const volumeField& Ux, 
                  const volumeField& Uy, 
                  const volumeField& Uz, 
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
    N_z = Uy_pred->dimz;

    // east and west coefficients for discretized Poisson equation
    double aW, aN, aS, aE, aP, aF, aB;


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
            for (int k= 1; k <= N_z; k++) {

               xF_i  = Uy.mesh[ i ][j][k].xFace(1);
               xC_i  = Uy.mesh[ i ][j][k].xCentroid();
               xC_i1 = Uy.mesh[i+1][j][k].xCentroid();
               lx_f  = (xF_i - xC_i)/(xC_i1 - xC_i);

               xF_i  = Uy.mesh[ i ][j][k].xFace(0);
               xC_i  = Uy.mesh[i-1][j][k].xCentroid();
               xC_i1 = Uy.mesh[ i ][j][k].xCentroid();
               lx_b = (xF_i - xC_i)/(xC_i1 - xC_i);   

               zF_k  = Uz.mesh[i][j][k].zFace(1);
               zC_k  = Uz.mesh[i][j][k].zCentroid();
               zC_k1 = Uz.mesh[i][j][k+1].zCentroid();
               lz_f =  (zF_j - zC_j)/(zC_j1 - zC_j);

               zF_k  = Uz.mesh[i][j][ k ].zFace(0);
               zC_k  = Uz.mesh[i][j][k-1].zCentroid();
               zC_k1 = Uz.mesh[i][j][ k ].zCentroid();
               lz_b =  (zF_k - zC_k)/(zC_k1 - zC_k);                              

               yF_j  = Uy.mesh[i][ j ][k].yFace(1);
               yC_j  = Uy.mesh[i][ j ][k].yCentroid();
               yC_j1 = Uy.mesh[i][j+1][k].yCentroid();
               ly =  (yF_j - yC_j)/(yC_j1 - yC_j);
              
               aE =   rho *    lx_f   * (ly*Ux.mesh[ i ][j+1][k].cellVal() + (1-ly)*Ux.mesh [i ][j][k].cellVal() ) / Uy.mesh[i][j][k].dx();

               aW = - rho * (1.-lx_b) * (ly*Ux.mesh[i-1][j+1][k].cellVal() + (1-ly)*Ux.mesh[i-1][j][k].cellVal() ) / Uy.mesh[i][j][k].dx();
 
               aN =  0.25* rho * (Uy.mesh[i][j+1][k].cellVal() + Uy.mesh[i][ j ][k].cellVal()) / (Uy.mesh[i][j+1][k].yCentroid() - Uy.mesh[i][j][k].yCentroid());

               aS = -0.25* rho * (Uy.mesh[i][ j ][k].cellVal() + Uy.mesh[i][j-1][k].cellVal()) / (Uy.mesh[i][j+1][k].yCentroid() - Uy.mesh[i][j][k].yCentroid());

               aF =   rho *    lz_f   * (ly*Uz.mesh[ i ][j+1][k].cellVal() + (1-ly)*Uz.mesh [i ][j][k].cellVal() ) / Uy.mesh[i][j][k].dz();

               aB = - rho * (1.-lz_b) * (ly*Uz.mesh[i][j+1][k-1].cellVal() + (1-ly)*Ux.mesh[i][j][k-1].cellVal() ) / Uy.mesh[i][j][k].dz();               

               aP =   rho *    lx_f   * (ly*Ux.mesh[ i ][j+1][k].cellVal() + (1.-ly)*Ux.mesh [i ][j][k].cellVal() ) / Uy.mesh[i][j][k].dx() 
                    - rho * (1.-lx_b) * (ly*Ux.mesh[i-1][j+1][k].cellVal() + (1.-ly)*Ux.mesh[i-1][j][k].cellVal() ) / Uy.mesh[i][j][k].dx()
                    + 0.25* rho * (Uy.mesh[i][j+1][k].cellVal() + Uy.mesh[i][ j ][k].cellVal()) / (Uy.mesh[i][j+1][k].yCentroid() - Uy.mesh[i][j][k].yCentroid())
                    - 0.25* rho * (Uy.mesh[i][ j ][k].cellVal() + Uy.mesh[i][j-1][k].cellVal()) / (Uy.mesh[i][j+1][k].yCentroid() - Uy.mesh[i][j][k].yCentroid())
                    + rho *    lz_f   * (ly*Uz.mesh[ i ][j+1][k].cellVal() + (1.-ly)*Uz.mesh [i ][j][k].cellVal() ) / Uy.mesh[i][j][k].dz()
                    - rho * (1.-lz_b) * (ly*Uz.mesh[i][j+1][k-1].cellVal() + (1.-ly)*Ux.mesh[i][j][k-1].cellVal() ) / Uy.mesh[i][j][k].dz();

                Fn = rho * Uy.mesh[i][j][k].cellVal() - 0.5 * dt*conv_v_expl.mesh[i][j][k].cellVal() + dt*diff_v_expl.mesh[i][j][k].cellVal();

                Uy_pred->mesh[i][j][k].setVal( 
                                     Fn/(rho + 0.5*dt*aP) +
                              - (0.5*dt/(rho + 0.5*dt*aP)) *
                                        ( aE*Uy_pred_star.mesh[i-1][ j ][k].cellVal() +
                                          aW*Uy_pred_star.mesh[i+1][ j ][k].cellVal() + 
                                          aN*Uy_pred_star.mesh[ i ][j+1][k].cellVal() + 
                                          aS*Uy_pred_star.mesh[ i ][j-1][k].cellVal() +      
                                          aF*Uy_pred_star.mesh[ i ][j][k+1].cellVal() + 
                                          aB*Uy_pred_star.mesh[ i ][j][k-1].cellVal()  ) ) ;
                                        

                };
            };
        };
        setBCvFV(Uy_pred);
        l1norm =  l1_norm (*Uy_pred);
        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old );
        l1norm_old =  l1norm; 
    };
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void  predictor_w (volumeField* Uz_pred, 
                  const volumeField& Ux, 
                  const volumeField& Uy, 
                  const volumeField& Uz, 
                  const volumeField& conv_w_expl, 
                  const volumeField& diff_w_expl, 
                  const double rho, const  double dt, const  double mu, const double tol, int &niter )
{

    /// Uz      velocity at old time step
    /// Uz_pred velocity at new time step

    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y;
    N_x = Uz_pred->dimx;
    N_y = Uz_pred->dimy;
    N_z = Uz_pred->dimz;

    // east and west coefficients for discretized Poisson equation
    double aW, aN, aS, aE, aP, aF, aB;


    // initialize mesh pointer to input pointer to predicted velocity 
    volumeField Uz_pred_star;
    Uz_pred_star = *Uz_pred;
 
    setBCwFV(&Uz_pred_star);
    niter = 0;

    double xF_i, xC_i, xC_i1;
    double yF_j, yC_j, yC_j1;
    double lx_b, lx_f, ly_b, ly_f, lz;

    double Fn;

    while (ratio > tol)
    {
        Uy_pred_star.mesh = Uy_pred->mesh;

        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
          for (int j= 1; j <= N_y; j++) {
            for (int k= 1; k <= N_z; k++) {

               xF_i  = Uz.mesh[ i ][j][k].xFace(1);
               xC_i  = Uz.mesh[ i ][j][k].xCentroid();
               xC_i1 = Uz.mesh[i+1][j][k].xCentroid();
               lx_f  = (xF_i - xC_i)/(xC_i1 - xC_i);

               xF_i  = Uz.mesh[ i ][j][k].xFace(0);
               xC_i  = Uz.mesh[i-1][j][k].xCentroid();
               xC_i1 = Uz.mesh[ i ][j][k].xCentroid();
               lx_b = (xF_i - xC_i)/(xC_i1 - xC_i);   

               yF_j  = Uy.mesh[i][ j ][k].yFace(1);
               yC_j  = Uy.mesh[i][ j ][k].yCentroid();
               yC_j1 = Uy.mesh[i][j+1][k].yCentroid();
               ly_f =  (yF_j - yC_j)/(yC_j1 - yC_j);
                              
               yF_j  = Uy.mesh[i][ j ][k].yFace(0);
               yC_j  = Uy.mesh[i][ j ][k].yCentroid();
               yC_j1 = Uy.mesh[i][j+1][k].yCentroid();
               ly_b =  (yF_j - yC_j)/(yC_j1 - yC_j);

               zF_k  = Uz.mesh[i][j][k].zFace(1);
               zC_k  = Uz.mesh[i][j][k].zCentroid();
               zC_k1 = Uz.mesh[i][j][k+1].zCentroid();
               lz   =  (zF_j - zC_j)/(zC_j1 - zC_j);                       

             
               aE =   rho *    lx_f   * (lz*Ux.mesh[ i ][j][k+1].cellVal() + (1-lz)*Ux.mesh [ i ][j][k].cellVal() ) / Uz.mesh[i][j][k].dx();
               aW = - rho *  (1-lx_b) * (lz*Ux.mesh[i-1][j][k+1].cellVal() + (1-lz)*Ux.mesh [i-1][j][k].cellVal() ) / Uz.mesh[i][j][k].dx();
 
               aN =   rho *    ly_f   * (lz*Uy.mesh[ i ][j][k+1].cellVal() + (1-lz)*Uy.mesh [i][ j ][k].cellVal() ) / Uz.mesh[i][j][k].dy();
               aS = - rho *  (1-ly_b) * (lz*Uy.mesh[i][j-1][k+1].cellVal() + (1-lz)*Uy.mesh [i][j-1][k].cellVal() ) / Uz.mesh[i][j][k].dy();

               aF =   rho * 0.25 * (Uz.mesh[i][j][k+1].cellVal() + Uz.mesh [i][j][k].cellVal() ) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid());
               aB = - rho * 0.25 * (Uz.mesh[i][j][k].cellVal() + Uz.mesh [i][j][k-1].cellVal() ) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid());

               aP =   rho *    lx_f   * (lz*Ux.mesh[ i ][j][k+1].cellVal() + (1-lz)*Ux.mesh [ i ][j][k].cellVal() ) / Uz.mesh[i][j][k].dx()
                    - rho *  (1-lx_b) * (lz*Ux.mesh[i-1][j][k+1].cellVal() + (1-lz)*Ux.mesh [i-1][j][k].cellVal() ) / Uz.mesh[i][j][k].dx()
                    + rho *    ly_f   * (lz*Uy.mesh[ i ][j][k+1].cellVal() + (1-lz)*Uy.mesh [i][ j ][k].cellVal() ) / Uz.mesh[i][j][k].dy()
                    - rho *  (1-ly_b) * (lz*Uy.mesh[i][j-1][k+1].cellVal() + (1-lz)*Uy.mesh [i][j-1][k].cellVal() ) / Uz.mesh[i][j][k].dy()
                    + rho * 0.25 * (Uz.mesh[i][j][k+1].cellVal() + Uz.mesh [i][j][k].cellVal() ) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid())
                    - rho * 0.25 * (Uz.mesh[i][j][k].cellVal() + Uz.mesh [i][j][k-1].cellVal() ) / (Uz.mesh[i][j][k+1].zCentroid() - Uz.mesh[i][j][k].zCentroid());

                Fn = rho * Uz.mesh[i][j][k].cellVal() - 0.5 * dt*conv_w_expl.mesh[i][j][k].cellVal() + dt*diff_w_expl.mesh[i][j][k].cellVal();

                Uz_pred->mesh[i][j][k].setVal( 
                                     Fn/(rho + 0.5*dt*aP) +
                              - (0.5*dt/(rho + 0.5*dt*aP)) *
                                        ( aE*Uz_pred_star.mesh[i-1][ j ][k].cellVal() +
                                          aW*Uz_pred_star.mesh[i+1][ j ][k].cellVal() + 
                                          aN*Uz_pred_star.mesh[ i ][j+1][k].cellVal() + 
                                          aS*Uz_pred_star.mesh[ i ][j-1][k].cellVal() +      
                                          aF*Uz_pred_star.mesh[ i ][j][k+1].cellVal() + 
                                          aB*Uz_pred_star.mesh[ i ][j][k-1].cellVal()  ) ) ;
                                        

                };
            };
        };
        setBCwFV(Uz_pred);
        l1norm =  l1_norm (*Uz_pred);
        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old );
        l1norm_old =  l1norm; 
    };
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////