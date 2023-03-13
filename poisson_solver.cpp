#include <iostream>
#include <cstring>

#include "mesh.h"
#include <cmath>



///////////////////////////////////////////////////////////////////////////
void poisson (      volumeField* p, 
             const volumeField div, 
            double rho, double dt, double tol, double omega, int& niter)
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
    niter = 0;

//    cout << "start poisson iterations....\n";

    while (ratio > tol)
    {
        p_star.mesh = p->mesh;
  //      cout << "update old pressure....\n";

      //  cout << "computing l1- norm....\n";

        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
            for (int j= 1; j <= N_y; j++) {

                aP = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() )
                     - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );

                aW =  1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() );
                aE =  1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() );

                aS =  1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );
                aN =  1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() );

                if (i==1)
                {
                  aP  = aP + aW;
                  aW  = 0; 
                };

                if (i==N_x)
                {
                  aP  = aP + aE;
                  aE  = 0; 
                };

                if (j==1)
                {
                  aP  = aP + aS;
                  aS  = 0; 
                };

                if (j==N_y)
                {
                  aP  = aP + aN;
                  aN  = 0; 
                };


                p->mesh[i][j].setVal( (1.-omega)*p_star.mesh[i][j].cellVal() +
                                omega/aP * (-aE*p_star.mesh[i+1][j].cellVal() - aW*p_star.mesh[i-1][j].cellVal() 
                                            -aN*p_star.mesh[i][j+1].cellVal() - aS*p_star.mesh[i][j-1].cellVal() +
                                              (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
       // setBCpFV(p);
        l1norm =  l1_norm (*p);
        ratio = fabs( l1norm -  l1norm_old ) / fabs( l1norm_old + 1e-30);
        l1norm_old =  l1_norm (p_star);
         //cout << "ratio "<< ratio<<"\n";

    };
     setBCpFV(p);

    //cout << "convergence achieved in "<< niter<<" iterations\n";
return;
};

 volumeField residual (const volumeField& u, const volumeField& b,  double Lx, double Ly)
 //////////////////////////////////////////////////////////////////////////////////////
 // compute residual of Au = b
 // A is the coefficients' matrix which which discretizes the laplacian operator
 //////////////////////////////////////////////////////////////////////////////////////
 {
   int N_x = b.dimx;
   int N_y = b.dimy;

   volumeField  res       = buildFVfield(N_x, N_y, Lx, Ly);
   initFieldVal (res , 0.0);   

   // east and west coefficients for discretized Poisson equation
   double aE, aW, aP;
   double aN, aS;   

   for (int i=1; i<=N_x; i++){
    for (int j=1; j<=N_y; j++){

            aP = - 1./ ( (u.mesh[i+1][j].xCentroid() - u.mesh[ i ][j].xCentroid()) * u.mesh[i][j].dx() )
                 - 1./ ( (u.mesh[ i ][j].xCentroid() - u.mesh[i-1][j].xCentroid()) * u.mesh[i][j].dx() )
                 - 1./ ( (u.mesh[i][j+1].yCentroid() - u.mesh[i][ j ].yCentroid()) * u.mesh[i][j].dy() )
                 - 1./ ( (u.mesh[i][ j ].yCentroid() - u.mesh[i][j-1].yCentroid()) * u.mesh[i][j].dy() );

            aW = - 1./ ( (u.mesh[ i ][j].xCentroid() - u.mesh[i-1][j].xCentroid()) * u.mesh[i][j].dx() );
            aE = - 1./ ( (u.mesh[i+1][j].xCentroid() - u.mesh[ i ][j].xCentroid()) * u.mesh[i][j].dx() );

            aS = - 1./ ( (u.mesh[i][ j ].yCentroid() - u.mesh[i][j-1].yCentroid()) * u.mesh[i][j].dy() );
            aN = - 1./ ( (u.mesh[i][j+1].yCentroid() - u.mesh[i][ j ].yCentroid()) * u.mesh[i][j].dy() );
                
            res.mesh[i][j].setVal( -aE*u.mesh[i+1][j].cellVal() 
                                   -aW*u.mesh[i-1][j].cellVal()  
                                   +aP*u.mesh[i ][j ].cellVal() 
                                   -aN*u.mesh[i][j+1].cellVal()  
                                   -aS*u.mesh[i][j-1].cellVal()  
                                      -b.mesh[i][ j ].cellVal() );
    }
   }
  return res;
  }


 void ssor (volumeField* p, const volumeField& b, double omega, double Lx, double Ly, double tol, const char BC [1])
  //////////////////////////////////////////////////////////////////////////////////////
 // compute the solution of Au = b within niter_max iterations
 // A is the coefficients' matrix which which discretizes the laplacian operator
 //////////////////////////////////////////////////////////////////////////////////////
 {
    int N_x, N_y;
    N_x = p->dimx;
    N_y = p->dimy;

    double ratio = 1;
    double l1norm, l1norm_old;

    // east and west coefficients for discretized Poisson equation
    double aE, aW, aP;
    double aN, aS;

    // initialize mesh pointer to input pointer to pressure
    volumeField p_star;
    p_star = *p;
 
    setBCpFV(&p_star);
    int niter = 0;
    
    //// execute pre-smoother  iterations ///    
    while (ratio>tol && niter <10000)
    {
        p_star.mesh = p->mesh;
        niter+= 1;

        for (int i = 1; i <= N_x; i++) {
            for (int j= 1; j <= N_y; j++) {

                aP = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() )
                     - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() )
                     - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );

                aW = - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() );
                aE = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() );

                aS = - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );
                aN = - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() );

                p->mesh[i][j].setVal( (1.-omega)*p_star.mesh[i][j].cellVal() +
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p->mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p->mesh[i][j-1].cellVal() + b.mesh[i][j].cellVal()) );
            };
        };

        if ( strcmp(BC,"N") ) {
            setBCpFV(p);
        } 
        if(strcmp(BC,"H") ) {     setBCerrFV(p);}
        else if (!strcmp(BC,"N") && !strcmp(BC,"H")){  std::cout<<"error on SSOR BC type";
        return;}
   
        l1norm =  l1_norm (*p);
        ratio = fabs(( l1norm -  l1norm_old )) / fabs( l1norm_old );
        l1norm_old =  l1_norm (p_star); 
             // cout << "ratio "<< ratio<<"\n";
  
    };
    return;
    
 }
////////////////////////////////////////////////////////////////////////////////////////////////////

void vCycle (volumeField* u, volumeField* b, double omega, double Lx, double Ly, double tol, const int num_levels, int level)
{
                  std:: cout << "here "<< "\n";

 if (level == num_levels)
 {
  // void ssor (volumeField* p, const volumeField& b, double omega, double Lx, double Ly, double tol, const char BC [1])
    ssor( u,*b, omega, Lx, Ly, tol, "H");
                std::cout << "here "<< "\n";

    return;
 }

 // Step 1: smoothing of solution
 char BC[2];
 if      (level ==1){  strcpy(BC,"N");}
 else if (level >1 ){  strcpy(BC,"H");}

 ssor( u,*b, omega, Lx, Ly, 1e-3, BC);
 //cout<< "step 1 done ..";

 // Step 2: restrict residual to coarse grid
 volumeField  res         = buildFVfield(b->dimx  , b->dimy  , Lx, Ly);
 volumeField  res_coarser = buildFVfield(b->dimx/2, b->dimy/2, Lx, Ly);

 res              = residual (*u, *b,   Lx,  Ly);
 interpolateFieldVal (res , &res_coarser);  
  //cout<< "step 2 done ..";

 // Step 3:Solve A delta_coarser = res_coarser the coarse grid. (Recursively)
 volumeField delta_coarser = res_coarser;
  //initFieldVal(delta_coarser, 0);
 vCycle(&delta_coarser, &res_coarser,  omega,   Lx,   Ly,   tol,   num_levels, ++level);

  //cout<< "step 3 done ..";

 //Step 4: Prolongate delta_coarser to fine grid and add to u
 volumeField  delta_finer         = buildFVfield(delta_coarser.dimx*2  , delta_coarser.dimy*2 , Lx, Ly);
 prolongateFieldVal (delta_coarser, &delta_finer);  
 *u = *u + delta_finer;
  //cout<< "step 4 done ..";

 // Step 5: Relax A u=b on this grid
 ssor( u,*b, omega, Lx, Ly, 1e-6, "H");

 return;   
}
/*
 void FMG(volumeField* u, volumeField* b, double omega, double Lx, double Ly, double tol, const int num_levels, int level=1)
 {

 if (level == num_levels)
 {
    ssor( u,*b, omega, Lx, Ly, tol, "N");
    return;
 }

  // Step 1: Restrict the rhs to a coarse grid
  volumeField  b_coarser         = buildFVfield(b->dimx/2  , b->dimy/2 , Lx, Ly);
  interpolateFieldVal (*b , &b_coarser);  

  // Step 2: Solve the coarse grid problem using FMG
  volumeField  u_coarser         = buildFVfield(b->dimx/2  , b->dimy/2 , Lx, Ly);
  FMG(&u_coarser, &b_coarser,   omega,   Lx,   Ly,   tol,  num_levels, level++);

  // Step 3: Interpolate u_c to the fine grid
  u_finer = buildFVfield(u_coarser.dimx*2  , u_coarser.dimy*2 , Lx, Ly);
  prolongateFieldVal (u_coarser, &u_finer);  

  //step 4: Execute '6' V-cycles
  for (int i = 1; i<=6; i++)
  { 
   vCycle(&u_finer, b,  omega,   Lx,   Ly,   tol,   num_levels-level, level++);
  }
  *u =u_finer; //,res
  return;
 }

 */


 ///////////////////////////////////////////////////////////////////////////
 /*
void poisson_direct_solver(volumeField* p, const volumeField b)
{
 int N_x = p->dimx;
 int N_y = p->dimy;

 int n = N_x*N_y;

 double A[n][n];
 double B[n];

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {

    A[i*n][j*n]

return;
};
 */
