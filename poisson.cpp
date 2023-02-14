#include <iostream>
#include "mesh.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
double l1_norm(volumeField const field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;
    int N_z = field.dimz;

    double accum = 0.;

    for (int i = 1; i <= N_x; i++)
    {
     for (int j= 1; j <= N_y; j++)
     {
      for (int k= 1; k <= N_z; k++)
      {
          accum += abs( field.mesh[i][j][k].cellVal() ) ;
      };
     };
    };
    return sqrt(accum);
};
///////////////////////////////////////////////////////////////////
double sum (const volumeField& field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;
    int N_z = field.dimz;

    double accum = 0.;

    for (int i = 1; i <= N_x; i++)
    {
     for (int j= 1; j <= N_y; j++)
     {
      for (int k= 1; k <= N_z; k++)
      {
          accum += ( field.mesh[i][j][k].cellVal() ) ;
      };
     };
    };
    return (accum);
};
///////////////////////////////////////////////////////////////////////////
void poisson (      volumeField* p, 
             const volumeField div, 
            double rho, double dt, double tol, double omega, int& niter)
{
    double ratio = 1;
    double l1norm, l1norm_old;

    int N_x, N_y, N_z;
    N_x = p->dimx;
    N_y = p->dimy;
    N_z = p->dimz;

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

                aW = - 1./ ( (p_star.mesh[ i ][j].xCentroid() - p_star.mesh[i-1][j].xCentroid()) * p_star.mesh[i][j].dx() );
                aE = - 1./ ( (p_star.mesh[i+1][j].xCentroid() - p_star.mesh[ i ][j].xCentroid()) * p_star.mesh[i][j].dx() );

                aS = - 1./ ( (p_star.mesh[i][ j ].yCentroid() - p_star.mesh[i][j-1].yCentroid()) * p_star.mesh[i][j].dy() );
                aN = - 1./ ( (p_star.mesh[i][j+1].yCentroid() - p_star.mesh[i][ j ].yCentroid()) * p_star.mesh[i][j].dy() );

                p->mesh[i][j].setVal( (1.-omega)*p_star.mesh[i][j].cellVal() +
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal() + (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
        setBCpFV(p);
        l1norm =  l1_norm (*p);
        ratio = ( l1norm -  l1norm_old ) / ( l1norm_old );
        l1norm_old =  l1_norm (p_star);

    };
    //cout << "convergence achieved in "<< niter<<" iterations\n";
return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
void multigrid(volumeField* p, const volumeField& div,
            const double tol, const double rho  ,
            const double dt , const double omega,
            const double Lx , const double Ly ,
            const int ncycle)
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
    
    //// execute pre-smoother 10 iterations ///    
    while (niter < 20)
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
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal() + (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
        setBCpFV(p);
    };

    cout << "pre-smoothing performed in "<< niter <<" iterations\n";


    //// compute residual ////
    volumeField  res      = buildFVfield(N_x, N_y, Lx, Ly);
    initFieldVal(res , 0.0);

    for (int i = 1; i <= N_x; i++) {
        for (int j= 1; j <= N_y; j++) {

                aP = - 1./ ( (res.mesh[i+1][j].xCentroid() - res.mesh[ i ][j].xCentroid()) * res.mesh[i][j].dx() )
                     - 1./ ( (res.mesh[ i ][j].xCentroid() - res.mesh[i-1][j].xCentroid()) * res.mesh[i][j].dx() )
                     - 1./ ( (res.mesh[i][j+1].yCentroid() - res.mesh[i][ j ].yCentroid()) * res.mesh[i][j].dy() )
                     - 1./ ( (res.mesh[i][ j ].yCentroid() - res.mesh[i][j-1].yCentroid()) * res.mesh[i][j].dy() );

                aW = - 1./ ( (res.mesh[ i ][j].xCentroid() - res.mesh[i-1][j].xCentroid()) * res.mesh[i][j].dx() );
                aE = - 1./ ( (res.mesh[i+1][j].xCentroid() - res.mesh[ i ][j].xCentroid()) * res.mesh[i][j].dx() );

                aS = - 1./ ( (res.mesh[i][ j ].yCentroid() - res.mesh[i][j-1].yCentroid()) * res.mesh[i][j].dy() );
                aN = - 1./ ( (res.mesh[i][j+1].yCentroid() - res.mesh[i][ j ].yCentroid()) * res.mesh[i][j].dy() );

            res.mesh[i][j].setVal( -aE*p->mesh[i+1][j].cellVal() 
                                   -aW*p->mesh[i-1][j].cellVal()  
                                   +aP*p->mesh[i ][j ].cellVal() 
                                   -aN*p->mesh[i][j+1].cellVal()  
                                   -aS*p->mesh[i][j-1].cellVal()  
                                   + (rho/dt)*div.mesh[i][j].cellVal() );
        };
    };
    setBCpFV(&res);

    cout << "residual initialized ..\n";


    // interpolate res on the coarser grid //
    volumeField  res_coarser = buildFVfield(N_x/2, N_y/2, Lx, Ly);
    interpolateFieldVal(res , &res_coarser); 

    //// solve for the error delta on the coarser grid ////
    volumeField  delta       = buildFVfield(N_x/2, N_y/2, Lx, Ly);
    initFieldVal(delta , 0.0);    
    volumeField delta_star;
    delta_star = delta;

    niter = 0;
    ratio = 1;

    l1norm_old =  l1_norm (delta);

    while (ratio > 1e-8)
    {
        delta_star.mesh = delta.mesh;
        niter+= 1;

        for (int i = 1; i <= N_x/2; i++) {
            for (int j= 1; j <= N_y/2; j++) {

                aP = - 1./ ( (delta_star.mesh[i+1][j].xCentroid() - delta_star.mesh[ i ][j].xCentroid()) * delta_star.mesh[i][j].dx() )
                     - 1./ ( (delta_star.mesh[ i ][j].xCentroid() - delta_star.mesh[i-1][j].xCentroid()) * delta_star.mesh[i][j].dx() )
                     - 1./ ( (delta_star.mesh[i][j+1].yCentroid() - delta_star.mesh[i][ j ].yCentroid()) * delta_star.mesh[i][j].dy() )
                     - 1./ ( (delta_star.mesh[i][ j ].yCentroid() - delta_star.mesh[i][j-1].yCentroid()) * delta_star.mesh[i][j].dy() );

                aW = - 1./ ( (delta_star.mesh[ i ][j].xCentroid() - delta_star.mesh[i-1][j].xCentroid()) * delta_star.mesh[i][j].dx() );
                aE = - 1./ ( (delta_star.mesh[i+1][j].xCentroid() - delta_star.mesh[ i ][j].xCentroid()) * delta_star.mesh[i][j].dx() );

                aS = - 1./ ( (delta_star.mesh[i][ j ].yCentroid() - delta_star.mesh[i][j-1].yCentroid()) * delta_star.mesh[i][j].dy() );
                aN = - 1./ ( (delta_star.mesh[i][j+1].yCentroid() - delta_star.mesh[i][ j ].yCentroid()) * delta_star.mesh[i][j].dy() );

                delta.mesh[i][j].setVal( (1.-omega)*delta_star.mesh[i][j].cellVal() +
                                omega/aP * ( aE*delta_star.mesh[i+1][j].cellVal() + aW*delta_star.mesh[i-1][j].cellVal() + 
                                             aN*delta_star.mesh[i][j+1].cellVal() + aS*delta_star.mesh[i][j-1].cellVal() + res_coarser.mesh[i][j].cellVal()) );
            };
        };
        setBCerrFV(&delta);
        l1norm =  l1_norm (delta);
        ratio = abs(( l1norm -  l1norm_old )) / abs( l1norm_old );
        l1norm_old =  l1norm;

       ///cout << "ratio of successive delta "<< ratio<<"\n";

    };
      cout << "convergence achieved on delta in "<< niter <<" iterations\n";

    
    // prolongate back the error delta on the finer grid //
    volumeField  delta_finer   = buildFVfield(N_x, N_y, Lx, Ly);  
    prolongateFieldVal(delta , &delta_finer);     
 
    // add the the error delta on the smoothed solution //
    *p = *p + delta_finer;

      cout << "error added back on solutions\n";


    niter = 0;
    ratio = 1;
    /// perfom the final post-smoothing iterations
    while (ratio > 1e-7)
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
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal() + (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
        setBCpFV(p);
        l1norm =  l1_norm (*p);
        ratio = abs(( l1norm -  l1norm_old )) / abs( l1norm_old );
       //           cout << "ratio  in post smoothing "<<ratio<<"\n";

        l1norm_old =  l1_norm (p_star);
    };
          cout << "post smoothhing done in "<<niter<<" iterations\n";

    return;
};



////////////////////////////////////////////////////////////////////////////////////////////////////
void multigrid_rec(volumeField* p, const volumeField& div,
            const double tol, const double rho  ,
            const double dt , const double omega,
            const double Lx , const double Ly ,
            int& dimx_coarse, int& dimy_coarse, int& ncycle)
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
    
    //// execute pre-smoother 10 iterations ///    
    while (niter < 20)
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
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal() + (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
        setBCpFV(p);
    };

    cout << "pre-smoothing performed in "<< niter <<" iterations\n";


    //// compute residual ////
    volumeField  res      = buildFVfield(N_x, N_y, Lx, Ly);
    initFieldVal(res , 0.0);

    for (int i = 1; i <= N_x; i++) {
        for (int j= 1; j <= N_y; j++) {

                aP = - 1./ ( (res.mesh[i+1][j].xCentroid() - res.mesh[ i ][j].xCentroid()) * res.mesh[i][j].dx() )
                     - 1./ ( (res.mesh[ i ][j].xCentroid() - res.mesh[i-1][j].xCentroid()) * res.mesh[i][j].dx() )
                     - 1./ ( (res.mesh[i][j+1].yCentroid() - res.mesh[i][ j ].yCentroid()) * res.mesh[i][j].dy() )
                     - 1./ ( (res.mesh[i][ j ].yCentroid() - res.mesh[i][j-1].yCentroid()) * res.mesh[i][j].dy() );

                aW = - 1./ ( (res.mesh[ i ][j].xCentroid() - res.mesh[i-1][j].xCentroid()) * res.mesh[i][j].dx() );
                aE = - 1./ ( (res.mesh[i+1][j].xCentroid() - res.mesh[ i ][j].xCentroid()) * res.mesh[i][j].dx() );

                aS = - 1./ ( (res.mesh[i][ j ].yCentroid() - res.mesh[i][j-1].yCentroid()) * res.mesh[i][j].dy() );
                aN = - 1./ ( (res.mesh[i][j+1].yCentroid() - res.mesh[i][ j ].yCentroid()) * res.mesh[i][j].dy() );

            res.mesh[i][j].setVal( -aE*p->mesh[i+1][j].cellVal() 
                                   -aW*p->mesh[i-1][j].cellVal()  
                                   +aP*p->mesh[i ][j ].cellVal() 
                                   -aN*p->mesh[i][j+1].cellVal()  
                                   -aS*p->mesh[i][j-1].cellVal()  
                                   + (rho/dt)*div.mesh[i][j].cellVal() );
        };
    };
    setBCpFV(&res);

    cout << "residual initialized ..\n";

    dimx_coarse = dimx_coarse/2;
    dimy_coarse = dimy_coarse/2;

    // interpolate res on the coarser grid //
    volumeField  res_coarser = buildFVfield(dimx_coarse, dimy_coarse, Lx, Ly);
    interpolateFieldVal(res , &res_coarser); 

    //// solve for the error delta on the coarser grid ////
    volumeField  delta       = buildFVfield(dimx_coarse, dimy_coarse, Lx, Ly);
    initFieldVal(delta , 1.0);    
    volumeField delta_star;
    delta_star = delta;

    niter = 0;
    ratio = 1;

    l1norm_old =  l1_norm (delta);

    while (ratio > 1e-14)
    {
        delta_star.mesh = delta.mesh;
        niter+= 1;

        for (int i = 1; i <= dimx_coarse; i++) {
            for (int j= 1; j <= dimy_coarse; j++) {

                aP = - 1./ ( (delta_star.mesh[i+1][j].xCentroid() - delta_star.mesh[ i ][j].xCentroid()) * delta_star.mesh[i][j].dx() )
                     - 1./ ( (delta_star.mesh[ i ][j].xCentroid() - delta_star.mesh[i-1][j].xCentroid()) * delta_star.mesh[i][j].dx() )
                     - 1./ ( (delta_star.mesh[i][j+1].yCentroid() - delta_star.mesh[i][ j ].yCentroid()) * delta_star.mesh[i][j].dy() )
                     - 1./ ( (delta_star.mesh[i][ j ].yCentroid() - delta_star.mesh[i][j-1].yCentroid()) * delta_star.mesh[i][j].dy() );

                aW = - 1./ ( (delta_star.mesh[ i ][j].xCentroid() - delta_star.mesh[i-1][j].xCentroid()) * delta_star.mesh[i][j].dx() );
                aE = - 1./ ( (delta_star.mesh[i+1][j].xCentroid() - delta_star.mesh[ i ][j].xCentroid()) * delta_star.mesh[i][j].dx() );

                aS = - 1./ ( (delta_star.mesh[i][ j ].yCentroid() - delta_star.mesh[i][j-1].yCentroid()) * delta_star.mesh[i][j].dy() );
                aN = - 1./ ( (delta_star.mesh[i][j+1].yCentroid() - delta_star.mesh[i][ j ].yCentroid()) * delta_star.mesh[i][j].dy() );

                delta.mesh[i][j].setVal( (1.-omega)*delta_star.mesh[i][j].cellVal() +
                                omega/aP * ( aE*delta_star.mesh[i+1][j].cellVal() + aW*delta_star.mesh[i-1][j].cellVal() + 
                                             aN*delta_star.mesh[i][j+1].cellVal() + aS*delta_star.mesh[i][j-1].cellVal() + res_coarser.mesh[i][j].cellVal()) );
            };
        };
        setBCerrFV(&delta);
        l1norm =  l1_norm (delta);
        ratio = abs(( l1norm -  l1norm_old )) / abs( l1norm_old );
        l1norm_old =  l1norm;

       ///cout << "ratio of successive delta "<< ratio<<"\n";

    };
      cout << "convergence achieved on delta in "<< niter <<" iterations\n";

    
    // prolongate back the error delta on the finer grid //
    volumeField  delta_finer   = buildFVfield(N_x, N_y, Lx, Ly);  
    prolongateFieldVal(delta , &delta_finer);     
 
    // add the the error delta on the smoothed solution //
    *p = *p + delta_finer;

      cout << "error added back on solutions\n";


    niter = 0;
    ratio = 1;
    /// perfom the final post-smoothing iterations
    while (ratio > 1e-10)
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
                                omega/aP * ( aE*p_star.mesh[i+1][j].cellVal() + aW*p_star.mesh[i-1][j].cellVal() + 
                                             aN*p_star.mesh[i][j+1].cellVal() + aS*p_star.mesh[i][j-1].cellVal() + (rho/dt)*div.mesh[i][j].cellVal()) );


            };
        };
        setBCpFV(p);
        l1norm =  l1_norm (*p);
        ratio = abs(( l1norm -  l1norm_old )) / abs( l1norm_old );
       //           cout << "ratio  in post smoothing "<<ratio<<"\n";

        l1norm_old =  l1_norm (p_star);
    };
          cout << "post smoothing done in "<<niter<<" iterations\n";

    return;
};