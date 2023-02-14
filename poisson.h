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
///////////////////////////////////////////////////////////////////
double sum (const volumeField& field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;

    double accum = 0.;

    for (int i = 1; i <= N_x; i++)
    {
         for (int j= 1; j <= N_y; j++)
         {
          accum = accum + ( field.mesh[i][j].cellVal() ) ;
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