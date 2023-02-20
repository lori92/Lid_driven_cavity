/*
 * write_paraview_output.c
 *
 *  Created on: Feb 6, 2023
 *      Author: Lorenzo Sufrà
 */

#include <stdio.h>
#include <stdlib.h>
//#include "finiteVolume.C"

void write_output(const volumeField& vel_x, const volumeField& vel_y,  const volumeField& p, std::string namefile)
{

    FILE *fid;
    int Nx = vel_x.dimx;
    int Ny = vel_x.dimy; 
    int Nz = 1;  
/*
    int N = Nx*Ny*Nz;

    vector<double> x(N), y(N), z(N);
    vector<double> u(N), v(N), w(N), pr(N);  

   int n = 0;

   for ( int k = 0; k < Nz; k++ )
   {
      for ( int j = 0; j < Ny; j++ )
      {
         for ( int i = 0; i < Nx; i++ )
         {
            x[n] = vel_x.mesh[i+1][j+1][k+1].xCentroid();
            y[n] = vel_x.mesh[i+1][j+1][k+1].yCentroid();
            z[n] = vel_x.mesh[i+1][j+1][k+1].zCentroid();

            u[n] = 0.5 * (vel_x.mesh[i+1][j+1][k+1].cellVal() + vel_x.mesh[i][j+1][k+1].cellVal()); 
            v[n] = 0.5 * (vel_y.mesh[i+1][j+1][k+1].cellVal() + vel_y.mesh[i+1][j][k+1].cellVal());
            w[n] = 0.5 * (vel_z.mesh[i+1][j+1][k+1].cellVal() + vel_z.mesh[i+1][j+1][k].cellVal());

            pr[n] = p.mesh[i+1][j+1][k+1].cellVal();
            n++;
         }
      }
   }
*/

    int N_p = (Nx+1)*(Ny+1);
    int N   = (Nx)*(Ny);

    vector<double> x(N_p), y(N_p), z(N_p);
    vector<double> u(N), v(N), w(N), pr(N);  

    int n = 0;


      for ( int j = 0; j <= Ny; j++ )
      {
         for ( int i = 0; i <= Nx; i++ )
         {
            x[n] = vel_x.mesh[i+1][j+1].xFace(0);
            y[n] = vel_x.mesh[i+1][j+1].yFace(0);
            z[n] = 1; //vel_x.mesh[i+1][j+1][k+1].zFace(0);

            //u[n] = 0.5 * (vel_x.mesh[i+1][j+1][k+1].cellVal() + vel_x.mesh[i][j+1][k+1].cellVal()); 
            //v[n] = 0.5 * (vel_y.mesh[i+1][j+1][k+1].cellVal() + vel_y.mesh[i+1][j][k+1].cellVal());
            //w[n] = 0.5 * (vel_z.mesh[i+1][j+1][k+1].cellVal() + vel_z.mesh[i+1][j+1][k].cellVal());

            //pr[n] = p.mesh[i+1][j+1][k+1].cellVal();
            n++;
         }
      }


     n = 0;


      for ( int j = 0; j < Ny; j++ )
      {
         for ( int i = 0; i < Nx; i++ )
         {
            u[n] = 0.5 * (vel_x.mesh[i+1][j+1].cellVal() + vel_x.mesh[i][j+1].cellVal()); 
            v[n] = 0.5 * (vel_y.mesh[i+1][j+1].cellVal() + vel_y.mesh[i+1][j].cellVal());
            //w[n] = 0.5 * (vel_z.mesh[i+1][j+1][k+1].cellVal() + vel_z.mesh[i+1][j+1][k].cellVal());

            pr[n] = p.mesh[i+1][j+1].cellVal();
            n++;
         }
      }




    std::string file = namefile + ".vtk";
    fid = fopen (file.c_str(), "w");
    fprintf (fid, "# vtk DataFile Version 3.0\n");
    fprintf (fid, "Solution file\n");
    fprintf (fid, "ASCII\n");
    fprintf (fid, "DATASET STRUCTURED_GRID\n");
    fprintf (fid, "DIMENSIONS %d %d %d\n",Nx+1,Ny+1,1);

    fprintf (fid,"POINTS %d double\n",N_p);
    for ( int i = 0; i < N_p; i++ )
    {
        fprintf(fid,"%16.7e %16.7e %16.7e\n",x[i], 
                                             y[i],1.);
    };

  //  fprintf(fid,"POINT_DATA %d\n",(Nx)*(Ny)*(Nz));

    fprintf(fid,"CELL_DATA %d\n",(Nx)*(Ny));

    // Write VEL to the file---------------------------
    fprintf (fid,"VECTORS VEL double\n");


    for (int i=0; i<N; i++)
    {
            fprintf (fid,"%16.7e %16.7e %16.7e\n", u[i],
                                                   v[i], 0.0);
    };

    
    // Write PRES to the file---------------------------
    fprintf (fid,"SCALARS PRES double 1\n");
    fprintf (fid,"LOOKUP_TABLE default\n");

    for (int i=0; i<N; i++)
    {
              fprintf (fid,"%16.7e\n", pr[i] );
    };
    
    fclose(fid);

}