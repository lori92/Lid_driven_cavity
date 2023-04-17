/*
 * write_paraview_output.c
 *
 *  Created on: Feb 6, 2023
 *      Author: Lorenzo Sufrà
 */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>      
#include <mpi.h>


#include "mesh.h"
#include "write_output.h"

 
void write_output(const volumeField& vel_x, const volumeField& vel_y,  const volumeField& p, const double time, std::string namefile)
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

    std::vector<double> x(N_p), y(N_p), z(N_p);
    std::vector<double> u(N), v(N), w(N), pr(N);  

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
            //u[n] = double(vel_x.mesh[i+1][j+1].IDCELL().EAST ); 
            //v[n] = double(vel_x.mesh[i+1][j+1].IDCELL().WEST ); 
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

void write_data(const volumeField& vel_x, const volumeField& vel_y,  const volumeField& p, const double time) 
{
  int N_x = vel_x.dimx;
  int N_y = vel_y.dimy;

  int n   = N_x*N_y;

  FILE* fid = fopen("data.dat", "w");
  
  if (fid)
  {
   fprintf(fid, "%16.7e\n", time);
   
   for ( int j = 0; j <= N_y; j++ )
   {
    for ( int i = 0; i <= N_x; i++ )
    {
      fprintf(fid, "%16.7e  %16.7e  %16.7e\n", 
                   vel_x.mesh[i][j].cellVal(), vel_y.mesh[i][j].cellVal(), p.mesh[i][j].cellVal() ); 
   

            n++;
    }
   }
  }
  else
   {
    std::cout << "error while writing data at time "<< time <<"\n";
   };
     fclose(fid);


 FILE *OutFile_meshx = fopen("x.dat","w++");
 
   for(int i=1; i<=N_x; i++){
       double x = (vel_x.mesh[i][0].xFace(1) );
      fprintf(OutFile_meshx,"\n %f ", x); }  

fclose(OutFile_meshx);

 FILE *OutFile_meshy = fopen("y.dat","w++");
 
   for(int j=1; j<=N_y; j++){
       double y = (vel_x.mesh[0][j].yFace(1) );
      fprintf(OutFile_meshy,"\n %f ", y); }  

fclose(OutFile_meshy);

return;
}


void read_data(volumeField* vel_x,  volumeField* vel_y,   volumeField* p, double& time) 
{
  int N_x = vel_x->dimx;
  int N_y = vel_y->dimy;

  int n   = N_x*N_y;

  FILE* fid = fopen("data.dat", "r");

  double velx, vely, pr;
  
  if (fid)
  {
   fscanf(fid, "%lf\n", &time);
   
   for ( int j = 0; j <= N_y; j++ )
   {
    for ( int i = 0; i <= N_x; i++ )
    {
      fscanf(fid, "%lf  %lf  %lf\n", 
                   &velx, &vely, &pr ); 

      vel_x->mesh[i][j].setVal(velx);             
      vel_y->mesh[i][j].setVal(vely);             
      p->mesh[i][j].setVal(pr);             
   

       n++;
    }
   }
  }


   else
   {
    std::cout << "error while writing data at time "<< time <<"\n";
   };
       fclose(fid);
return;
}

void write_k(const volumeField vel_x, const volumeField vel_y, const double time) 
{

 double k = 0;
 double totalK = 0;

 double cellVol;
 double totalCellVol = 0;

 int N_x = vel_x.dimx;
 int N_y = vel_x.dimy;

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
        k =  0.5* ( pow( (vel_x.mesh[i][j].cellVal() + vel_x.mesh[i-1][j].cellVal()) ,2) +  
                       pow( (vel_x.mesh[i][j].cellVal() + vel_x.mesh[i][j-1].cellVal()) ,2) ); 
         cellVol = vel_x.mesh[i][j].dx() * vel_x.mesh[i][j].dy() ;
         totalCellVol = totalCellVol + cellVol;

         totalK = totalK + k * cellVol;
  };
};
  totalK = totalK / totalCellVol;

  std::ofstream outfile;
  outfile.open("k.dat", std::ios_base::app); // append instead of overwrite
  outfile << std::scientific << std::setprecision(5) << time  <<"    "<< std::scientific << std::setprecision(5) << totalK <<"\n"; 

return;
}

void write_centre_line(const volumeField vel_x, const volumeField vel_y, const double time) 
{

 double k = 0;
 double totalK = 0;

 double cellVol;
 double totalCellVol = 0;

 int N_x = vel_x.dimx;
 int N_y = vel_x.dimy;

 std::ofstream outfile;
 outfile.open("u.dat", std::ios_base::out); // append instead of overwrite
 outfile << std::scientific << std::setprecision(5) << time  <<"\n"; 

 for (int j = 1; j<=N_y; j++) 
 {
  outfile << std::scientific << std::setprecision(5) << vel_x.mesh[N_x/2][j].cellVal() <<"    "<< std::scientific << std::setprecision(5) << vel_x.mesh[N_x/2][j].yCentroid()  <<"\n"; 
 };

 std::ofstream outfile_v;
 outfile_v.open("v.dat", std::ios_base::out); // append instead of overwrite
 outfile_v << std::scientific << std::setprecision(5) << "#" << time  <<"\n"; 

 for (int i = 1; i<=N_x;i++) 
 {
  outfile_v << std::scientific << std::setprecision(5) << vel_y.mesh[i][N_y/2].cellVal()  <<"    "<< std::scientific << std::setprecision(5) << vel_y.mesh[i][N_y/2].xCentroid()  <<"\n"; 
 };

return;
}

void read_input(const int rank, int& restart, double &t0, double& tend, double& dt_write, double& mu, double& rho, double& Lx, double& Ly)
{

 MPI_Comm comm;
 
 if (rank == 0)
 {
   FILE* fid;
   fid = fopen("input.in","r");
   if (fid)
   {
      fscanf(fid, "%d\n",  &restart);
      fscanf(fid, "%lf\n", &t0);
      fscanf(fid, "%lf\n", &tend);
      fscanf(fid, "%lf\n", &dt_write);

      fscanf(fid, "%lf\n", &mu);
      fscanf(fid, "%lf\n", &rho);
      fscanf(fid, "%lf\n", &Lx);
      fscanf(fid, "%lf\n", &Ly);  
   }
   else
   {
      std::cout << "***** input.in not found *****\n";
      return;
   };
   fclose(fid);
 }

 MPI_Bcast( &restart, 1, MPI_INT    , 0, MPI_COMM_WORLD );
 MPI_Bcast( &t0, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &tend, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &dt_write, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &mu, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &rho, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &Lx, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );
 MPI_Bcast( &Ly, 1, MPI_DOUBLE , 0, MPI_COMM_WORLD );





return;
};