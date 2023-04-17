#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>      
#include <mpi.h>

#include "FVIO.h"
////////////////////////////////////////////
void readInput(int rank, int& restart, 
                         double& t0, 
                         double& dt, 
                         double& tend, 
                         double& dt_write,
                         double& mu,
                         double& rho,
                         double& Lx,
                         double& Ly)
{
 if (rank == 0)
 {
  FILE* fid;
  fid = fopen("input.in","r");
  if (fid)
  {
    fscanf(fid, "%d\n",  &restart);
    fscanf(fid, "%lf\n", &t0);
    fscanf(fid, "%lf\n", &dt);
    fscanf(fid, "%lf\n", &tend);
    fscanf(fid, "%lf\n", &dt_write);

    fscanf(fid, "%lf\n", &mu);
    fscanf(fid, "%lf\n", &rho);
    fscanf(fid, "%lf\n", &Lx);
    fscanf(fid, "%lf\n", &Ly);  

    std:: cout << "------- Lid Cavity Flow 2D ---------"<<"\n";
    std:: cout << "\n";
    std:: cout << "simulation parameters:"<<"\n";
    std:: cout << "\n";
    std:: cout << "-----------------------------------"<<"\n";
    std:: cout << "tend "<<tend <<"\n";
    std:: cout << "dt "<<dt <<"\n";
    std:: cout << "mu "<<mu <<"\n";
    std:: cout << "rho "<<rho <<"\n";
    std:: cout << "Lx "<<Lx <<"\n";
    std::cout << "Ly "<<Ly <<"\n";
    std:: cout << "-----------------------------------"<<"\n";
    std:: cout << "\n";
  }
  else
  {
    std:: cout << "***** input.in not found *****\n";
    return;
  };
  fclose(fid);
  }
  MPI_Bcast(&restart,  1, MPI_INT,   0,MPI_COMM_WORLD);
  MPI_Bcast(&t0,       1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&dt,       1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&tend,     1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&dt_write, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&mu,       1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&rho,      1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lx,       1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Ly,       1, MPI_DOUBLE,0,MPI_COMM_WORLD);


};
////////////////////////////////////////////
void meshToVTK(FVmesh& mesh, int rank)
{
  int Nx = mesh.get_dimx();
  int Ny = mesh.get_dimy();

  int Nx_tot = mesh.get_dimx_tot();
  int Ny_tot = mesh.get_dimy_tot();

  int nProcs = Nx_tot/Nx;

  std::vector<double> local_mesh_xc;
  std::vector<double> local_mesh_yc;

  std::vector<bool> local_mesh_id; 

  local_mesh_xc.reserve(Nx );
  local_mesh_yc.reserve(Ny);

  local_mesh_id.reserve(Nx*Ny);

  std::vector<double> global_mesh_xc;
  std::vector<double> global_mesh_yc;

  std::vector<bool> global_mesh_id;

  if (rank==0)
  {
   global_mesh_xc.reserve(Nx_tot );
   global_mesh_yc.reserve( Ny_tot);

   global_mesh_id.reserve( (Ny_tot)*(Nx_tot) );
   }

  // assemble local data buffer to be sent to master:
  // only internal nodes of mesh are considered
  // centroids
  for (int i=0;i<Nx;i++)
  {
    local_mesh_xc[i] =mesh.getFVcell(i+1,1).xCentroid() ;
  }

  for (int j=0;j<Ny;j++)
  {
    local_mesh_yc[j] =mesh.getFVcell(1,j+1).yCentroid() ;   
    
  }

  //assemble vector if cells ID for each rank
  for(int i=0; i<Nx; i++)
  {
    for (int j=0; j<Ny; j++)
    {
      local_mesh_id[i*Ny+j] = mesh.getFVcell(i+1,j+1).IDCELL().INTERNAL ;
    }
  }

 
  // send the data buffer to master process
  MPI_Gather(&local_mesh_xc.front(), (Nx)   , MPI_DOUBLE, &global_mesh_xc.front(), (Nx)   , MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Gather(&local_mesh_id.front(), (Nx*Ny), MPI_C_BOOL, &global_mesh_id.front(), (Nx*Ny), MPI_C_BOOL, 0, MPI_COMM_WORLD);

 // MPI_Gather(&local_mesh_xf.front(), (Nx+1)*(Ny+1), MPI_DOUBLE, &global_mesh_xf.front(), (Nx+1)*(Ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Gather(&local_mesh_yf.front(), (Nx+1)*(Ny+1), MPI_DOUBLE, &global_mesh_yf.front(), (Nx+1)*(Ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
 

  if (rank ==0)
  {
 
    std::string file = "mesh.vtk";

    FILE *fid;
    fid = fopen (file.c_str(), "w");
    fprintf (fid, "# vtk DataFile Version 3.0\n");
    fprintf (fid, "Solution file\n");
    fprintf (fid, "ASCII\n");
    fprintf (fid, "DATASET RECTILINEAR_GRID\n");
    fprintf (fid, "DIMENSIONS %d %d %d\n",Nx_tot,Ny_tot,1);

    fprintf (fid,"X_COORDINATES  %d double\n",Nx_tot);
    for ( int i = 0; i < Nx_tot; i++ )
    {

       fprintf(fid,"%16.7e \n",global_mesh_xc[i]);

    };

    fprintf (fid,"Y_COORDINATES  %d double\n",Ny_tot);
    for ( int j = 0; j < Ny_tot; j++ )
     {
       fprintf(fid,"%16.7e \n",local_mesh_yc[j]);
    };

    fprintf (fid,"Z_COORDINATES  %d double\n",1);
    fprintf (fid,"%16.7e \n",1.);
 
    
    fclose(fid);
  }
 }