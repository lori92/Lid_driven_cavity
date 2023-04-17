#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iomanip>
//#include <mpi.h>

#include <petscksp.h>

#include "global_variables.h"
#include "FVmesh.h"
#include "FVIO.h"
#include "FVmatrix.h"

/// total number of internal cells (without ghost nodes) ///

static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";


int main( int argc, char **args)
{
  
  // Initialize the MPI environment.
  PetscMPIInt size;
  PetscInt  rstart, rend;
  PetscErrorCode ierr;
  int rank;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 4, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!" );
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


  // Define number of grid points (centroids) per MPI process
  int Nx = Nx_tot/num_procs;
  int Ny = Ny_tot;

  // create the grid
  FVmesh grid(Nx, Ny, Nx_tot, Ny_tot,  1., 1., rank);
  grid.writeFVmesh(rank);
  meshToVTK(grid, rank);

  //define the matrix
  Mat         A;     

  // n: local size of rows/colums for each rank
  // N: global size
  int N = (Nx_tot)*(Ny_tot);
  int n = N/num_procs;     
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatSetBlockSize(A, 1));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));
  //PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));

  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));
 // PetscCall(MatGetLocalSize(A, &n, &n));

  std::cout<<"rank "<<rank<<" rstart "<<rstart<<std::endl;
  std::cout<<"rank "<<rank<<" rend "<<rend<<std::endl;  


  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));                      

  PetscCall(assemble_matrix(  A,  grid,   rank, rstart, rend));
  /* Assemble the matrix */
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  PetscViewer viewer; 
  PetscObjectSetName((PetscObject)A, "A");
  /* ... populate matrix A and vector v ... */
  /* use built-in stdout viewer */
  viewer = PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD);
  /* view matrix/vector */
  PetscCall(MatView(A, viewer));


  // read simulation parameters
  readInput(rank, restart, t0, dt, tend,dt_write, mu, rho,Lx, Ly);

  // Finalize the MPI environment.
  PetscCall(PetscFinalize());

  return 0;
}
 