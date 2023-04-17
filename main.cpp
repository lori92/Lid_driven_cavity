#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "global_variables.h"
#include "mesh.h"
#include "poisson_solver.h"
#include "finiteVolume.h"
#include "numerics.h"
#include "predictor_step_new.h"
#include "write_output.h"
#include "BoundaryConditions.h"
#include <chrono>

#include <mpi.h>

using namespace std;

/// total number of internal cells (without ghost nodes) ///
#define Nx_tot 64
#define Ny_tot 64

   


int main( int argc, char *argv[])
{
  int rank, num_procs;
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank(  MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
  MPI_Status status;
  MPI_Request request;

  printf("rank= %d\n", rank);

  int Nx = Nx_tot/num_procs;
  int Ny = Ny_tot;



