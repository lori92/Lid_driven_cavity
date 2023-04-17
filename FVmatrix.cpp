#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>      
#include <mpi.h>

#include "FVmatrix.h"

inline int return_i(int n,  int N_y)
{
  return (n-1)/N_y+1;
}
 
inline int return_j(int i, int n, int N_y)
{
  return n - (i-1)*N_y;
}


PetscErrorCode assemble_matrix(Mat A,  FVmesh& mesh, int rank, PetscInt rstart, PetscInt rend)
{

 PetscErrorCode ierr;  

 // number of physical nodes for the local rank
 int Nx = mesh.get_dimx();
 int Ny = mesh.get_dimy();

 int Nx_tot = mesh.get_dimx_tot();
 int Ny_tot = mesh.get_dimy_tot();

 int nProcs = Nx_tot/Nx;
 int N = Nx_tot*Ny_tot;


 int i_global;

 // declare local indexes of mesh points
 int n, nn;
 int index_sud  ;
 int index_nord ;
 int index_west ;
 int index_east ;

 // discretization coefficients for ddp/dx2 in (i,j) as 
 // ddp/dx2 = aP*P + aP_W*P_W + aE*P_E + aN*P_N + aS*P_S
 PetscScalar aP[Nx][Ny];
 PetscScalar aS[Nx][Ny];
 PetscScalar aN[Nx][Ny];
 PetscScalar aE[Nx][Ny];
 PetscScalar aW[Nx][Ny];

 for (int i = 1; i<=Nx; i++)
 {
  for (int j = 1; j<=Ny; j++)
  {
     aP[i][j] = - 1./ ( (mesh.getFVcell(i+1,j).xCentroid() - mesh.getFVcell( i ,j).xCentroid()) * 
                         mesh.getFVcell(i,j).dx() )
                - 1./ ( (mesh.getFVcell( i ,j).xCentroid() - mesh.getFVcell(i-1,j).xCentroid()) * 
                         mesh.getFVcell(i,j).dx() )
                - 1./ ( (mesh.getFVcell(i,j+1).yCentroid() - mesh.getFVcell(i, j ).yCentroid()) * 
                         mesh.getFVcell(i,j).dy() )
                - 1./ ( (mesh.getFVcell(i, j ).yCentroid() - mesh.getFVcell(i,j-1).yCentroid()) * 
                         mesh.getFVcell(i,j).dy() );

     aW[i][j] =  1./ ( (mesh.getFVcell( i ,j).xCentroid() - mesh.getFVcell(i-1,j).xCentroid()) * 
                        mesh.getFVcell(i,j).dx() );
     aE[i][j] =  1./ ( (mesh.getFVcell(i+1,j).xCentroid() - mesh.getFVcell( i ,j).xCentroid()) * 
                        mesh.getFVcell(i,j).dx() );

     aS[i][j] =  1./ ( (mesh.getFVcell(i, j ).yCentroid() - mesh.getFVcell(i,j-1).yCentroid()) * 
                        mesh.getFVcell(i,j).dy() );
     aN[i][j] =  1./ ( (mesh.getFVcell(i,j+1).yCentroid() - mesh.getFVcell(i, j ).yCentroid()) * 
                        mesh.getFVcell(i,j).dy() );
  }
 }
/////////////////////////////////////////////////////////////////////////////////
 int i;
 int j;

// nn row counter starting from rstart to rend
 for (int nn=rstart; nn<rend; nn++)
 {
    // nn global mesh counter 
    n = nn+1;
    i = return_i(n,Ny);
    j = return_j(i,n,Ny);

    //std::cout << "cell n = "<<n<<", INTERNAL = "<< mesh.getFVcell(i,j).IDCELL().INTERNAL<<std::endl;


    // if n corresponds to an internal cell
     if (mesh.getFVcell(i,j).IDCELL().INTERNAL)
     {
      /*  
      A_[ n ][ n ] = aP[i][j];
      A_[ n ][n-1] = aS[i][j];
      A_[ n ][n+1] = aN[i][j];
      A_[ n ][n-N_y]= aW[i][j];
      A_[ n ][n+N_y]= aE[i][j];
      */
        //  std::cout << "internal cell n = "<<n<<", (i,j) = ("<< i<<","<<j <<") on rank = "<<rank<<std::endl;

       //   std::cout << "aP[i][j] = "<<aP[i][j]<<", (i,j) = ("<< i<<","<<j <<") on rank = "<<rank<<std::endl;

      ierr = (MatSetValues(A, 1, &nn, 1, &nn, &aP[i][j], INSERT_VALUES));
     }   

 }
////////////////////////////////////////////////////////////////////////////////////////////////
  return ierr;
}