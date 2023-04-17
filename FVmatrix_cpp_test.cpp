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

    std::cout << "cell n = "<<n<<", INTERNAL = "<< mesh.getFVcell(i,j).IDCELL().INTERNAL<<std::endl;


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























#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>      
#include <mpi.h>

#include "FVmatrix.h"

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
////////////////////////////////////////////////////////////////////////////////////////////////
 for (int i=1; i<=Nx; i++){
   for (int j=1; j<=Ny; j++){

    i_global  = (i) +rank*Nx;

    n = (i_global-1)*Ny + (j);
    index_sud = n-1;
    index_nord= n+1;
    index_west= n-Ny;
    index_east= n+Ny;

     std:: cout << "for rank ="<<rank<<",  node nn= "<<nn<<"ID_south= "<< mesh.getFVcell(i,j).IDCELL().SUD<< std::endl;

    // if n corresponds to an internal cell
    if (n-1<rend && n-1>=rstart && mesh.getFVcell(i,j).IDCELL().INTERNAL ==1)
    {
      /*
      A_[ n ][ n ] = aP[i][j];
      A_[ n ][n-1] = aS[i][j];
      A_[ n ][n+1] = aN[i][j];
      A_[ n ][n-N_y]= aW[i][j];
      A_[ n ][n+N_y]= aE[i][j];
     */

      nn = n-1;
      index_sud = nn-1;
      index_nord= nn+1;
      index_west= nn-Ny;
      index_east= nn+Ny; 
          
      ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &aP[i][j], INSERT_VALUES));
      ierr =(MatSetValues(A, 1, &nn, 1, &index_sud , &aS[i][j], INSERT_VALUES));
      ierr =(MatSetValues(A, 1, &nn, 1, &index_nord, &aN[i][j], INSERT_VALUES));
      ierr =(MatSetValues(A, 1, &nn, 1, &index_west, &aW[i][j], INSERT_VALUES));
      ierr =(MatSetValues(A, 1, &nn, 1, &index_east, &aE[i][j], INSERT_VALUES));

      std:: cout << "for rank ="<<rank<<", internal node nn= "<<nn<<std::endl;


    } 
    
    // if n does not correspond to an internal cell: check on which boundary I am
    else if (n-1<rend && n-1>=rstart && mesh.getFVcell(i,j).IDCELL().INTERNAL ==0)
    {
       nn = n-1;
       index_sud = nn-1;
      index_nord= nn+1;
      index_west= nn-Ny;
      index_east= nn+Ny; 
      //A_[ n ][ n ] = aP[i][j];
      ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &aP[i][j], INSERT_VALUES));

      if (mesh.getFVcell(i,j).IDCELL().NORD ==1)   
      {
       //A_[ n ][ n ] = A_[ n ][ n ] + aN[i][j] ; 
       PetscScalar val =   aN[i][j];
        //ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &val, ADD_VALUES));
      }

      if (mesh.getFVcell(i,j).IDCELL().SUD ==1)  
      { 
          //    std:: cout << "for rank ="<<rank<<", south node nn= "<<nn<<std::endl;
              

       //A_[ n ][ n ] = A_[ n ][ n ] + aS[i][j];
       PetscScalar val =  aS[i][j];
        //ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &val, ADD_VALUES));    
          }

      if (mesh.getFVcell(i,j).IDCELL().EAST ==1)  
      { 
       //A_[ n ][ n ] = A_[ n ][ n ] + aE[i][j];
       PetscScalar val =   aE[i][j];
        //ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &val, ADD_VALUES));     
         }

      if (mesh.getFVcell(i,j).IDCELL().WEST ==1)  
      { 
       //A_[ n ][ n ] = A_[ n ][ n ] + aW[i][j];
       PetscScalar val =   aW[i][j];
       //ierr =(MatSetValues(A, 1, &nn, 1, &nn         , &val, ADD_VALUES)); 
            }

      // assemble anyway the ghost slots, which will
      // be discarded when the final matrix will be
      // obtained by restricting A_ 
      
      if (!mesh.getFVcell(i,j).IDCELL().SUD  )
      { 
        //A_[ n ][n -  1]= (aS[i][j]); 
        ierr =(MatSetValues(A, 1, &nn, 1, &index_sud         , &aS[i][j], INSERT_VALUES));
      }
      if (!mesh.getFVcell(i,j).IDCELL().NORD )
      {
        // A_[ n ][n +  1]= (aN[i][j]); 
        ierr =(MatSetValues(A, 1, &nn, 1, &index_nord         , &aN[i][j], INSERT_VALUES));
      }
      if (!mesh.getFVcell(i,j).IDCELL().WEST )
      {
         //A_[ n ][n -N_y]= (aW[i][j]); 
        ierr =(MatSetValues(A, 1, &nn, 1, &index_west         , &aW[i][j], INSERT_VALUES));
      }
      if (!mesh.getFVcell(i,j).IDCELL().EAST )
      {
        // A_[ n ][n +N_y]= (aE[i][j]); 
        ierr =(MatSetValues(A, 1, &nn, 1, &index_east         , &aE[i][j], INSERT_VALUES));
      }
	}

  }
 }
/*
 n          = N-1;
 index_nord = n + 1;

 if (n+1<rend && n+1>=rstart)
 {
   //A_[ n ][ n ] = A_[ n ][ n ] - aN[N_x][N_y];
   PetscScalar val = aP[Nx][Ny] - aN[Nx][Ny];
   ierr =(MatSetValues(A, 1, &n, 1, &n         ,  &val, INSERT_VALUES));

   //A_[ n ][n+1] = aN[N_x][N_y];
   ierr =(MatSetValues(A, 1, &n, 1, &index_nord         ,  &aN[Nx][Ny], INSERT_VALUES));
 }

 n          = N;
 index_sud = n - 1;

 if (n<rend && n>=rstart)
 {
   // A_[N][N] = 1.;
   PetscScalar val = 1.0;
   ierr =(MatSetValues(A, 1, &n, 1, &n         ,  &val, INSERT_VALUES));

   // A_[N][N-1] = 1.;
   ierr =(MatSetValues(A, 1, &n, 1, &index_sud ,  &val, INSERT_VALUES));
 }
*/
////////////////////////////////////////////////////////////////////////////////////////////////
  return ierr;
}