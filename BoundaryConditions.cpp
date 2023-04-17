#include <iostream>
#include <cmath>
#include <petscksp.h>

#include "mesh.h"

void shift(const int rank, const int nprocs, volumeField* u)
{
 const int Nx = u->dimx;
 const int Ny = u->dimy;

 double u_send_left [Ny+2]; 
 double u_send_right[Ny+2]; 

 double u_recv_left [Ny+2]; 
 double u_recv_right[Ny+2]; 

 int dest;
 int send;
 int tag;


 MPI_Status status;


 switch(rank)
 {

   // first rank (master) sends to right and receive from right 
   case(0):
   {
    for (int j = 0; j<=Ny+1; j++)
    {
      u_send_right[j]=u->mesh[Nx][j].cellVal();  
    }
    tag = 1;

    dest = rank+1;
    MPI_Send ( u_send_right, Ny+2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD );

    send = rank + 1;
    MPI_Recv ( u_recv_right, Ny+2, MPI_FLOAT, send, tag, MPI_COMM_WORLD, &status );
    for (int j = 0; j<=Ny+1; j++)
    {
      u->mesh[Nx+1][j].setVal( u_recv_right[j] );  
    }
   break;
   };

   // last rank (master) sends to left (rank-1) and receives from left (rank-1)
   case(4):
   {
    for (int j = 0; j<=Ny+1; j++)
    {
      u_send_left[j]=u->mesh[1][j].cellVal();  
    }
    tag = 1;

    dest = rank-1;
    MPI_Send ( u_send_left, Ny+2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD );

    send = rank-1;
    MPI_Recv ( u_recv_left, Ny+2, MPI_FLOAT, send, tag, MPI_COMM_WORLD, &status );
    for (int j = 0; j<=Ny+1; j++)
    {
      u->mesh[0][j].setVal( u_recv_left[j] );  
    }
   break;
   };   


   // internal rank sends to left (rank-1) and right (rank+1) and receive from left (rank-1) and right (rank+1)
   default:
   {
    for (int j = 0; j<=Ny+1; j++)
    {
      u_send_left [j]=u->mesh[1 ][j].cellVal();  
      u_send_right[j]=u->mesh[Nx][j].cellVal();  
    }
    tag = 1;

    // sending to left 
    dest = rank-1;
    MPI_Send ( u_send_left, Ny+2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD );

    // sending to right 
    dest = rank+1;
    MPI_Send ( u_send_right, Ny+2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD );

    // receiving from left
    send = rank-1;
    MPI_Recv ( u_send_left, Ny+2, MPI_FLOAT, send, tag, MPI_COMM_WORLD, &status );

    // receiving from right
    send = rank+1;
    MPI_Recv ( u_send_right, Ny+2, MPI_FLOAT, send, tag, MPI_COMM_WORLD, &status );
    for (int j = 0; j<=Ny+1; j++)
    {
      u->mesh[  0 ][j].setVal(u_recv_left [j]);  
      u->mesh[Nx+1][j].setVal(u_recv_right[j]);  
    }
   break;
   };      


 };

 return;
}
