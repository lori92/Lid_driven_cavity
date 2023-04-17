#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iomanip>

#include "FVmesh.h"

///////////////////////////////////////////////////////////////
void FVcell :: setCellPos(double x, double y, double dx, double dy)
{
	// x, z coordinates of the cell center of a given cell corresponding to index i, j
        xC = x; 
        xFaces[0] = xC - dx/2.;
        xFaces[1] = xC + dx/2.;
 
        yC = y; 
	yFaces[0] = yC - dy/2.;
	yFaces[1] = yC + dy/2.;
        	
};      

void FVcell :: setCellLabel(int i, int j, int N_x, int N_y, int rank, int num_procs)
{
  // setting the ID for a given cell according to its indexing (i,j)

  // default: internal node
  this->ID.INTERNAL = 1;
  this->ID.EAST     = 0;
  this->ID.WEST     = 0;
  this->ID.NORD     = 0;
  this->ID.SUD      = 0;

  if (i <= 1 && rank ==0) { // west node
    this->ID.INTERNAL = 0;
    this->ID.WEST     = 1;
  }

  if (i >= N_x && rank ==num_procs-1) { // east node
    this->ID.INTERNAL = 0;
    this->ID.EAST     = 1;
  }  

  if (j >= N_y) { // nord node
    this->ID.INTERNAL = 0;
    this->ID.NORD     = 1;
  }  

  if (j <= 1) { // sud node
    this->ID.INTERNAL = 0;
    this->ID.SUD      = 1;
  }
};

ID_CELL FVcell:: IDCELL()
{
  return this->ID;
}

void FVcell :: print_cellCentre(){
        std::cout << "(" <<xC<<","<<yC<<") ";};

double FVcell :: dx(){
       return this->xFaces[1]-this->xFaces[0];};
       
double FVcell :: dy(){
       return this->yFaces[1]-this->yFaces[0];};

double FVcell :: xFace(int face ){
       return this->xFaces[face];};

double FVcell :: yFace(int face){
       return this->yFaces[face];};

double FVcell :: xCentroid(){
       return this->xC;};

double FVcell :: yCentroid(){
       return this->yC;};

///////////////////////////////////////////////////////////////
FVmesh::FVmesh(int Nx, int Ny, int Nx_tot, int Ny_tot, double Lx, double Ly, int rank)
   {      
      double xf[Nx+2];
      double xc[Nx+2];
      double yf[Ny+2];
      double yc[Ny+2];

      double gamma= 1;
      
      this->dimx = Nx;
      this->dimy = Ny;
      this->dimx_tot = Nx_tot;
      this->dimy_tot = Ny_tot;
      this->Lx = Lx;
      this->Ly = Ly;
      this->mesh.reserve( (Nx+2)*(Ny+2) );
        
      gamma  = 0.0000000001 ;

      int nProcs  = Nx_tot /  Nx;
      int Lx_rank = Lx     / (nProcs);

      /// assemble on x the faces of cell
      /// 0 is the left one on the wall
      /// Nx is the right one on the other wall
       for (int i = 0; i <= Nx; i++)
      {  
        xf[i] = Lx * 0.5*(1. - tanh(gamma*(1.-2.*(rank*Nx+i)/(Nx_tot)))/tanh(gamma));
       }

      for (int i = 1; i <= Nx; i++)
      {  
       xc[i] =  xf[i-1] + 0.5*(xf[i]-xf[i-1]);
      }

       double xf_prev;
       xf_prev = Lx * 0.5*(1. - tanh(gamma*(1.-2.*(rank*Nx-1)/(Nx_tot)))/tanh(gamma));
       xc[0]   =  xf[0]- 0.5*(xf[0]-xf_prev);  
    
       double xf_next;
       xf_next = Lx * 0.5*(1. - tanh(gamma*(1.-2.*(rank*Nx+Nx+1)/(Nx_tot)))/tanh(gamma));   
       xc[Nx+1]=  xf[Nx] + 0.5*(xf_next-xf[Nx]); 

       xf[Nx+1]= xf_next;

      for (int j = 0; j <= Ny; j++)
      {
       yf[j] = Ly * 0.5*(1. - tanh(gamma*(1.-2.*j/(Ny_tot)))/tanh(gamma));
      }

      yf[Ny+1] = Ly + (Ly - yf[Ny-1]);
  

      for (int j = 1; j <= Ny; j++)
      {  
       yc[j] =  yf[j-1] + 0.5*(yf[j]-yf[j-1]);
      }

      yc[0]   =  0 - yc[1];
      yc[Ny+1]=  Ly + (Ly-yc[Ny]);

      for (int i = 0; i <= Nx+1; i++)
  {       
   for(int j = 0; j<= Ny+1; j++)
   {       
    double dx = 2.*(xf[i]  - xc[i]);   
    double dy = 2.*(yf[j]  - yc[j]);     
  
    this->mesh[i*(Ny+2) + j].setCellPos( xc[i], yc[j], dx, dy  );
    this->mesh[i*(Ny+2) + j].setCellLabel( i, j, Nx, Ny, rank,  nProcs);
    
    if (rank ==3)
    {
    std::cout << "----------"<<std::endl;
    std::cout << "cell at (i,j) = "<<i<<","<<j<<std::endl;
    std::cout << this->mesh[i*(Ny+2) + j].IDCELL().INTERNAL<<std::endl;
    std::cout << this->mesh[i*(Ny+2) + j].IDCELL().SUD     <<std::endl;
    std::cout << this->mesh[i*(Ny+2) + j].IDCELL().NORD    <<std::endl;
    std::cout << this->mesh[i*(Ny+2) + j].IDCELL().EAST    <<std::endl;
    std::cout << this->mesh[i*(Ny+2) + j].IDCELL().WEST    <<std::endl;
    std::cout << "----------"<<std::endl;
    }

   }     
 }
}; 
///////////////////////////////////////////////////////////////
std::vector<FVcell> FVmesh::mesh;
///////////////////////////////////////////////////////////////
FVcell FVmesh:: getFVcell(int i, int j)
{
  return this->mesh[i*(this->dimy+2) + j ];
};
///////////////////////////////////////////////////////////////
int FVmesh::get_dimx()
{
  return this->dimx;
};
///////////////////////////////////////////////////////////////
int FVmesh::get_dimy()
{
  return this->dimy;
};
///////////////////////////////////////////////////////////////
int FVmesh::get_dimx_tot()
{
  return this->dimx_tot;
};
///////////////////////////////////////////////////////////////
int FVmesh::get_dimy_tot()
{
  return this->dimy_tot;
}
///////////////////////////////////////////////////////////////
void FVmesh :: writeFVmesh(int rank)
{
 int rank_dimx = this->dimx;
 int rank_dimy = this->dimy;

 std:: string filename = std::to_string(rank);
 filename = filename+"_x.dat";

 std::ofstream file;
 file.open(filename);
 int j =0;

 file << std::setw(8)<<"# x_f@left  "<< std::setw(8)<<"x centroid   "<< std::setw(8)<< "x_f@right"<<std::endl;
 for (int i = 0 ; i<=rank_dimx+1; i++)
 {  
  file << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[i*(dimy+2)+j].xFace(0)   <<"\t"
       << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[i*(dimy+2)+j].xCentroid()<<"\t"
       << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[i*(dimy+2)+j].xFace(1)   <<  std::endl;
 
 }
 file.close();



 std:: string filename2 = std::to_string(rank);
 filename2 = filename2+"_y.dat";

 std::ofstream file2;
 file2.open(filename2);

 file2 << std::setw(8)<<"# y_f@down  "<< std::setw(8)<<"y centroid   "<< std::setw(8)<< "y_f@up"<<std::endl;
 for (j = 0 ; j<=rank_dimy+1; j++)
 {  
  file2 << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[3*(dimy+2)+j].yFace(0)   <<"\t"
       << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[3*(dimy+2)+j].yCentroid()<<"\t"
       << std::setw(8)<<std::setprecision(5) << std::left<<std::scientific<<this->mesh[3*(dimy+2)+j].yFace(1)   <<  std::endl;
 
 }
 file2.close();





}

///////////////////////////////////////////////////////////////
void FVfield :: setVal(int i, int j, double val)
{
 cellValue[i*dimy + j] = val;
};
 
void FVfield :: allocateFVfield()
{
 cellValue.reserve( ( dimx+2)*( dimy+2)  );
};
///////////////////////////////////////////////////////////////

