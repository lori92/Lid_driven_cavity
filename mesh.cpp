#include <iostream>
using namespace std;

class cell
{       
 private:
  double xC;
  double yC;
  double xFaces[2];
  double yFaces[2];
  double value;

 public:

 void setCellPos   (double, double, double, double);
 void setVal (double);
 void print_cellCentre ();
 void print_cellVal ();
 double cellVal ();
 double dx ();
 double dy ();



};

void cell :: setCellPos(double x, double y, double dx, double dy)
{
	// x, z coordinates of the cell center of a given cell corresponding to index i, j
        xC = x; 
        xFaces[0] = xC - dx/2;
        xFaces[1] = xC + dx/2;
 
        yC = y; 
	yFaces[0] = yC - dy/2;
	yFaces[1] = yC + dy/2;
        	
};      

void cell :: setVal(double value_in)
{
        //  value of the volumeField for given cell corresponding to index i, j
	value = value_in;
};


void cell :: print_cellCentre(){
        cout << "(" <<xC<<","<<yC<<") ";};

void cell :: print_cellVal(){
        cout << value<<" ";};

double cell :: cellVal(){
       return this->value;};
       
double cell :: dx(){
       return this->xFaces[1]-this->xFaces[0];};
       
double cell :: dy(){
       return this->yFaces[1]-this->yFaces[0];};
////////////////////////////////////////////////////////////////
cell** buildMesh (int N_x, int N_y, double L_x, double L_y)
{
 // index i spans the y-direction
 // index j spans the x-direction

 cell **mesh_p = new cell* [N_x+2];
 for (int i = 0; i <= N_x+1; i++)
 {       
   mesh_p[i] = new cell [N_y+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
    mesh_p[i][j].setCellPos( L_y/N_y/2. +  (L_y/N_y)*(i-1), L_x/N_x/2. + (L_x/N_x)*(j-1), L_y/N_y, L_x/N_x  );
   }     
 }       

 return mesh_p;
};

////////////////////////////////////////////////////////////////
void initFieldVal (cell** mesh_p, int N_x, int N_y, double val_in)
{

 for (int i = 0; i <= N_x+1; i++)
 {
   for(int j = 0; j<= N_y+1; j++)
   {
    mesh_p[i][j].setVal(val_in);
   }
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCu (cell** mesh_p, int N_x, int N_y)
{
 // set wall boundary condition for staggered Ux velocity field  

// along y-axis     
 for (int j= 0; j <= N_y+1; j++)
 {
    mesh_p[ 0   ][j].setVal(0);
    mesh_p[N_x  ][j].setVal(0);
    mesh_p[N_x+1][j].setVal(0);
 }

// along x-axis     
 for (int i = 0; i <= N_x+1; i++)
 {
    mesh_p[i][  0  ].setVal( -mesh_p[i][ 1 ].cellVal() )  ;
    mesh_p[i][N_x+1].setVal( -mesh_p[i][N_x].cellVal() )  ;
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCv (cell** mesh_p, int N_x, int N_y)
{
 // set wall boundary condition for staggered Uy velocity field  

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 { 
    mesh_p[  0  ][j].setVal( - mesh_p[ 1 ][j].cellVal()  );
    mesh_p[N_x+1][j].setVal( - mesh_p[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    mesh_p[i][ N_y ].setVal(0);
    mesh_p[i][N_y+1].setVal(0);
    mesh_p[i][  0  ].setVal(0);

 }
 return;
};