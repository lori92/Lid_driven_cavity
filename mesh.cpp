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
        cout << value;};
////////////////////////////////////////////////////////////////
cell** buildMesh (int N_x, int N_y, double L_x, double L_y)
{

 cell **mesh_p = new cell* [N_x+2];
 for (int i = 0; i <= N_x+1; i++)
 {       
   mesh_p[i] = new cell [N_y+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
    mesh_p[i][j].setCellPos( L_x/N_x/2. +  (L_x/N_x)*(i-1), L_y/N_y/2. + (L_y/N_y)*(j-1), L_x/N_x, L_y/N_y  );
   }     
 }       

 return mesh_p;
};

////////////////////////////////////////////////////////////////
void initMeshVal (cell** mesh_p, int N_x, int N_y, double val_in)
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


///////////////////////////////////////////////////////////////
// generic class for finite volume field //
class volumeField
{
 private:
  cell** mesh;
   
 public:
  
};
