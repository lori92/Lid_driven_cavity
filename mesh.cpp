#include <iostream>
using namespace std;



class FVcell
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
 double xCentroid ();
 double yCentroid ();
 double xFace (int);
 double yFace (int) ;

};

void FVcell :: setCellPos(double x, double y, double dx, double dy)
{
	// x, z coordinates of the cell center of a given cell corresponding to index i, j
        xC = x; 
        xFaces[0] = xC - dx/2;
        xFaces[1] = xC + dx/2;
 
        yC = y; 
	yFaces[0] = yC - dy/2;
	yFaces[1] = yC + dy/2;
        	
};      

void FVcell :: setVal(double value_in)
{
        //  value of the volumeField for given cell corresponding to index i, j
	value = value_in;
};


void FVcell :: print_cellCentre(){
        cout << "(" <<xC<<","<<yC<<") ";};

void FVcell :: print_cellVal(){
        cout << value<<" ";};

double FVcell :: cellVal(){
       return this->value;};
       
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
struct volumeField
{
 FVcell** mesh;
 int dimx;
 int dimy;
};
////////////////////////////////////////////////////////////////
FVcell** buildMesh (int N_x, int N_y, double L_x, double L_y)
{
 // index i spans the y-direction
 // index j spans the x-direction

 FVcell **mesh_p = new FVcell* [N_x+2];
 for (int i = 0; i <= N_x+1; i++)
 {       
   mesh_p[i] = new FVcell [N_y+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
    mesh_p[i][j].setCellPos( L_y/N_y/2. +  (L_y/N_y)*(i-1), L_x/N_x/2. + (L_x/N_x)*(j-1), L_y/N_y, L_x/N_x  );
   }     
 }       

 return mesh_p;
};
////////////////////////////////////////////////////////////////
volumeField buildFVfield (int N_x, int N_y, double L_x, double L_y)
{
 // index i spans the y-direction
 // index j spans the x-direction
 volumeField field;
 FVcell **mesh_p = new FVcell* [N_x+2];
 for (int i = 0; i <= N_x+1; i++)
 {       
   mesh_p[i] = new FVcell [N_y+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
    mesh_p[i][j].setCellPos( L_x/N_x/2. +  (L_x/N_x)*(i-1), L_y/N_y/2. + (L_y/N_y)*(j-1), L_x/N_x, L_y/N_y  );
   }     
 }       
 field.mesh = mesh_p;
 field.dimx = N_x;
 field.dimy = N_y;
 return field;
};

////////////////////////////////////////////////////////////////
void initFieldVal (volumeField field, double val_in) //FVcell** mesh_p, int N_x, int N_y, double val_in)
{

 int N_x = field.dimx;
 int N_y = field.dimy;

 for (int i = 0; i <= N_x+1; i++)
 {
   for(int j = 0; j<= N_y+1; j++)
   {
    field.mesh[i][j].setVal(val_in);
   }
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCuFV (volumeField* field)
{
 // set wall boundary condition for staggered Ux velocity field  
 int N_x, N_y;

 N_x = field->dimx;
 N_y = field->dimy;

// along y-axis     
 for (int j= 0; j <= N_y+1; j++)
 {
    field->mesh[ 0   ][j].setVal(0);
    field->mesh[N_x  ][j].setVal(0);
    field->mesh[N_x+1][j].setVal(0);
 }

// along x-axis     
 for (int i = 0; i <= N_x+1; i++)
 {
    field->mesh[i][  0  ].setVal( - field->mesh[i][ 1 ].cellVal() )  ;
    field->mesh[i][N_x+1].setVal( 100 - field->mesh[i][N_x].cellVal() )  ;
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCvFV (volumeField* field)
{
 // set wall boundary condition for staggered Uy velocity field  
 // set wall boundary condition for staggered Ux velocity field  
 int N_x, N_y;

 N_x = field->dimx;
 N_y = field->dimy;

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 {
    field->mesh[  0  ][j].setVal( - field->mesh[ 1 ][j].cellVal()  );
    field->mesh[N_x+1][j].setVal( - field->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    field->mesh[i][ N_y ].setVal(0);
    field->mesh[i][N_y+1].setVal( -field->mesh[i][ N_y ].cellVal() );
    field->mesh[i][  0  ].setVal(0);
     }
 return;
};
////////////////////////////////////////////////////////////////
void setBCpFV (volumeField* p)
{
 // set wall boundary condition for pressure  
 int N_x, N_y;

 N_x = p->dimx;
 N_y = p->dimy;

 double val = 10;

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 { 
    p->mesh[  0  ][j].setVal(0.1*    val - p->mesh[ 1 ][j].cellVal()  );
    p->mesh[N_x+1][j].setVal( 2.*val - p->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    p->mesh[i][  0  ].setVal( 2.*val - p->mesh[i][ 1   ].cellVal() );
    p->mesh[i][N_y+1].setVal( 2.*val - p->mesh[i][ N_y ].cellVal() );
     }
 return;
};
