#include <iostream>
#include "mesh.h"
#include <cmath>
using namespace std;

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

void FVcell :: setCellLabel(int i, int j, int N_x, int N_y)
{
  // setting the ID for a given cell according to its indexing (i,j)

  // default: internal node
  this->ID.INTERNAL = 1;
  this->ID.EAST     = 0;
  this->ID.WEST     = 0;
  this->ID.NORD     = 0;
  this->ID.SUD      = 0;

  if (i <= 1) { // west node
    this->ID.INTERNAL = 0;
    this->ID.WEST     = 1;
  }

  if (i >= N_x) { // east node
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
volumeField volumeField::operator+(const volumeField& rhs)
 {
   volumeField tmp = *this; 

   for (int i= 0; i<=this->dimx+1; i++){
   for (int j= 0; j <=this->dimy+1; j++){
     tmp.mesh[i][j].setVal(tmp.mesh[i][j].cellVal() + 
                             rhs.mesh[i][j].cellVal() );
    };
   };
   return tmp;
 };

volumeField volumeField::operator-(const volumeField& rhs)
 {

   for (int i= 0; i<=this->dimx+1; i++){
    for (int j= 0; j <=this->dimy+1; j++){
     this->mesh[i][j].setVal(mesh[i][j].cellVal() - 
                            rhs.mesh[i][j].cellVal() );
    };
   };
   return *this;
 };


volumeField volumeField::operator*(double scalar_val)
 {

   volumeField tmp = *this; 

   for (int i= 0; i<=this->dimx+1; i++){
   for (int j= 0; j <=this->dimy+1; j++){
     tmp.mesh[i][j].setVal(tmp.mesh[i][j].cellVal() * scalar_val );
    };
   };
   return tmp;
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
     // yj=1−tanh[γ(1–2jN2)]tanh(γ)               //(j=0,…,N2),(3)
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

 double xf[N_x+2];
 double xc[N_x+2];
 double yf[N_y+2];
 double yc[N_y+2];

 volumeField field;
 double gamma;
 gamma  = 1. ;
 FVcell **mesh_p = new FVcell* [N_x+2];

  for (int i = 0; i <= N_x; i++)
  {  
   xf[i] = L_x * 0.5*(1. - tanh(gamma*(1.-2.*i/(N_x)))/tanh(gamma));
   //xf[i] = L_x *i/(N_x);
  }
  xf[N_x+1] = L_x + (L_x - xf[N_x-1]);


  for (int i = 1; i <= N_x; i++)
  {  
    xc[i] =  xf[i-1] + 0.5*(xf[i]-xf[i-1]);
  }

  
  xc[0]   =  0 - xc[1];
  xc[N_x+1]=  L_x + (L_x-xc[N_x]);

  for (int j = 0; j <= N_y; j++)
  {
    yf[j] = L_y * 0.5*(1. - tanh(gamma*(1.-2.*j/(N_y)))/tanh(gamma));
    //yf[j] = L_y * j/(N_y);
  }

  yf[N_x+1] = L_y + (L_y - yf[N_y-1]);
  

  for (int j = 1; j <= N_y; j++)
  {  
    yc[j] =  yf[j-1] + 0.5*(yf[j]-yf[j-1]);
  }

  yc[0]   =  0 - yc[1];
  yc[N_y+1]=  L_y + (L_y-yc[N_y]);

  for (int i = 0; i <= N_x+1; i++)
  {       
   mesh_p[i] = new FVcell [N_y+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
    double dx = 2.*(xf[i]  - xc[i]);   
    double dy = 2.*(yf[j]  - yc[j]);     
    mesh_p[i][j].setCellPos( xc[i], yc[j], dx, dy  );
    mesh_p[i][j].setCellLabel( i, j, N_x, N_y);

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
void updateFieldVal (volumeField* field_to_update, volumeField field) //FVcell** mesh_p, int N_x, int N_y, double val_in)
{
 int N_x = field.dimx;
 int N_y = field.dimy;

 for (int i = 0; i <= N_x+1; i++)
 {
   for(int j = 0; j<= N_y+1; j++)
   {
    field_to_update->mesh[i][j].setVal(field.mesh[i][j].cellVal());
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
    field->mesh[N_x+1][j].setVal( -field->mesh[N_x-1][j].cellVal()  );
 }

// along x-axis     
 for (int i = 0; i <= N_x+1; i++)
 {
    field->mesh[i][  0  ].setVal(     - field->mesh[i][ 1 ].cellVal() )  ;
    field->mesh[i][N_y+1].setVal( 2 - field->mesh[i][N_y].cellVal() ) ;
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCvFV (volumeField* field)
{
 // set wall boundary condition for staggered Uy velocity field  
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
    field->mesh[i][N_y+1].setVal(- field->mesh[i][N_y-1].cellVal() );
    field->mesh[i][0].setVal(0 );
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

 double val = 0;

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 { 
    p->mesh[  0  ][j].setVal(  p->mesh[ 1 ][j].cellVal()  );
    p->mesh[N_x+1][j].setVal(    p->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    p->mesh[i][  0  ].setVal(    p->mesh[i][ 1   ].cellVal() );
    p->mesh[i][N_y+1].setVal(    p->mesh[i][ N_y ].cellVal() );
     }
 return;
};
//////////////////////////////////////////////////////////////////////
void divergence (volumeField* div, const volumeField  Ux_pred, const volumeField  Uy_pred)
{
 int N_x = Ux_pred.dimx;
 int N_y = Ux_pred.dimy;

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
        div->mesh[i][j].setVal( ( Ux_pred.mesh[i][ j ].cellVal() - Ux_pred.mesh[i-1][j].cellVal() ) / Ux_pred.mesh[i][j].dx() +
                               ( Uy_pred.mesh[ i ][j].cellVal() - Uy_pred.mesh[i][j-1].cellVal() ) / Uy_pred.mesh[i][j].dy() ); 
  };
};
 return;
};
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void setBCneumannFV (volumeField* p)
{
 // set wall boundary condition for pressure  
 int N_x, N_y;

 N_x = p->dimx;
 N_y = p->dimy;

 double val = 0;

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 { 
    p->mesh[  0  ][j].setVal(  p->mesh[ 1 ][j].cellVal()  );
    p->mesh[N_x+1][j].setVal(   p->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    p->mesh[i][  0  ].setVal( p->mesh[i][ 1   ].cellVal() );
    p->mesh[i][N_y+1].setVal( p->mesh[i][ N_y ].cellVal() );
     }
 return;
};
void setBCerrFV (volumeField* p)
{
 // set wall boundary condition for pressure  
 int N_x, N_y;

 N_x = p->dimx;
 N_y = p->dimy;

 double val = 0;

 // along y-axis        
 for (int j = 0; j <= N_y+1; j++)
 { 
    p->mesh[  0  ][j].setVal( - p->mesh[ 1 ][j].cellVal()  );
    p->mesh[N_x+1][j].setVal( -  p->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    p->mesh[i][  0  ].setVal( -p->mesh[i][ 1   ].cellVal() );
    p->mesh[i][N_y+1].setVal(- p->mesh[i][ N_y ].cellVal() );
     }
 return;
};
//////////////////////////////////////////////////////////////////////
void correct(volumeField* Ux_corr, volumeField* Uy_corr, const volumeField  Ux_pred, const volumeField  Uy_pred, const volumeField p, const double rho, const double dt)
{
 int N_x = Ux_pred.dimx;
 int N_y = Ux_pred.dimy;

 double gradp_x;
 double gradp_y;

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {

        gradp_x = (p.mesh[i+1][j].cellVal() - p.mesh[i][j].cellVal())/(p.mesh[i+1][j].xCentroid()-p.mesh[i][j].xCentroid()   );
        gradp_y = (p.mesh[i][j+1].cellVal() - p.mesh[i][j].cellVal())/(p.mesh[i][j+1].yCentroid()-p.mesh[i][j].yCentroid()   );

        Ux_corr->mesh[i][j].setVal( Ux_pred.mesh[i][j].cellVal() - (dt/rho) * gradp_x );  
        Uy_corr->mesh[i][j].setVal( Uy_pred.mesh[i][j].cellVal() - (dt/rho) * gradp_y );  
  };
};
 return;
};
//////////////////////////////////////////////////////////////////////////////
void interpolateFieldVal(const volumeField& field_fine, volumeField* field_coarse )
{
  int Nx_c = field_coarse->dimx;
  int Ny_c = field_coarse->dimy;

  int Nx_f = field_fine.dimx;
  int Ny_f = field_fine.dimy;

   double local_value;

  for (int i_c = 1; i_c<=Nx_c; i_c++) {
   for (int j_c = 1; j_c<=Ny_c; j_c++) {

        int i_f = i_c*2 - 1;
        int j_f = j_c*2 - 1;

        local_value = 0.25* ( field_fine.mesh [i_f ][ j_f ].cellVal() +
                              field_fine.mesh[ i_f ][j_f+1].cellVal() + 
                              field_fine.mesh[i_f+1][ j_f ].cellVal() +
                              field_fine.mesh[i_f+1][j_f+1].cellVal() );

        field_coarse->mesh[i_c][j_c].setVal(local_value); 
  };
}; 
return;
}
//////////////////////////////////////////////////////////////////////////////
void prolongateFieldVal(const volumeField& field_coarse, volumeField* field_fine )
{
  int Nx_c = field_coarse.dimx;
  int Ny_c = field_coarse.dimy;

  int Nx_f = field_fine->dimx;
  int Ny_f = field_fine->dimy;

   double local_value;

  for (int i_c = 1; i_c<=Nx_c; i_c++) {
   for (int j_c = 1; j_c<=Ny_c; j_c++) {

        int i_f = i_c*2 - 1;
        int j_f = j_c*2 - 1;

        local_value = field_coarse.mesh[i_c][j_c].cellVal(); 
        
        field_fine->mesh [i_f ][ j_f ].setVal(local_value); 
        field_fine->mesh[ i_f ][j_f+1].setVal(local_value);  
        field_fine->mesh[i_f+1][ j_f ].setVal(local_value); 
        field_fine->mesh[i_f+1][j_f+1].setVal(local_value);

  };
}; 
return;
}

///////////////////////////////////////////////////////////////////////////////////////////
double l1_norm(volumeField const field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;

    double accum = 0.;

    for (int i = 1; i <= N_x; ++i)
    {
         for (int j= 1; j <= N_y; ++j)
         {
          accum += fabs( field.mesh[i][j].cellVal() ) ;
         };
    };
    return sqrt(accum);
};
///////////////////////////////////////////////////////////////////
double sum (const volumeField& field)
{

    int N_x = field.dimx;
    int N_y = field.dimy;

    double accum = 0.;

    for (int i = 1; i <= N_x; i++)
    {
         for (int j= 1; j <= N_y; j++)
         {
          accum = accum + ( field.mesh[i][j].cellVal() ) ;
         };
    };
    return (accum);
};
