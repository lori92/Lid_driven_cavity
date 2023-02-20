#include <iostream>
#include "mesh.h"
using namespace std;


////////////////////////////////////////////////////////////////
void initFieldVal (volumeField field, double val_in) 
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
    field->mesh[i][N_y+1].setVal( 1. - field->mesh[i][N_y].cellVal() ) ;
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
    p->mesh[N_x+1][j].setVal(   p->mesh[N_x][j].cellVal()  );
 }

 // along x-axis        
 for (int i = 0; i <= N_x+1; i++)
 {
    p->mesh[i][  0  ].setVal(   p->mesh[i][ 1   ].cellVal() );
    p->mesh[i][N_y+1].setVal( p->mesh[i][ N_y ].cellVal() );
     }
 return;
};