#include <iostream>
#include "mesh.h"
using namespace std;

void FVcell :: setCellPos(double x, double y, double z, double dx, double dy, double dz)
{
	// x, z coordinates of the cell center of a given cell corresponding to index i, j
        xC = x; 
        xFaces[0] = xC - dx/2;
        xFaces[1] = xC + dx/2;
 
        yC = y; 
	yFaces[0] = yC - dy/2;
	yFaces[1] = yC + dy/2;

          zC = z; 
	zFaces[0] = zC - dz/2;
	zFaces[1] = zC + dz/2;
        	
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

double FVcell :: dz(){
       return this->zFaces[1]-this->zFaces[0];};

double FVcell :: xFace(int face ){
       return this->xFaces[face];};

double FVcell :: yFace(int face){
       return this->yFaces[face];};

double FVcell :: zFace(int face){
       return this->zFaces[face];};

double FVcell :: xCentroid(){
       return this->xC;};

double FVcell :: yCentroid(){
       return this->yC;};

double FVcell :: zCentroid(){
       return this->zC;};
///////////////////////////////////////////////////////////////
volumeField volumeField::operator+(const volumeField& rhs)
 {
   volumeField tmp = *this; 

   for (int i= 0; i<=this->dimx+1; i++){
     for (int j= 0; j <=this->dimy+1; j++){
       for (int k= 0; k <=this->dimz+1; k++){

     tmp.mesh[i][j][k].setVal(tmp.mesh[i][j][k].cellVal() + 
                             rhs.mesh[i][j][k].cellVal() );
     };
    };
   };
   return tmp;
 };

volumeField volumeField::operator-(const volumeField& rhs)
 {

   for (int i= 0; i<=this->dimx+1; i++){
    for (int j= 0; j <=this->dimy+1; j++){
     for (int k= 0; k <=this->dimz+1; k++){

     this->mesh[i][j][k].setVal(mesh[i][j][k].cellVal() - 
                            rhs.mesh[i][j][k].cellVal() );
     };
    };
   };
   return *this;
 };


volumeField volumeField::operator*(double scalar_val)
 {
   volumeField tmp = *this; 

   for (int i= 0; i<=this->dimx+1; i++){
    for (int j= 0; j <=this->dimy+1; j++){
     for (int k= 0; k <=this->dimz+1; k++){

     tmp.mesh[i][j][k].setVal(tmp.mesh[i][j][k].cellVal() * scalar_val );
     };
    };
   };
   return tmp;
 };

////////////////////////////////////////////////////////////////
FVcell** buildMesh (int N_x, int N_y, int N_z, double L_x, double L_y, double L_z)
{
 // index i spans the y-direction
 // index j spans the x-direction
 // index k spans the z-direction

     // Allocate memory blocks of size
    // x i.e., no of 2D Arrays
  FVcell*** mesh_p = new FVcell**[N_x+2];
 
  for (int i = 0; i <= N_x+1; i++) {
        // Allocate memory blocks for
        // rows of each 2D array
        mesh_p[i] = new FVcell*[N_y+2];

        for (int j = 0; j <= N_y+1; j++) {
 
            // Allocate memory blocks for
            // columns of each 2D array
            mesh_p[i][j] = new int[N_z+2];
        }
    }

    for (int i = 0; i < N_x; i++) {
        for (int j = 0; j < N_y; j++) {
            for (int k = 0; k < N_z; k++) {
 
                // Assign values to the
                // memory blocks created
                 mesh_p[i][j][k].setCellPos( L_x/N_x/2. + (L_x/N_x)*(i-1), 
                                             L_y/N_y/2. + (L_y/N_y)*(j-1), 
                                             L_z/N_z/2. + (L_z/N_z)*(k-1),
                                             L_y/N_y, L_x/N_x, L_z/N_z  );
            }
        }
    }


/*
 FVcell ***mesh_p = new FVcell** [N_x+2];
 for (int i = 0; i <= N_x+1; i++)
 {       
   *mesh_p[i] = new FVcell [N_y+2][N_z+2];
   for(int j = 0; j<= N_y+1; j++)
   {       
     // yj=1−tanh[γ(1–2jN2)]tanh(γ)               //(j=0,…,N2),(3)
    mesh_p[i][j].setCellPos( L_y/N_y/2. +  (L_y/N_y)*(i-1), L_x/N_x/2. + (L_x/N_x)*(j-1), L_y/N_y, L_x/N_x  );
   }     
 }       
*/
 return mesh_p;
};
////////////////////////////////////////////////////////////////
volumeField buildFVfield (int N_x, int N_y, int N_z, double L_x, double L_y, double Lz)
{
 // index i spans the y-direction
 // index j spans the x-direction

 double xf[N_x+1];
 double xc[N_x+2];
 double yf[N_y+1];
 double yc[N_y+2];
 double zf[N_z+1];
 double zc[N_z+2];

 volumeField field;
 double gamma;
 gamma  = 1.;
 FVcell **mesh_p = new FVcell* [N_x+2];

  for (int i = 0; i <= N_x+1; i++)
  {  
   //xf[i] = L_x * 0.5*(1. - tanh(gamma*(1.-2.*i/(N_x)))/tanh(gamma));
   xf[i] = L_x *i/(N_x);

  }
  for (int i = 1; i <= N_x; i++)
  {  
    xc[i] =  xf[i-1] + 0.5*(xf[i]-xf[i-1]);
  }

  xc[0]   =  0 - xc[1];
  xc[N_x+1]=  L_x + (L_x-xc[N_x]);

  for (int j = 0; j <= N_y+1; j++)
    {
      //yf[j] = L_y * 0.5*(1. - tanh(gamma*(1.-2.*j/(N_y)))/tanh(gamma));
      yf[j] = L_y * j/(N_y);
    }

  for (int j = 1; j <= N_y+1; j++)
  {  
    yc[j] =  yf[j-1] + 0.5*(yf[j]-yf[j-1]);
  }

  yc[0]   =  0 - yc[1];
  yc[N_y+1]=  L_y + (L_y-yc[N_y]);


  for (int k = 0; k <= N_z+1; k++)
  {  
   //zf[i] = L_z * 0.5*(1. - tanh(gamma*(1.-2.*k/(N_z)))/tanh(gamma));
   zf[k] = L_z *k/(N_z);
  }

  for (int k = 1; k <= N_z+1; k++)
  {  
    zc[k] =  zf[k-1] + 0.5*(zf[k]-zf[k-1]);
  }

  zc[0]   =  0 - zc[1];
  zc[N_z+1]=  L_z + (L_z-zc[N_z]);

  FVcell*** mesh_p = new FVcell**[N_x+2];
  for (int i = 0; i <= N_x+1; i++) {
        // Allocate memory blocks for
        // rows of each 2D array
        mesh_p[i] = new FVcell*[N_y+2];

        for (int j = 0; j <= N_y+1; j++) {
 
            // Allocate memory blocks for
            // columns of each 2D array
            mesh_p[i][j] = new int[N_z+2];
        }
    }

    for (int i = 0; i < N_x; i++) {
        for (int j = 0; j < N_y; j++) {
            for (int k = 0; k < N_z; k++) {
 
                // Assign values to the
                // memory blocks created
                 mesh_p[i][j][k].setCellPos( L_x/N_x/2. + (L_x/N_x)*(i-1), 
                                             L_y/N_y/2. + (L_y/N_y)*(j-1), 
                                             L_z/N_z/2. + (L_z/N_z)*(k-1),
                                             L_y/N_y, L_x/N_x, L_z/N_z  );
            }
        }
    }

 field.mesh = mesh_p;
 field.dimx = N_x;
 field.dimy = N_y;
 field.dimz = N_z;

 return field;
};
////////////////////////////////////////////////////////////////
void initFieldVal (volumeField field, double val_in) //FVcell** mesh_p, int N_x, int N_y, double val_in)
{

 int N_x = field.dimx;
 int N_y = field.dimy;
 int N_z = field.dimz;

 for (int i = 0; i <= N_x+1; i++)
 {
   for(int j = 0; j<= N_y+1; j++)
   {
    for(int k = 0; k<= N_z+1; k++)
    {
     field.mesh[i][j][k].setVal(val_in);
    }
   } 
 }
 return;
};
////////////////////////////////////////////////////////////////
void updateFieldVal (volumeField* field_to_update, volumeField field) //FVcell** mesh_p, int N_x, int N_y, double val_in)
{
 int N_x = field.dimx;
 int N_y = field.dimy;
 int N_z = field.dimz;

 for (int i = 0; i <= N_x+1; i++)
 {
   for(int j = 0; j<= N_y+1; j++)
   {
    for(int k = 0; k<= N_z+1; k++)
    {
     field_to_update->mesh[i][j][k].setVal(field.mesh[i][j][k].cellVal());
    }
   }
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCuFV (volumeField* field)
{
 // set wall boundary condition for staggered Ux velocity field  
 int N_x, N_y, N_z;

 N_x = field->dimx;
 N_y = field->dimy;
 N_z = field->dimz;

// along x-planes     
 for (int j= 0; j <= N_y+1; j++)
 {
  for (int k= 0; k <= N_z+1; k++)
  {
    field->mesh[ 0   ][j][k].setVal(0);
    field->mesh[N_x  ][j][k].setVal(0);
    field->mesh[N_x+1][j][k].setVal( -field->mesh[N_x-1][j][k].cellVal()  );
  }
 }

// along y-planes     
 for (int i = 0; i <= N_x+1; i++)
 {
  for (int k = 0; k <= N_z+1; k++)
  {
    field->mesh[i][  0  ][k].setVal(     - field->mesh[i][ 1 ][k].cellVal() )  ;
    field->mesh[i][N_y+1][k].setVal( 1.  - field->mesh[i][N_y][k].cellVal() ) ;
  }
 }

// along z-planes     
 for (int i = 0; i <= N_x+1; i++)
 {
  for (int j = 0; j <= N_y+1; j++)
  {
    field->mesh[i][  j  ][0].setVal(     - field->mesh[i][ j ][1].cellVal() )  ;
    field->mesh[i][j][N_z+1].setVal(     - field->mesh[i][j][N_z].cellVal() ) ;
  }
 }

 return;
};
////////////////////////////////////////////////////////////////
void setBCvFV (volumeField* field)
{
 // set wall boundary condition for staggered Uy velocity field  
 // set wall boundary condition for staggered Ux velocity field  
 int N_x, N_y, N_z;

 N_x = field->dimx;
 N_y = field->dimy;
 N_z = field->dimz;

 // along x-planes       
 for (int j = 0; j <= N_y+1; j++)
 {
   for (int k = 0; k <= N_z+1; k++)
    {
     field->mesh[  0  ][j][k].setVal( - field->mesh[ 1 ][j][k].cellVal()  );
     field->mesh[N_x+1][j][k].setVal( - field->mesh[N_x][j][k].cellVal()  );
    }
 }

 // along y-planes        
 for (int i = 0; i <= N_x+1; i++)
 {
   for (int k = 0; k <= N_z+1; k++)
    {
     field->mesh[i][ N_y ][k].setVal(0);
     field->mesh[i][N_y+1][k].setVal(- field->mesh[i][N_y-1][k].cellVal() );
     field->mesh[i][0][k].setVal(0 );
    }
 }

// along z-planes     
 for (int i = 0; i <= N_x+1; i++)
 {
  for (int j = 0; j <= N_y+1; j++)
  {
    field->mesh[i][  j  ][0].setVal(     - field->mesh[i][ j ][1].cellVal() )  ;
    field->mesh[i][j][N_z+1].setVal(     - field->mesh[i][j][N_z].cellVal() ) ;
  }
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCwFV (volumeField* field)
{
 // set wall boundary condition for staggered Uy velocity field  
 // set wall boundary condition for staggered Ux velocity field  
 int N_x, N_y, N_z;

 N_x = field->dimx;
 N_y = field->dimy;
 N_z = field->dimz;

 // along x-planes       
 for (int j = 0; j <= N_y+1; j++)
 {
   for (int k = 0; k <= N_z+1; k++)
    {
     field->mesh[  0  ][j][k].setVal( - field->mesh[ 1 ][j][k].cellVal()  );
     field->mesh[N_x+1][j][k].setVal( - field->mesh[N_x][j][k].cellVal()  );
    }
 }

 // along y-planes        
 for (int i = 0; i <= N_x+1; i++)
 {
   for (int k = 0; k <= N_z+1; k++)
    {
     field->mesh[i][N_y+1][k].setVal(- field->mesh[i][N_y-1][k].cellVal() );
     field->mesh[i][0][k].setVal(0 );
    }
 }

// along z-planes     
 for (int i = 0; i <= N_x+1; i++)
 {
  for (int j = 0; j <= N_y+1; j++)
  {
    field->mesh[i][ j ][N_z].setVal(0);
    field->mesh[i][  j  ][0].setVal(0);
    field->mesh[i][j][N_z+1].setVal(     - field->mesh[i][j][N_z-1].cellVal() ) ;
  }
 }
 return;
};
////////////////////////////////////////////////////////////////
void setBCpFV (volumeField* p)
{
 // set wall boundary condition for pressure  
 int N_x, N_y, N_z;

 N_x = p->dimx;
 N_y = p->dimy;
 N_z = p->dimz;

 double val = 0;

 // along x-plane        
 for (int j = 0; j <= N_y+1; j++)
 {
   for (int k = 0; k <= N_z+1; k++)
   { 
    p->mesh[  0  ][j][k].setVal(  p->mesh[ 1 ][j][k].cellVal()  );
    p->mesh[N_x+1][j][k].setVal(  p->mesh[N_x][j][k].cellVal()  );
   }
 }

 // along y-plane        
 for (int i = 0; i <= N_x+1; i++)
 {
    for (int k = 0; k <= N_z+1; k++)
    { 
     p->mesh[i][  0  ][k].setVal(   p->mesh[i][ 1   ][k].cellVal() );
     p->mesh[i][N_y+1][k].setVal( p->mesh[i][ N_y ][k].cellVal() );
    }
 }

 // along z-plane        
 for (int i = 0; i <= N_x+1; i++)
 {
   for (int j = 0; j <= N_y+1; j++)
    { 
     p->mesh[i][  j  ][0].setVal( p->mesh[i][ j   ][1].cellVal() );
     p->mesh[i][j][N_z+1].setVal( p->mesh[i][ j ][N_z].cellVal() );
    }
 }

 return;
};
//////////////////////////////////////////////////////////////////////
void divergence (volumeField* div, const volumeField  Ux_pred, const volumeField  Uy_pred, const volumeField  Uz_pred)
{
 int N_x = Ux_pred.dimx;
 int N_y = Ux_pred.dimy;
 int N_z = Ux_pred.dimz;

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
   for (int k = 1; k<=N_z; k++) {
        div->mesh[i][j].setVal( ( Ux_pred.mesh[i][ j ][k].cellVal() - Ux_pred.mesh[i-1][j][k].cellVal() ) / Ux_pred.mesh[i][j][k].dx() +
                                ( Uy_pred.mesh[ i ][j][k].cellVal() - Uy_pred.mesh[i][j-1][k].cellVal() ) / Uy_pred.mesh[i][j][k].dy() +
                                ( Uz_pred.mesh[ i ][j][k].cellVal() - Uz_pred.mesh[i][j-1][k].cellVal() ) / Uz_pred.mesh[i][j][k].dz() ); 
    };
  };
};
 return;
};
////////////////////////////////////////////////////////////////
void setBCneumannFV (volumeField* p)
{
 // set wall boundary condition for pressure  
 int N_x, N_y, N_z;

 N_x = p->dimx;
 N_y = p->dimy;
 N_z = p->dimz;

 double val = 0;

 // along x-planes        
 for (int k = 1; i<=N_z; k++) {
  for (int j = 1; j<=N_y; j++)  {

    p->mesh[  0  ][j][k].setVal(  p->mesh[ 1 ][j][k].cellVal()  );
    p->mesh[N_x+1][j][k].setVal(  p->mesh[N_x][j][k].cellVal()  );
 }
 }
 // along y-planes        
 for (int k = 1; i<=N_z; k++) {
  for (int i = 1; i<=N_x; i++) {
    p->mesh[i][  0  ][k].setVal( p->mesh[i][ 1   ][k].cellVal() );
    p->mesh[i][N_y+1][k].setVal( p->mesh[i][ N_y ][k].cellVal() );
  }
 }

 // along z-planes        
 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
    p->mesh[i][j][   0 ].setVal( p->mesh[i][ j ][ 1 ].cellVal() );
    p->mesh[i][j][N_z+1].setVal( p->mesh[i][ j ][N_z].cellVal() );
  }
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

 // along x-planes        
 for (int k = 1; i<=N_z; k++) {
  for (int j = 1; j<=N_y; j++)  {

    p->mesh[  0  ][j][k].setVal(-p->mesh[ 1 ][j][k].cellVal()  );
    p->mesh[N_x+1][j][k].setVal(-p->mesh[N_x][j][k].cellVal()  );
  }
 }
 // along y-planes        
 for (int k = 1; i<=N_z; k++) {
  for (int i = 1; i<=N_x; i++) {
    p->mesh[i][  0  ][k].setVal(-p->mesh[i][ 1   ][k].cellVal() );
    p->mesh[i][N_y+1][k].setVal(-p->mesh[i][ N_y ][k].cellVal() );
  }
 }

 // along z-planes        
 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
    p->mesh[i][j][   0 ].setVal(-p->mesh[i][ j ][ 1 ].cellVal() );
    p->mesh[i][j][N_z+1].setVal(-p->mesh[i][ j ][N_z].cellVal() );
  }
 }
 return;
};
//////////////////////////////////////////////////////////////////////
void correct(     volumeField* Ux_corr,       volumeField* Uy_corr,       volumeField* Uz_corr,
            const volumeField  Ux_pred, const volumeField  Uy_pred, const volumeField  Uz_pred,
            const volumeField p, const double rho, const double dt)
{
 int N_x = Ux_pred.dimx;
 int N_y = Ux_pred.dimy;
 int N_z = Uz_pred.dimy;

 double gradp_x;
 double gradp_y;
 double gradp_z;

 for (int i = 1; i<=N_x; i++) {
  for (int j = 1; j<=N_y; j++) {
   for (int k = 1; k<=N_z; k++) {

        gradp_x = (p.mesh[i+1][j][k].cellVal() - p.mesh[i][j][k].cellVal())/(p.mesh[i+1][j][k].xCentroid()-p.mesh[i][j][k].xCentroid()   );
        gradp_y = (p.mesh[i][j+1][k].cellVal() - p.mesh[i][j][k].cellVal())/(p.mesh[i][j+1][k].yCentroid()-p.mesh[i][j][k].yCentroid()   );
        gradp_z = (p.mesh[i][j][k+1].cellVal() - p.mesh[i][j][k].cellVal())/(p.mesh[i][j][k+1].zCentroid()-p.mesh[i][j][k].zCentroid()   );

        Ux_corr->mesh[i][j][k].setVal( Ux_pred.mesh[i][j][k].cellVal() - (dt/rho) * gradp_x );  
        Uy_corr->mesh[i][j][k].setVal( Uy_pred.mesh[i][j][k].cellVal() - (dt/rho) * gradp_y );  
        Uz_corr->mesh[i][j][k].setVal( Uz_pred.mesh[i][j][k].cellVal() - (dt/rho) * gradp_z );  

    };
  };
};
 return;
};
//////////////////////////////////////////////////////////////////////////////
void interpolateFieldVal(const volumeField& field_fine, volumeField* field_coarse )
{
  int Nx_c = field_coarse->dimx;
  int Ny_c = field_coarse->dimy;
  int Nz_c = field_coarse->dimz;

  int Nx_f = field_fine.dimx;
  int Ny_f = field_fine.dimy;
  int Nz_f = field_fine.dimz;

   double local_value;

  for (int i_c = 1; i_c<=Nx_c; i_c++) {
   for (int j_c = 1; j_c<=Ny_c; j_c++) {
     for (int k_c = 1; k_c<=Nz_c; k_c++) {


        int i_f = i_c*2 - 1;
        int j_f = j_c*2 - 1;
        int k_f = k_c*2 - 1;


        local_value = 0.125* ( field_fine.mesh[i_f  ][ j_f ][ k_f ].cellVal() +
                               field_fine.mesh[ i_f ][j_f+1][ k_f ].cellVal() + 
                               field_fine.mesh[i_f+1][ j_f ][ k_f ].cellVal() +
                               field_fine.mesh[i_f+1][j_f+1][ k_f ].cellVal() +  
                               field_fine.mesh[ i_f ][ j_f ][k_f+1].cellVal() +
                               field_fine.mesh[ i_f ][ j_f ][k_f+1].cellVal() );

        field_coarse->mesh[i_c][j_c][k_c].setVal(local_value); 
  };
}; 
return;
}
//////////////////////////////////////////////////////////////////////////////
void prolongateFieldVal(const volumeField& field_coarse, volumeField* field_fine )
{
  int Nx_c = field_coarse.dimx;
  int Ny_c = field_coarse.dimy;
  int Nz_c = field_coarse.dimz;

  int Nx_f = field_fine->dimx;
  int Ny_f = field_fine->dimy;
  int Nz_f = field_fine->dimz;

   double local_value;

  for (int i_c = 1; i_c<=Nx_c; i_c++) {
   for (int j_c = 1; j_c<=Ny_c; j_c++) {
    for (int k_c = 1; k_c<=Nz_c; k_c++) {

        int i_f = i_c*2 - 1;
        int j_f = j_c*2 - 1;
        int k_f = k_c*2 - 1;

        local_value = field_coarse.mesh[i_c][j_c][k_c].cellVal(); 
        
        field_fine->mesh [i_f ][ j_f ] [k_f ].setVal(local_value); 
        field_fine->mesh[ i_f ][j_f+1][ k_f ].setVal(local_value);  
        field_fine->mesh[i_f+1][ j_f ][ k_f ].setVal(local_value); 
        field_fine->mesh[i_f+1][j_f+1][ k_f ].setVal(local_value);

        field_fine->mesh [i_f ][ j_f ][k_f+1].setVal(local_value); 
        field_fine->mesh[ i_f ][j_f+1][k_f+1].setVal(local_value);  
        field_fine->mesh[i_f+1][ j_f ][k_f+1].setVal(local_value); 
        field_fine->mesh[i_f+1][j_f+1][k_f+1].setVal(local_value);       
    };
  };
}; 
return;
}