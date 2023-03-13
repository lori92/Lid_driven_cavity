#ifndef FVCELL_H
#define FVCELL_H

struct ID_CELL
{ 
   bool INTERNAL;
   bool NORD;
   bool SUD;
   bool EAST;
   bool WEST;
};

class FVcell
{       
 private:
  double xC;
  double yC;
  double xFaces[2];
  double yFaces[2];
  double value;
  ID_CELL ID;

 public:

 void setCellPos   (double, double, double, double);
 void setCellLabel (int, int, int, int);
 ID_CELL IDCELL();
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

///////////////////////////////////////////////////////////////
struct volumeField 
{
 FVcell** mesh;
 int dimx;
 int dimy;

  volumeField operator+(const volumeField& rhs);
  volumeField operator-(const volumeField& rhs);
  volumeField operator*(double scalar_val);

};
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////

volumeField buildFVfield (int  , int  , double  , double  );

////////////////////////////////////////////////////////////////

void initFieldVal (volumeField  , double  );

////////////////////////////////////////////////////////////////

void updateFieldVal (volumeField*  , volumeField  );

////////////////////////////////////////////////////////////////
void setBCuFV (volumeField*  );
////////////////////////////////////////////////////////////////
void setBCvFV (volumeField*  );
////////////////////////////////////////////////////////////////
void setBCpFV (volumeField*  );
//////////////////////////////////////////////////////////////////////
void divergence (volumeField*  , const volumeField   , const volumeField   );
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void setBCneumannFV (volumeField*  );
void setBCerrFV (volumeField*  );
//////////////////////////////////////////////////////////////////////
void correct(volumeField*  , volumeField*  , const volumeField   , const volumeField   , const volumeField  , const double  , const double  );
//////////////////////////////////////////////////////////////////////////////
void interpolateFieldVal(const volumeField&  , volumeField*   );
//////////////////////////////////////////////////////////////////////////////
void prolongateFieldVal(const volumeField&  , volumeField*   );
///////////////////////////////////////////////////////////////////
double l1_norm(volumeField const  );
///////////////////////////////////////////////////////////////////
double sum (const volumeField&  );

#endif