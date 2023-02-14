#ifndef FVCELL_H
#define FVCELL_H

class FVcell
{       
 private:
  double xC;
  double yC;
  double zC;

  double xFaces[2];
  double yFaces[2];
  double zFaces[2];

  double value;

 public:

 void setCellPos   (double, double, double, double, double, double);
 void setVal (double);
 void print_cellCentre ();
 void print_cellVal ();
 double cellVal ();
 double dx ();
 double dy ();
 double dz ();

 double xCentroid ();
 double yCentroid ();
 double zCentroid ();

 double xFace (int);
 double yFace (int) ;
 double zFace (int) ;

};

///////////////////////////////////////////////////////////////
struct volumeField 
{
 FVcell*** mesh;
 int dimx;
 int dimy;
 int dimz;

  volumeField operator+(const volumeField& rhs);
  volumeField operator-(const volumeField& rhs);
  volumeField operator*(double scalar_val);

};

#endif