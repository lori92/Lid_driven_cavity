#ifndef FVCELL_H
#define FVCELL_H

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

#endif