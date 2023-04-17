#ifndef FVCELL_H
#define FVCELL_H

struct ID_CELL
{ 
   int INTERNAL;
   int NORD;
   int SUD;
   int EAST;
   int WEST;
};
///////////////////////////////////////////////////////////////
class FVcell
{       
 private:
  double xC;
  double yC;
  double xFaces[2];
  double yFaces[2];
  ID_CELL ID;

 public:
  void setCellPos   (double, double, double, double);
  void setCellLabel (int, int, int, int, int,int);
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
class FVmesh 
{
   protected:
    int dimx;
    int dimy;
    int dimx_tot;
    int dimy_tot;
    double Lx;
    double Ly;

   private:
    static std::vector<FVcell> mesh;
  
   public:
    FVmesh(int, int, int, int, double, double, int);
    void writeFVmesh (int);
    FVcell getFVcell(int, int);
    int get_dimx();
    int get_dimy();
    int get_dimx_tot();
    int get_dimy_tot();

};

///////////////////////////////////////////////////////////////
class FVfield:   public FVmesh
{
   private:
    std::vector<double> cellValue;

   public:
    void allocateFVfield();
    void setVal(int, int, double);
};
///////////////////////////////////////////////////////////////

#endif