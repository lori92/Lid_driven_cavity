#include <iostream>
#include "mesh.h"
//#include <petscksp.h>




inline int return_n(int, int, int, int);
inline int return_i(int, int, int);
inline int return_j(int, int, int);

extern "C" void dgesv_( int * , int * , double  * , int * , int * , double * , int * , int *   );

extern "C" void dgetrf_(int*  , int * , double*  , int*  , int*  , int*  );

extern "C" void dgetrs_(char*  , int*  , int*  , double*  , int*  , int*  , double*  , int*  , int*  );

/****************************************/

double** assemble_matrix(double**,
                         int, int, 
                         double**,
                         double**,
                         double**,
                         double**,
                         double**,
                         const volumeField& );

/****************************************/

void assemble_coeff(const volumeField&, 
                         double**,
                         double**,
                         double**,
                         double**,
                         double**);

/****************************************/

void reshape(int , double*  , double**   );

/****************************************/

double* forcing_term(const volumeField&  , 
                         double** ,
                         double**  ,
                         double**  ,
                         double**  ,
                         double**  );

/****************************************/

void updateForcingTerm(  double [],
                         const volumeField&);

/**************************************/

void reshape(int  , double*  , double**   );

/**************************************/

void getVolumeField(volumeField*  ,  double*  );

/**************************************/

void getArray(int, double[], const FVcell**, int );

/**************************************/


void solvePressure( double**  ,
                    double [] , 
                    double [],
                    volumeField*  , 
                    volumeField  , 
                    const double  , 
                    const double  ,
                    double**  ,
                    double**  ,
                    double**  ,
                    double**  ,
                    double**   );

/**************************************/

bool tdma(int  , const double [], const double [], const double [], const double [], 
           double [] );

/**************************************/

void solvePressure_test( 
                    double []   ,
                    double [],
                    volumeField* , 
                    volumeField  , 
                    const double , 
                    const double ,
                    int&,int&, int  [] );

/**************************************/

 void getMPI_Buffer(int  , buffer*  , FVcell**   , int,  int  , int  );

/**************************************/

 void getMesh      (int  , buffer*  ,  volumeField*,
                                        volumeField*,
                                        volumeField*,
                                        volumeField*,
                                        volumeField*, 
                                        volumeField*, int);



