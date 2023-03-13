#include "mesh.h"


void predictor_u(volumeField*  , 
                  const volumeField&  , 
                  const volumeField&  , 
                   const volumeField&  , 
                  const volumeField&  , 
                  const double  , const  double  , 
                  const double  ,  const double  , 
                  int&   );
////////////////////////////////////////////////////////////////////////////////////////////////////
void  predictor_v(volumeField*  , 
                  const volumeField&  , 
                  const volumeField&  , 
                   const volumeField&  , 
                  const volumeField&  , 
                  const double  , const  double  , 
                  const double  ,  const double  , 
                  int&   );
////////////////////////////////////////////////////////////////////////////////////////////////////