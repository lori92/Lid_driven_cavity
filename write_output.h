/*
 * write_paraview_output.c
 *
 *  Created on: Feb 6, 2023
 *      Author: Lorenzo Sufrà
 */
#include "mesh.h"

 
void write_output(const volumeField&  , const volumeField&  ,  const volumeField&  , const double  , std::string  );

void write_data(const volumeField&  , const volumeField&  ,  const volumeField&  , const double  );


void read_data(volumeField*  ,  volumeField*  ,   volumeField*  , double&  );

void write_k(const volumeField  , const volumeField  , const double  );

void write_centre_line(const volumeField  , const volumeField  , const double  ) ;

void read_input( const int , int& , double& , double& , double&, double&, double&, double&, double&);
