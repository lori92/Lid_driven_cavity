/*
 * write_paraview_output.c
 *
 *  Created on: Feb 6, 2023
 *      Author: Lorenzo Sufrà
 */

#include <stdio.h>
#include <stdlib.h>
//#include "finiteVolume.C"

void write_output(const volumeField& vel_x, const volumeField& vel_y, const volumeField& vel_z, const volumeField& p, std::string namefile)
{

    FILE *fid;
    int Nx = vel_x.dimx;
    int Ny = vel_x.dimy; 
    int Nz = vel_x.dimz;  

    std::string file = namefile + ".vtk";
    fid = fopen (file.c_str(), "w");
    fprintf (fid, "# vtk DataFile Version 2.0\n");
    fprintf (fid, "Solution file\n");
    fprintf (fid, "ASCII\n");
    fprintf (fid, "DATASET STRUCTURED_GRID\n");
    fprintf (fid, "DIMENSIONS %d %d %d\n",Nx,Ny,Nz);

    fprintf (fid,"POINTS %d double\n",Nx*Ny*Nz);
    for (int i=1;i<=Nx;i++)
    {
      for (int j=1;j<=Ny;j++)
      {
              for (int k=1;k<=Nz;k++)
      {
        fprintf(fid,"%16.7e %16.7e %16.7e\n",vel_x.mesh[i][j][k].xCentroid(), 
                                             vel_y.mesh[i][j][k].yCentroid(),
                                             vel_z.mesh[i][j][k].zCentroid());
      };
      };
    };

    fprintf(fid,"POINT_DATA %d\n",(Nx)*(Ny)*(Nz));

    //fprintf(fid,"CELL_DATA %d\n",(Nx)*(Ny)*(Nz));

    // Write VEL to the file---------------------------
    fprintf (fid,"VECTORS VEL double\n");


    for (int i=1; i<=Nx; i++)
    {
        for (int j=1; j<=Ny; j++)
        {
           for (int k=1; k<=Nz; k++)
           {
            fprintf (fid,"%16.7e %16.7e %16.7e\n", 0.5 * (vel_x.mesh[i][j][k].cellVal() + vel_x.mesh[i-1][ j ][k].cellVal()),
                                                   0.5 * (vel_y.mesh[i][j][k].cellVal() + vel_y.mesh[ i ][j-1][k].cellVal()),
                                                   0.5 *(vel_z.mesh[i][j][k].cellVal() + vel_z.mesh[ i ][j][k-1].cellVal()));
           };
        };
    };

    
    // Write PRES to the file---------------------------
    fprintf (fid,"SCALARS PRES double 1\n");
    fprintf (fid,"LOOKUP_TABLE default\n");

    for (int i=1;i<=Nx;i++)
    {
        for (int j=1;j<=Ny;j++)
        {
            for (int k=1; k<=Nz; k++)
            {
              fprintf (fid,"%16.7e\n", p.mesh[i][j][k].cellVal());
            }
        };
    };
    
    fclose(fid);

}