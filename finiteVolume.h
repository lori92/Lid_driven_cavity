#include <iostream>
#include "mesh.h"

///////////////// finite volume discretization of explicit convective terms //////////
void convection_u_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const double rho);

void convection_v_expl (volumeField* conv, const volumeField& Ux, const volumeField& Uy, const double rho);

///////////////// finite volume discretization of explicit diffusive terms //////////
void diffusion_u_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, const double mu);

void diffusion_v_expl (volumeField* diff, const volumeField& Ux, const volumeField& Uy, double mu);