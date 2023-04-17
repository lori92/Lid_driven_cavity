int n_iter_ssor_p;
int n_iter_ssor_u;
int n_iter_ssor_v;

// restart flag
int restart;

// fluid properties
double mu;
double rho;

// time parameters
double t0;
double tend;
double dt;
double t;
double dt_write;

//geometry parameters
const int Nx_tot = 4;
const int Ny_tot = 4;

double Lx;
double Ly;

//number of ranks 
int num_procs = 4;

