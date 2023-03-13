#include <iostream>
#include <cmath>

#include "mesh.h"


inline int return_n(int, int, int, int);
inline int return_i(int, int, int);
inline int return_j(int, int, int);

extern "C" void dgesv_( int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );

extern "C" void dgetrf_(int* M, int *N, double* A, int* lda, int* ipiv, int* INFO);
//// generate inverse of a matrix given its LU decomposition
//void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int*    lwork, int* INFO);
extern "C" void dgetrs_(char* C, int* N, int* NRHS, double* A, int* lda, int* ipiv, double* B, int* LDB, int* INFO);


/****************************************/

double** assemble_matrix(int, int, 
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


inline int return_i(int n,   int N_x, int N_y)
{
  return n/N_y;
}
 

inline int return_j(int n, int N_x, int N_y)
{
  return n%N_y;
}


inline int return_n(  int i, int j, int N_x, int N_y)
{
    return (i)*N_y + j;
}

/****************************************/

void assemble_coeff( const volumeField& p,
                         double** aP,
                         double** aN,
                         double** aS,
                         double** aW,
                         double** aE)
{
 int N_x = p.dimx;
 int N_y = p.dimy;

 for (int i = 1; i<=N_x; i++)
 {
  for (int j = 1; j<=N_y; j++)
  {
     aP[i][j] = - 1./ ( (p.mesh[i+1][j].xCentroid() - p.mesh[ i ][j].xCentroid()) * 
                         p.mesh[i][j].dx() )
                - 1./ ( (p.mesh[ i ][j].xCentroid() - p.mesh[i-1][j].xCentroid()) * 
                         p.mesh[i][j].dx() )
                - 1./ ( (p.mesh[i][j+1].yCentroid() - p.mesh[i][ j ].yCentroid()) * 
                         p.mesh[i][j].dy() )
                - 1./ ( (p.mesh[i][ j ].yCentroid() - p.mesh[i][j-1].yCentroid()) * 
                         p.mesh[i][j].dy() );

     aW[i][j] =  1./ ( (p.mesh[ i ][j].xCentroid() - p.mesh[i-1][j].xCentroid()) * 
                        p.mesh[i][j].dx() );
     aE[i][j] =  1./ ( (p.mesh[i+1][j].xCentroid() - p.mesh[ i ][j].xCentroid()) * 
                        p.mesh[i][j].dx() );

     aS[i][j] =  1./ ( (p.mesh[i][ j ].yCentroid() - p.mesh[i][j-1].yCentroid()) * 
                        p.mesh[i][j].dy() );
     aN[i][j] =  1./ ( (p.mesh[i][j+1].yCentroid() - p.mesh[i][ j ].yCentroid()) * 
                        p.mesh[i][j].dy() );
  }
 }

 return;
}

/****************************************/
void assemble_matrix(    double**A_,
                         int N_x, int N_y, 
                         double** aP,
                         double** aS,
                         double** aN,
                         double** aW,
                         double** aE,
                         const volumeField& p)
{
 // number of total physical nodes in the mesh
 int N = N_x*N_y;

 // std:: cout << "assembling A_ .."<<"\n";
    // A_ IS THE MATRIX COMPRISING THE GHOST NODES BEYOND THE BOUNDARIES
 
 double** A_tmp = new double* [N+1];
 for (int n=0;n<N+1;n++)
 {
    A_tmp[n] = new double[N+1];
 }
 
 for (int m = 0; m<N+1; m++)
 {
   for (int n = 0; n<N+1; n++) 
   {
       A_[m][n] = 0;
   }
 }

 for (int i=1; i<=N_x; i++){
   for (int j=1; j<=N_y; j++){

    int n = return_n(i-1,j-1,N_x,N_y);

    // if n corresponds to an internal cell
    if (p.mesh[i][j].IDCELL().INTERNAL ==1)
    {
      A_[ n ][ n ] = aP[i][j];
      A_[ n ][n-1] = aS[i][j];
      A_[ n ][n+1] = aN[i][j];
      A_[ n ][n-N_y]= aW[i][j];
      A_[ n ][n+N_y]= aE[i][j];
    }
   
    // if n does not correspond to an internal cell: check on which boundary I am
    else if (p.mesh[i][j].IDCELL().INTERNAL ==0)
    {
      A_[ n ][ n ] = aP[i][j];

      if (p.mesh[i][j].IDCELL().NORD ==1)   
      {
       A_[ n ][ n ] = A_[ n ][ n ] + aN[i][j] ; 
      }

      if (p.mesh[i][j].IDCELL().SUD ==1)  
      { 
       A_[ n ][ n ] = A_[ n ][ n ] + aS[i][j];
      }

      if (p.mesh[i][j].IDCELL().EAST ==1)  
      { 
       A_[ n ][ n ] = A_[ n ][ n ] + aE[i][j];
      }

      if (p.mesh[i][j].IDCELL().WEST ==1)  
      { 
       A_[ n ][ n ] = A_[ n ][ n ] + aW[i][j];
      }

      // assemble anyway the ghost slots, which will
      // be discarded when the final matrix will be
      // obtained by restricting A_ 
      
      if (!p.mesh[i][j].IDCELL().SUD  ){ A_[ n ][n -  1]= (aS[i][j]); }
      if (!p.mesh[i][j].IDCELL().NORD ){ A_[ n ][n +  1]= (aN[i][j]); }
      if (!p.mesh[i][j].IDCELL().WEST ){ A_[ n ][n -N_y]= (aW[i][j]); }
      if (!p.mesh[i][j].IDCELL().EAST ){ A_[ n ][n +N_y]= (aE[i][j]); }
	 }
   }
 }

 int n = N-1;
 A_[ n ][ n ] = A_[ n ][ n ] - aN[N_x][N_y];
 A_[ n ][n+1] = aN[N_x][N_y];


 // BOUNDARY NODE AT n = N set to atmosferic
/* 
 for (int ROW = 0; ROW<N; ROW++)
 {
  A_[ROW][N] = 0.;
 } 
 
 for (int COL = 0; COL<N; COL++)
 {
  A_[N][COL] = 0.;
 } 
 */
  A_[N][N] = 1.;
  A_[N][N-1] = 1.;
 


 FILE* fid = fopen("A_.dat","w");
 if (fid)
 {
    for (int i = 0; i<N+1; i++)
    {
     for (int j = 0; j<N+1; j++)
     {
       A_tmp[i][j] = A_[i][j];
       fprintf(fid, "%16.7e    ", A_[i][j]);
     }       
     fprintf(fid, "\n");
    }
  fclose(fid);
 }
 else 
 {
    std::cout<<"\n";
    std::cout<<"************ ERROR: matrix A_ not written.. ************"<<"\n";
    std::cout<<"\n";
 };

    std::cout<<"transposing A ...\n";

  for (int i = 0; i<N+1; i++)
  {
    for (int j = 0; j<N+1; j++)
    {
      A_[j][i] = A_tmp[i][j];
    }       
  }


 return;
}


/****************************************/

double* forcing_term(const volumeField& b, 
                         double** aP,
                         double** aN,
                         double** aS,
                         double** aW,
                         double** aE)
{
 int N_x = b.dimx;
 int N_y = b.dimy;

 int N = N_x*N_y;

 double* f = new double[N+1];
 
 for (int i=0; i<N_x; i++){
   for (int j=0; j<N_y; j++){

      int n = return_n(i,j,N_x,N_y);

      f[n] = b.mesh[i+1][j+1].cellVal();

      if (b.mesh[i+1][j+1].IDCELL().EAST ==1)   
      {
       //f[n] = f[n] - 2*5*aE[i+1][j+1] ; 
      }
	  
      if (b.mesh[i+1][j+1].IDCELL().WEST ==1)   
      {
       //f[n] = f[n] - 2*15*aW[i+1][j+1] ; 
      }
	  
	
   }
 }

 f[N] =   2.*101325.; 

 FILE* fid = fopen("b_.dat","w");
 if (fid)
 {
    
     for (int j = 0; j<N+1; j++)
     {
       fprintf(fid, "%16.7e    ", f[j]);
     }       
     fprintf(fid, "\n");
    
  fclose(fid);
 }
 else 
 {
    std::cout<<"\n";
    std::cout<<"************ ERROR: matrix A_ not written.. ************"<<"\n";
    std::cout<<"\n";
 };

 


 return f;
}
/****************************************/

void forcing_term2( double* f, const volumeField& b, 
                         double** aP,
                         double** aN,
                         double** aS,
                         double** aW,
                         double** aE)
{
 int N_x = b.dimx;
 int N_y = b.dimy;

 int N = N_x*N_y;

 
 for (int i=0; i<N_x; i++){
   for (int j=0; j<N_y; j++){

      int n = return_n(i,j,N_x,N_y);

      f[n] = b.mesh[i+1][j+1].cellVal();

      if (b.mesh[i+1][j+1].IDCELL().EAST ==1)   
      {
       //f[n] = f[n] - 2*5*aE[i+1][j+1] ; 
      }
	  
      if (b.mesh[i+1][j+1].IDCELL().WEST ==1)   
      {
       //f[n] = f[n] - 2*15*aW[i+1][j+1] ; 
      }
	  
	
   }
 }

 f[N] =   2.*101325.; 

 FILE* fid = fopen("b_.dat","w");
 if (fid)
 {
    
     for (int j = 0; j<N+1; j++)
     {
       fprintf(fid, "%16.7e    ", f[j]);
     }       
     fprintf(fid, "\n");
    
  fclose(fid);
 }
 else 
 {
    std::cout<<"\n";
    std::cout<<"************ ERROR: matrix A_ not written.. ************"<<"\n";
    std::cout<<"\n";
 };

 


 return;// f;
}

/**************************************/

void reshape(int N, double* array, double** matrix )
{
    //double array[N];
    int k = 0;

    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
        {
            *(array+k)= matrix[i][j];
                        k = k+1;
        } 
    }
 return;
};

/**************************************/

void getVolumeField(volumeField* field,  double* array)
{
   int N_x = field->dimx;
   int N_y = field->dimy;

   int n;

    for (int i=0;i<N_x;i++)
    {
        for (int j=0;j<N_y;j++)
        {
         n = return_n (i,j,N_x,N_y);
         if (n>=N_x*N_y)
         {std::cout<<"ERROR: n="<< n<<" >N"<<"\n";}
         field->mesh[i+1][j+1].setVal(array[n]);
        } 
    }
 return;
};

/**************************************/

void solvePressure( double** LaplaceOperator,
                    double L[],
                    volumeField* p, 
                    volumeField div, 
                    const double dt, 
                    const double rho,
                    double** aP,
                    double** aN,
                    double** aS,
                    double** aW,
                    double** aE ) 
{
   int N_x = p->dimx;
   int N_y = p->dimy;
   int dim = N_x*N_y + 1;
   int nrhs = 1;
   int lda  = dim;
   int ipiv[dim-1];
   int ldb = dim;
   int info;



   // reorder in a vector all the elements from the field div*(rho/dt)
   double* b_reshaped;
   b_reshaped          = new double [dim];
   b_reshaped = forcing_term( div*(rho/dt) ,aP,aS,aN,aW,aE );

   // reshape the laplacian operator in a dimxdim vector so as to have contiguos allocation of memor
   // (needed for dgesv routine)    
   reshape(dim, &L[0], LaplaceOperator);
   
   std:: cout<<"call dgesv_"<<"\n";
   // solve for the pressure vectorp p: L p = b and impose BC on the reshapded field
   dgesv_(&dim, &nrhs,  &L[0], &lda, ipiv, b_reshaped, &ldb, &info);
   getVolumeField(p,  b_reshaped);  
   setBCpFV(p);

   //delete(L);
   delete(b_reshaped);

   return;
}

/**************************************/

bool tdma(int N, const double a[], const double b[], const double c[], const double d[], 
           double X[] )
//*********************************************************************************
// Solve, using the Thomas Algorithm (TDMA), the tri-diagonal system              *
//     a_i X_i-1 + b_i X_i + c_i X_i+1 = d_i,     i = 0, N - 1                    *
//                                                                                *
// Effectively, this is the N x N matrix equation.                                *
// a[i], b[i], c[i] are the non-zero diagonals of the matrix and d[i] is the rhs. *
// a[0] and c[n-1] aren't used.                                                   *
//*********************************************************************************
{
   
   double P[N];
   double Q[N];
   X = P;

   // Forward pass
   int i = 0;
   double denominator = b[i];
   P[i] = -c[i] / denominator;
   Q[i] =  d[i] / denominator;
   for ( i = 1; i < N; i++ )
   {
      denominator = b[i] + a[i] * P[i-1];
      if ( fabs( denominator ) < 1e-30 ) return false;
      P[i] =  -c[i]                   / denominator;
      Q[i] = ( d[i] - a[i] * Q[i-1] ) / denominator;
   }

   // Backward pass
   X[N-1] = Q[N-1];
   for ( i = N - 2; i >= 0; i-- ) X[i] = P[i] * X[i+1] + Q[i];
   
   return true;

}

/**************************************/

void solvePressure_test( 
                    double*  L,
                    volumeField* p, 
                    volumeField div, 
                    const double dt, 
                    const double rho,
                    double** aP,
                    double** aN,
                    double** aS,
                    double** aW,
                    double** aE,
                    int& factorization,
                    int& info,
                    int ipiv[] ) 
{
   int N_x = p->dimx;
   int N_y = p->dimy;
   int dim = N_x*N_y + 1;
   int nrhs = 1;
   int lda  = dim;
   //int ipiv[dim];
   int ldb = dim;
   //int info;
   char trans = 'N';


   // reorder in a vector all the elements from the field div*(rho/dt)
   //double* b_reshaped;
   //b_reshaped          = new double [dim];
   //forcing_term2(b_reshaped, div*(rho/dt) ,aP,aS,aN,aW,aE );

   double* b_reshaped;
   b_reshaped          = new double [dim];
   b_reshaped = forcing_term( div*(rho/dt) ,aP,aS,aN,aW,aE );
   // reshape the laplacian operator in a dimxdim vector so as to have contiguos allocation of memor
   // (needed for dgesv routine)    
   //reshape(dim, &L[0], LaplaceOperator);
   
   //std:: cout<<"call dgesv_"<<"\n";
   // solve for the pressure vectorp p: L p = b and impose BC on the reshapded field
   //dgesv_(&dim, &nrhs,  &L[0], &lda, ipiv, b_reshaped, &ldb, &info);
  // dgetrs_(&trans, &dim, &nrhs, &L[0], &lda, ipiv, b_reshaped, &ldb, &info);


  /* if no factorization is given, calculate it */
  if (factorization == 0)
  {

    dgetrf_(&dim, &dim, &L[0], &dim, ipiv, &info);
    factorization = 1;
 

    dgetrs_(&trans, &dim, &nrhs, &L[0], &lda,  ipiv, b_reshaped, &ldb, &info);
 
  }
  else
  {

 

     dgetrs_(&trans, &dim, &nrhs, &L[0], &lda,  ipiv, b_reshaped, &ldb, &info);
 
  }
   getVolumeField(p,  b_reshaped);  
   setBCpFV(p);

   //delete(L);
   //delete(b_reshaped);('N',neqn,1,coeff,neqn,ipiv,lhs,neqn,info2)

   return;
}

/**************************************/

