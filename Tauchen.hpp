#include <iostream>
#include </usr/include/eigen3/Eigen/Dense>
#include "Normal.hpp"

using namespace std;
using Eigen::MatrixXd;

/* A função Tauchen é a principal deste
   header file, e ela implementa o método de Tauchen. */

MatrixXd tauchen(double     sigma,
		 double     rho,
		 const int  N,
		 double     mu,
		 double     m=1)
{
    /* Step 1:
        "Compute the points of the grid"  */
    double theta_N = m*sigma/sqrt(1 - rho*rho);
    double theta_1 = -theta_N;
    double grid[N];
    for(size_t z{0}; z < N; z++)
    {  // Grid Linear
        grid[z] = (z + 1)*((theta_N -theta_1)/N) + theta_1;
    }
    /* Para definir o intervalo D_theta
       a ser utilizado no cômputo da CDF adiante: */
    double D_theta = grid[2] - grid[1];

    /* Step 2:
        "Compute each p_ij of the matrix"  */

    MatrixXd P(N,N); // Aloca espaço na memória para a matriz P.

    for(size_t i{0}; i < N; i++)
    {
        for(size_t j{0}; j < N; j++)
        {
            if(j==1)        // first corner points
            {
	      P(i,j) = normCDF((theta_1 + D_theta/2 - (1 - rho)*mu - rho*grid[i])/sigma);
            } else if(j==N) // last corner points
            {
	      P(i,j) = 1 - normCDF((grid[j] - D_theta/2 - (1 - rho)*mu - rho*grid[i])/sigma);
            } else          // middle points
            {
	      P(i,j) = normCDF((grid[j] + D_theta/2 - (1 - rho)*mu - rho*grid[i])/sigma)
		     - normCDF((grid[j] - D_theta/2 - (1 - rho)*mu - rho*grid[i])/sigma);
            }
        }
    }
    /* Now the function returns the desired markov transition matrix.*/
    return P;
}
