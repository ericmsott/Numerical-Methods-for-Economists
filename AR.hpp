#include <iostream>
#include </usr/include/eigen3/Eigen/Dense> // Para implementar Alg Lin.
#include "Tauchen.hpp"

using namespace std;
using Eigen::MatrixXd; // matriz dinamica do pacote 'Eigen'.
using Eigen::VectorXd; // vetor dinamico do pacote 'Eigen'.
using namespace Eigen;

VectorXd residuals(double mu,
		   double sigma,
		   const int T)
{
    VectorXd e(T + 100);

    for(int t = 0; t <= (T - 1); t++)
    {
        e(t) = normDraw(mu,sigma,t);
    }

    return e;
}

VectorXd AR1(double mu,
	     double rho,
	     double sigma,
	     const int T,
	     VectorXd e)
{
    VectorXd Y(T + 100);
    Y(0) = mu;

    for(int t = 1; t <= (T+99); t++)
    {
      Y(t) = (1 - rho)*mu + rho*Y[t-1] + e[t];
    }

    return Y.tail(T);
}

VectorXd MarkovChainSim(double sigma,
			double rho,
			const int N,
			double mu,
			double m,
			const int T,
			VectorXd e,
			MatrixXd P)
{
    VectorXd u = normCDFvec(e/sigma);

    double theta_N = m*sigma/sqrt(1-rho*rho);
    double theta_1 = -theta_N;
    VectorXd grid(N);

    for(int z = 0; z < N; z++)
    {  // Grid Linear
        grid(z) = (z + 1)*((theta_N -theta_1)/N) + theta_1;
    }

    VectorXd Y(T+100);
    int state = ceil(N/2) -1;

    Y(0) = grid(state);

    VectorXd csum(N);
    int index;

    for(int t=1; t < (T+100); t++)
    {
        csum = P.row(state).colwise().sum();
	index = 0;
	for(int j=0; j<(N-1); j++)
	{
	  if(u(t) >= csum(j) && u(t) < csum(j+1)  && index <= (N-2))
	    {
	        index = index + 1;
	    }
	}
	state = index + 1;
        Y(t) = grid(state);
    }

    return Y.tail(T);
}
