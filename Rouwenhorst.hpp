#include <iostream>
#include </usr/include/eigen3/Eigen/Dense> // Para implementar Alg Lin utilizando o LAPACK.

using namespace std;
using Eigen::MatrixXd; // namespace do pacote 'Armadillo'.

/* A função rouwenhorst é a principal deste
   header file, e ela implementa o método de Rouwenhorst. */

MatrixXd rouwenhorst(double sigma, double rho, const int N)
{
    double p = (1 + rho)/2;
    MatrixXd P(2,2); // Aloca espaço na memória para a matriz P.
    MatrixXd eye = Matrix<double, 2, 2>::Identity();
    P = p*eye + (1-p)*eye.transpose(); // Constroi a matriz P para N=2

    for(int i=2; i < N; i++)
    {
	MatrixXd P_A = MatrixXd::Zero(i+1,i+1);
	MatrixXd P_B = P_A;
	MatrixXd P_C = P_A;
	MatrixXd P_D = P_A;

	P_A.block(0,0,i,i) =     p*P;
	P_B.block(1,1,i,i) =     p*P;
	P_C.block(1,0,i,i) = (1-p)*P;
	P_D.block(0,1,i,i) = (1-p)*P;

	P = P_A + P_B + P_C + P_D;

	for(size_t j{1}; j < i; j++)
        {
            P.row(j) = P.row(j)/2;
	}
    }
    /* Now the function returns the desired markov transition matrix.*/
    return P;
}
