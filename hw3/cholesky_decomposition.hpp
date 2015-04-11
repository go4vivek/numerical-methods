#ifndef CHOLESKY_DECOMPOSITION_HPP_
#define CHOLESKY_DECOMPOSITION_HPP_

#include <Eigen/Dense>
#include <boost/tuple/tuple.hpp>
#include "lu_decomposition.hpp"

using namespace Eigen;
using namespace std;

MatrixXd cholesky(MatrixXd A){
	int n = sqrt(A.size());
	MatrixXd U(n,n);
	U.setZero();
	for (int i = 0; i <= n-2; i++) {
		U(i,i) = sqrt(A(i,i));
		for (int k = i+1; k <= n-1; k++)
			U(i,k) = A(i,k) / U(i,i);
		for (int j = i+1; j <= n-1; j++)
			for (int k = j; k <= n-1; k++)
				A(j,k) = A(j,k) - U(i,j)*U(i,k);
	}
	U(n-1,n-1) = sqrt(A(n-1,n-1));
	return U;
}

VectorXd linear_solve_cholesky(MatrixXd A, VectorXd b){
	int n = sqrt(A.size());
	MatrixXd U = cholesky(A);
	VectorXd y = forward_subst(U.transpose(), b);
	return backward_subst(U, y);
}

#endif /* CHOLESKY_DECOMPOSITION_HPP_ */
