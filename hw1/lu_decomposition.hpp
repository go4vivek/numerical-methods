#ifndef LU_DECOMPOSITION_HPP_
#define LU_DECOMPOSITION_HPP_

#include <Eigen/Dense>
#include <boost/tuple/tuple.hpp>

using namespace Eigen;

VectorXd forward_subst(MatrixXd L, VectorXd b){
	/*
	@summary: Forward substitution
	@param L: nonsingular lower triangular matrix of size n
	@param b: column vector of size n
	@return x: solution to Lx=b 
	*/
	int n = sqrt(L.size());
	VectorXd x(n);
	x(0) = b(0) / L(0,0);
	for (int j = 1; j < n; ++j){
		double sum = 0;
		for (int k = 0; k < j; ++k)
			sum += L(j,k) * x(k);
		x(j) = (b(j) - sum) / L(j,j);
	}
	return x;
}

VectorXd backward_subst(MatrixXd U, VectorXd b){
	/*
	@summary: Backward substitution
	@param U: nonsingular upper triangular matrix of size n
	@param b: column vector of size n
	@return x: solution to Ux=b 
	*/
	int n = sqrt(U.size());
	VectorXd x(n);
	x(n-1) = b(n-1) / U(n-1,n-1);
	for (int j = n-2; j >= 0; --j){
		double sum = 0.0;
		for (int k = j+1; k < n; ++k)
			sum += U(j,k) * x(k);
		x(j) = (b(j) - sum) / U(j,j);
	}
	return x;
}

boost::tuple<MatrixXd, MatrixXd> lu_no_pivoting(MatrixXd A) {
	/*
	@summary: LU decomposition without pivoting
	@param A: nonsingular matrix of size n with LU decomposition
	@return L: lower triangular matrix with entries 1 on main diagonal
	@return U: upper triangular matrix such that A=LU 
	*/
	int n = sqrt(A.size());
	MatrixXd L(n,n), U(n,n);
	L.setZero();
	U.setZero();
	for (int i = 0; i < n-1; ++i) {
		for (int k = i; k < n; ++k) {
			U(i,k) = A(i,k);
			L(k,i) = A(k,i) / U(i,i);
		}
		A.block(i+1,i+1,n-i-1,n-i-1) -= L.block(i+1,i,n-i-1,1)*U.block(i,i+1,1,n-i-1);
	}
	L(n-1,n-1) = 1;
	U(n-1,n-1) = A(n-1, n-1);
	return boost::make_tuple(L, U);
}

VectorXd linear_solve_LU_no_pivoting(MatrixXd A, VectorXd b) {
	/*
	@summary: linear solver using LU decomposition without pivoting
	@param A: nonsingular square matrix of size n with LU decomposition
	@param b: column vector of size n
	@return x: solution to Ax=b
	*/
	int n = sqrt(A.size());
	MatrixXd L(n,n), U(n,n);
	boost::tie(L, U) = lu_no_pivoting(A);
	VectorXd y(n), x(n);
	y = forward_subst(L, b);
	x = backward_subst(U, y);
	return x;
}

boost::tuple<VectorXi, MatrixXd, MatrixXd> lu_row_pivoting(MatrixXd A) {
	/*
	@summary: LU decomposition with row pivoting
	@param A: nonsingular matrix of size n
	@return P: permutation matrix, stored as vector of its diagonal entries
	@return L: lower triangular matrix with entries 1 on main diagonal
	@return U: upper triangular matrix such that PA=LU
	*/
	int n = sqrt(A.size());
	VectorXi P = VectorXi::LinSpaced(n,1,n);	// initialize P as an identity matrix
	MatrixXd L(n,n), U(n,n);
	L.setIdentity();							// initialize L as an identity matrix
	for (int i = 0; i < n; ++i){
		ArrayXd::Index i_row, i_column;
		A.block(i,i,n-i,1).array().abs().maxCoeff(&i_row, &i_column);
		ArrayXd::Index i_max = i_row + i;		// find i_max, index of largest entry in absolute value from vector A(i:n,i)
		A.row(i).swap(A.row(i_max));			// switch rows i and i_max of A
		P.row(i_max).swap(P.row(i));					// update matrix P
		if (i > 0)
			L.block(i,0,1,i).swap(L.block(i_max,0,1,i));		// switch rows i and i_max of L
		for (int j = i; j < n; j++) {
			L(j,i) = A(j,i) / A(i,i);			// compute column i of L
			U(i,j) = A(i,j);					// compute row i of U
		}
		A.block(i+1,i+1,n-i-1,n-i-1) -= L.block(i+1,i,n-i-1,1)*U.block(i,i+1,1,n-i-1);
	}
	L(n-1,n-1) = 1;
	U(n-1,n-1) = A(n-1,n-1);
	return boost::make_tuple(P,L,U);
}

#endif /* LU_DECOMPOSITION_HPP_ */
