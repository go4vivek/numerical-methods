#include <Eigen/Dense>
#include <iostream>
#include <boost/tuple/tuple.hpp>

using namespace Eigen;
using namespace std;
//using namespace boost;

VectorXd forward_subst(MatrixXd L, VectorXd b){
	int n = b.rows();
	VectorXd x(n);
	x(0) = b(0) / L(0,0);
	for(int j= 1; j< n ; j++){
		double sum = 0;
		for (int k = 0; k< j; k++){
			sum += L(j,k)*x(k);			
		}
		x(j) = (b(j)-sum)/L(j,j);
	}
	return x;
}

VectorXd forward_subst_banded(MatrixXd L, VectorXd b, int m){
	int n = b.rows();
	VectorXd x(n);
	x(0) = b(0) / L(0,0);
	for(int j= 1; j< n ; j++){
		double sum = 0;
		for (int k = fmax(0, j-m); k< j; k++){
			sum += L(j,k)*x(k);			
		}
		x(j) = (b(j)-sum)/L(j,j);
	}
	return x;
}

VectorXd backward_subst(MatrixXd L, VectorXd b){
	int n = b.rows();
	VectorXd x(n);
	x(n-1) = b(n-1) / L(n-1,n-1);
	for(int j= n-2; j>-1 ; j--){
		double sum = 0;
		for (int k = j+1; k<n; k++){
			sum += L(j,k)*x(k);			
		}
		x(j) = (b(j)-sum)/L(j,j);
	}
	return x;
}


VectorXd backward_subst_banded(MatrixXd L, VectorXd b, int m){
	int n = b.rows();
	VectorXd x(n);
	x(n-1) = b(n-1) / L(n-1,n-1);
	for(int j= n-2; j>-1 ; j--){
		double sum = 0;
		int upper = fmin(n, j+m+1);
		for (int k = j+1; k < upper; k++){
			sum += L(j,k)*x(k);			
		}
		x(j) = (b(j)-sum)/L(j,j);
	}
	return x;
}

boost::tuple<MatrixXd,MatrixXd> lu_no_pivoting(MatrixXd A){
	int n = A.rows();
	MatrixXd L(n,n), U(n,n);
	L.setZero();
	U.setZero();
	for(int i=0; i<n-1;i++){
		for(int k=i;k<n;k++){
			U(i,k) = A(i,k);
			L(k,i) = A(k,i)/U(i,i);
		}
		A.block(i+1,i+1,n-i-1,n-i-1) = A.block(i+1,i+1,n-i-1,n-i-1) - L.block(i+1,i, n-i-1, 1) * U.block(i, i+1, 1, n-i-1);
	}
	L(n-1,n-1) = 1; 
	U(n-1,n-1) = A(n-1,n-1);
	return boost::make_tuple(L,U);
}

boost::tuple<MatrixXd,MatrixXd> lu_no_pivoting_banded(MatrixXd A, int m){
	int n = A.rows();
	MatrixXd L(n,n), U(n,n);
	L.setZero();
	U.setZero();
	for (int i=0; i<n-1;i++) {
		int upper = fmin(i+m+1, n);
		for (int k=i; k<upper; k++){
			U(i,k) = A(i,k);
			L(k,i) = A(k,i)/U(i,i);
		}
		for (int j = i+1; j< upper; j++){
			for (int k = i+1; k< upper; k++){
				A(j,k) = A(j,k) - L(j,i) * U(i,k);
			}
		}
	}
	L(n-1,n-1) = 1; 
	U(n-1,n-1) = A(n-1,n-1);
	return boost::make_tuple(L,U);
}

boost::tuple<MatrixXd, MatrixXd, MatrixXd> lu_row_pivoting(MatrixXd A){
	int n = A.rows();
	MatrixXd L(n,n), U(n,n), P(n,n);
	L.setZero();
	U.setZero();
	P.setIdentity();
	ArrayXd::Index max_row, max_col;
	for(int i=0; i<n-1; i++){
		A.block(i,i,n-i,1).array().abs().maxCoeff(&max_row, &max_col);
		int i_max = max_row + i;
		A.row(i).swap(A.row(i_max));
		P.row(i).swap(P.row(i_max));
		if (i>0)
		{
		L.block(i,0,1,i).swap(L.block(i_max,0,1,i));
		}
		for(int j=i;j<n;j++){
			L(j,i) = A(j,i)/A(i,i);
			U(i,j) = A(i,j);
		}
		A.block(i+1,i+1,n-i-1,n-i-1) = A.block(i+1,i+1,n-i-1,n-i-1) - L.block(i+1,i, n-i-1, 1) * U.block(i, i+1, 1, n-i-1);
	}
	L(n-1,n-1) = 1; 
	U(n-1,n-1) = A(n-1,n-1);
	return boost::make_tuple(L,U,P);
}

VectorXd linear_solve_LU_row_pivoting(MatrixXd A, VectorXd b){
	MatrixXd L, U, P;
	VectorXd x,y;
	boost::tie(L,U, P) = lu_row_pivoting(A);
	y = forward_subst(L, P * b);
	x = backward_subst(U, y);
	return x;
}


MatrixXd cholesky(MatrixXd A){
	int n = A.rows();
	MatrixXd U(n,n);
	U.setZero();
	for(int i=0; i<n-1; i++){
		U(i,i) = sqrt(A(i,i));
		for(int j=i+1; j<n; j++){
			U(i,j) = A(i,j)/U(i,i);
		}
		A.block(i+1,i+1,n-i-1,n-i-1) = A.block(i+1,i+1,n-i-1,n-i-1) - U.block(i, i+1, 1, n-i-1).transpose() * U.block(i, i+1, 1, n-i-1);
	}
	U(n-1,n-1) = sqrt(A(n-1,n-1));
	return U;
}

MatrixXd cholesky_banded(MatrixXd A, int m){
	int n = A.rows();
	MatrixXd U(n,n);
	U.setZero();
	for(int i=0; i<n-1; i++){
		U(i,i) = sqrt(A(i,i));
		int upper = fmin(n, i+m+1);
		for(int j=i+1; j<upper; j++){
			U(i,j) = A(i,j)/U(i,i);
		}
		for (int j = i+1; j < upper; j++){
			for (int k = j; k < upper; k++){
				A(j,k) -= U(i,j) * U(i,k);
			}
		}
	}
	U(n-1,n-1) = sqrt(A(n-1,n-1));
	return U;
}

boost::tuple<MatrixXd,MatrixXd> lu_no_pivoting_tridiag(MatrixXd A){   // Op count 3m + O(1) , assignment doesn't count, only + - * / sqrt counts
	int n = A.rows();
	MatrixXd L(n,n), U(n,n);
	L.setZero();
	U.setZero();
	for(int i=0; i<n-1; i++){
		U(i,i) = A(i,i);
		U(i,i+1) = A(i,i+1);
		L(i,i) = 1;
		L(i+1,i) = A(i+1,i)/U(i,i);
		A(i+1,i+1) = A(i+1,i+1) - L(i+1,i)*U(i,i+1);
	}
	L(n-1,n-1) = 1; U(n-1,n-1) = A(n-1,n-1);
	return boost::make_tuple(L,U);
}

VectorXd linear_solve_LU_no_pivoting_tridiag(MatrixXd A, VectorXd b){
	MatrixXd L, U;
	VectorXd x,y;
	boost::tie(L,U) = lu_no_pivoting_tridiag(A);
	y = forward_subst(L, b);
	x = backward_subst(U, y);
	return x;
}

VectorXd linear_solve_LU_no_pivoting(MatrixXd A, VectorXd b){
	MatrixXd L, U;
	VectorXd x,y;
	boost::tie(L,U) = lu_no_pivoting(A);
	y = forward_subst(L, b);
	x = backward_subst(U, y);
	return x;
}

VectorXd linear_solve_cholesky(MatrixXd A, VectorXd b){
	MatrixXd U;
	VectorXd x, y;
	U = cholesky(A);
	y = forward_subst(U.transpose(), b);
	x = backward_subst(U, y);
	return x;
}

VectorXd linear_solve_cholesky_banded(MatrixXd A, VectorXd b, int m){
	MatrixXd U;
	VectorXd x, y;
	U = cholesky_banded(A, m);
	y = forward_subst_banded(U.transpose(), b, m);
	x = backward_subst_banded(U, y, m);
	return x;
}


MatrixXd natural_cubic_spline_interpolation(VectorXd nodes, VectorXd values){
	MatrixXd results, M;
	int n = nodes.rows() - 1;
	VectorXd w(n+1), q(n), r(n);
	w.setZero();
	results.resize(n, 4);
	M.resize(n-1, n-1);
	M.setZero();
	VectorXd z(n-1);
	for (int i = 0; i < n-1; i++){
		z(i) = 6 * ((values(i+2)-values(i+1)) / (nodes(i+2)-nodes(i+1)) - (values(i+1)-values(i)) / (nodes(i+1)-nodes(i)));
	}
	//cout <<z <<endl;
	for (int i=0; i < n-1; i++) M(i,i) = 2 * (nodes(i+2)-nodes(i));
	for (int i=0; i < n-2; i++) M(i,i+1) = nodes(i+2) - nodes(i+1);
	for (int i=1; i < n-1; i++) M(i,i-1) = nodes(i+1) - nodes(i);

	//cout << M<<endl;

	w.segment(1, n-1) = linear_solve_LU_no_pivoting_tridiag(M, z);

	//cout << w<<endl;

	for (int i = 0; i < n; i++){
		results(i,2) = (w(i)*nodes(i+1) - w(i+1)*nodes(i)) / (2 * (nodes(i+1)-nodes(i)));
		results(i,3) = (w(i+1) - w(i)) / (6 * (nodes(i+1) - nodes(i)));
	}
	for (int i = 0; i < n; i++){
		q(i) = values(i) - results(i,2) * pow(nodes(i), 2) - results(i,3)*pow(nodes(i), 3);
		r(i) = values(i+1) - results(i,2) * pow(nodes(i+1), 2) - results(i,3)*pow(nodes(i+1), 3);
	}
	for (int i = 0; i < n; i++){
		results(i,0) = (q(i)*nodes(i+1) - r(i)*nodes(i)) / (nodes(i+1) - nodes(i));
		results(i,1) = (r(i) - q(i)) / (nodes(i+1) - nodes(i));
	}
	return results;
}

VectorXd least_squares(MatrixXd X, VectorXd b){
	return linear_solve_cholesky(X.transpose()*X, X.transpose()*b);
}

boost::tuple <VectorXd, int> Jacobi_iteration(MatrixXd A, VectorXd b, VectorXd initial, bool rezidual_based=1, double tol=pow(10, -6)){
	VectorXd x, r, b_new;
	MatrixXd M, N;
	int ic = 0;
	M = A.diagonal().asDiagonal();
	N = A - M;
	b_new = forward_subst_banded(M, b, 0);
	
	if(rezidual_based == 1){
		x = initial;
		r = b - A * x;
		double stop_iter_resid = tol * r.norm();
		while(r.norm() > stop_iter_resid){
			x = - forward_subst_banded(M, N*x, 0) + b_new;
			r = b - A * x;
			ic ++;
			//if (ic<=3) cout << x << endl << endl;
		}
	}
	else {
		VectorXd x_old = initial;
		VectorXd diff(initial.rows());
		diff.setOnes();
		while(diff.norm()>tol){
			x = x = - forward_subst_banded(M, N*x_old, 0) + b_new;
			diff = x - x_old;
			x_old = x;
			ic ++;
		}

	}
	return boost::make_tuple(x, ic);
}

boost::tuple <MatrixXd, int> GS_iteration(MatrixXd A, VectorXd b, VectorXd initial, bool rezidual_based=1, double tol=pow(10, -6)){
	VectorXd x, r, b_new;
	MatrixXd M, N;
	int ic = 0;
	M = A.triangularView<Lower>();
	N = A - M;
	b_new = forward_subst(M, b);
	if(rezidual_based == 1){
		x = initial;
		r = b - A * x;
		double stop_iter_resid = tol * r.norm();
		while(r.norm() > stop_iter_resid){
			x = - forward_subst(M, N*x) + b_new;
			r = b - A * x;
			ic ++;
			//if (ic<=3) cout << x << endl << endl;
		}
	}
	else {
		VectorXd x_old = initial;
		VectorXd diff(initial.rows());
		diff.setOnes();
		while(diff.norm()>tol){
			x = - forward_subst(M, N*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic ++;

		}
	}
	return boost::make_tuple(x, ic);
}

boost::tuple <MatrixXd, int> SOR_iteration(MatrixXd A, VectorXd b, VectorXd initial, 
										   double w, bool rezidual_based=1, double tol=pow(10, -6)){
	VectorXd x, r, b_new;
	MatrixXd L, D, U, M, N;
	int ic = 0;
	D = A.diagonal().asDiagonal();
	L = A.triangularView<Lower>();
	L = L - D;
	U = A.triangularView<Upper>();
	U = U - D;
	M = D + w*L;
	N = (1-w)*D - w*U;
	b_new = w * forward_subst(M, b);
	if(rezidual_based == 1){
		x = initial;
		r = b - A * x;
		double stop_iter_resid = tol * r.norm();
		while(r.norm() > stop_iter_resid){
			x = forward_subst(M, N*x) + b_new;
			r = b - A * x;
			ic ++;
			if (ic<=3) cout << x << endl << endl;
		}
	}
	else{
		VectorXd x_old = initial;
		VectorXd diff(initial.rows());
		diff.setOnes();
		while(diff.norm()>tol){
			x = forward_subst(M, N*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic ++;
		}
	}
	return boost::make_tuple(x, ic);
}

boost::tuple <MatrixXd, int> GS_iteration_banded(MatrixXd A, VectorXd b, VectorXd initial, int band, 
												bool rezidual_based=1, double tol=pow(10, -6)){
	VectorXd x, r, b_new;
	MatrixXd M, N;
	int ic = 0;
	M = A.triangularView<Lower>();
	N = A - M;
	b_new = forward_subst_banded(M, b, band);
	if(rezidual_based == 1){
		x = initial;
		r = b - A * x;
		double stop_iter_resid = tol * r.norm();
		while(r.norm() > stop_iter_resid){
			x = - forward_subst_banded(M, N*x, band) + b_new;
			r = b - A * x;
			ic ++;
		}
	}
	else {
		VectorXd x_old = initial;
		VectorXd diff(initial.rows());
		diff.setOnes();
		while(diff.norm()>tol){
			x = - forward_subst_banded(M, N*x_old, band) + b_new;
			diff = x - x_old;
			x_old = x;
			ic ++;
		}
	}
	return boost::make_tuple(x, ic);
}

boost::tuple <MatrixXd, int> SOR_iteration_banded(MatrixXd A, VectorXd b, VectorXd initial, 
										   double w, int band, bool rezidual_based=1, double tol=pow(10, -6)){
	VectorXd x, r, b_new;
	MatrixXd L, D, U, M, N;
	int ic = 0;
	D = A.diagonal().asDiagonal();
	L = A.triangularView<Lower>();
	L = L - D;
	U = A.triangularView<Upper>();
	U = U - D;
	M = D + w*L;
	N = (1-w)*D - w*U;
	b_new = w * forward_subst_banded(M, b, band);
	if(rezidual_based == 1){
		x = initial;
		r = b - A * x;
		double stop_iter_resid = tol * r.norm();
		while(r.norm() > stop_iter_resid){
			x = forward_subst_banded(M, N*x, band) + b_new;
			r = b - A * x;
			ic ++;
		}
	}
	else{
		VectorXd x_old = initial;
		VectorXd diff(initial.rows());
		diff.setOnes();
		while(diff.norm()>tol){
			x = forward_subst_banded(M, N*x_old, band) + b_new;
			diff = x - x_old;
			x_old = x;
			ic ++;
		}
	}
	return boost::make_tuple(x, ic);
}



