#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "lu_decomposition.hpp"
#include <fstream> 
#include "eigen.h"
#include "cholesky_decomposition.hpp"
#include "linear_solver.hpp"
#include <iomanip>

using namespace std;
using namespace Eigen;

IOFormat CSVFormat(FullPrecision, 0, ", ", "\n");
IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

int main()
{
    /*
    MatrixXd A(8,8);
    for (int i = 0; i <= 7; ++i)
        A(i,i) = 9;
    for (int i = 0; i <= 5; ++i)
        A(i,i+2) = -2;
    for (int i = 2; i <= 7; ++i)
        A(i,i-2) = 3;
    for (int i = 0; i <= 4; ++i)
        A(i,i+3) = -1;
    for (int i = 3; i <= 7; ++i)
        A(i,i-3) = 1;
    VectorXd b(8);
    for (int i = 0; i <= 7; ++i)
        b(i) = (3*i-4)*1./(i*i+1.);
    VectorXd x0(8);
    x0.setZero();
    double tol = pow(10,-6);
    double omega1 = 0.9, omega2 = 1.15;

    VectorXd x;
    int i;
    boost::tie(x,i) = SOR_Iteration(A, b, x0, tol, omega2, 0);
    */
    
    MatrixXd A(9,9);
    A.setZero();
    for (int i = 0; i <= 8; ++i)
        A(i,i) = 3;
    for (int i = 1; i <= 8; ++i)
        A(i,i-1) = -2;
    for (int i = 2; i <= 8; ++i)
        A(i,i-2) = 4;
    for (int i = 0; i <= 6; ++i)
        A(i,i+2) = -1;
    /*
    MatrixXd C = A.transpose()+A;
    VectorXd b(9);
    for (int i = 0; i <= 8; ++i)
        b(i) = (2*i*i-5.)/(2*i+3);
    VectorXd x1 = linear_solve_LU_no_pivoting(C,b);
    cout << std::fixed << std::setprecision(50) << (b-C*x1).norm() << std::endl;
    */

    
    MatrixXd A2 = A.transpose()*A;
    MatrixXd U = cholesky(A);

    VectorXd b(9);
    for (int i = 0; i <= 8; ++i)
        b(i) = (2*i*i-5)/(2*i+3);

    VectorXd v1 = linear_solve_cholesky(A, b);
    cout << std::fixed << std::setprecision(50) << (b-A2*v1).norm();
    


    /*
    VectorXd b(9);
    for (int i = 0; i <= 8; ++i)
        b(i) = sqrt(i*i-2*i+5);
    */



    //VectorXd v1 = linear_solve_LU_no_pivoting(A,b);

    
    //int n = 8;
    // input
    //MatrixXd A(n,n);
    //ifstream ifile("in.csv");
    //ifile >> A;

    //cout << "The matrix A is:\n" << A.format(HeavyFmt) << endl;
    
    //VectorXd b(n);
    //ifstream ifile2("in2.csv");
    //ifile2 >> b;

    //VectorXd x0(n);
    //x0.setOnes();
    //double tol=1e-6;

    // 1
    //MatrixXd U;
    //U = cholesky(A);
    //cout << "U:\r" << U.format(HeavyFmt) << endl;
    //cout << U << endl;

    // 2
    //VectorXd x;
    //int i;
    //boost::tie(x,i) = Gauss_Siedel_Iteration(A, b, x0, tol, 0);
    //cout << "x:\r" << x.format(HeavyFmt) << endl;
    //cout << x << endl;
    /* 1
    MatrixXd Ab = A * b;
    */
    

    /* 2
    VectorXd x(7); 
    x = forward_subst(A, b);
    */
    
    // 3
    //VectorXd x = backward_subst(A, b);
    

    // 4
    //MatrixXd L, U;
    //boost::tie(L, U) = lu_no_pivoting(A);

    

    // 5
    //VectorXi P;
    //MatrixXd L, U;
    //boost::tie(P, L, U) = lu_row_pivoting(A);
    //cout << "P:\r" << P.format(HeavyFmt) << endl;
    //cout << "L:\r" << L.format(HeavyFmt) << endl;
    //cout << "U:\r" << U.format(HeavyFmt) << endl;

    

    // output
    //cout << i;
    //cout << (b-A*v1).norm();
    ofstream ofile("out.csv");
    ofile << (v1).format(CSVFormat);
    //VectorXd v2 = A.inverse() * b;
    //cout << (b - A*v2).norm();
    return 0;
}
