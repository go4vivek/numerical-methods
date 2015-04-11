#include <iostream>
#include <Eigen/Core>
#include "lu_decomposition.hpp"
#include <fstream> 
#include "eigen.h"
#include "cholesky_decomposition.hpp"

using namespace std;
using namespace Eigen;

IOFormat CSVFormat(11, 0, ", ", "\n");
IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

int main()
{
    int n = 7;
    // input
    MatrixXd A(n,n);
    ifstream ifile("in.csv");
    ifile >> A;

    //cout << "The matrix A is:\n" << A.format(HeavyFmt) << endl;
    
    VectorXd b(n);
    ifstream ifile2("in2.csv");
    ifile2 >> b;

    // 1
    MatrixXd U;
    U = cholesky(A);
    //cout << "U:\r" << U.format(HeavyFmt) << endl;
    //cout << U << endl;

    // 2
    VectorXd x = linear_solve_cholesky(A, b);
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
    ofstream ofile("out.csv");
    ofile << x.format(CSVFormat);
    return 0;
}
