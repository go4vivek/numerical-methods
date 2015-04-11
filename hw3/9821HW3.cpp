//============================================================================
// Name        : 9821HW3.cpp
// Author      : wz
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <boost/tuple/tuple.hpp>
#include "Jacobi_Iteration.hpp"
#include "Gauss_Siedel_Iteration.hpp"
#include "SOR_Iteration.hpp"
#include "Basic_Operation.hpp"
using namespace std;
using namespace Eigen;

boost::tuple<VectorXd,int> Jocobi_Iteration(MatrixXd A,VectorXd b,VectorXd x0,double tol=,int flag);

int main()
{
	VectorXd x(14);
	int ic;
	double tol=pow(10,-6);
	MatrixXd A(14,14);
	VectorXd b(14);
	VectorXd x0(14);

	x0.setOnes();

	for(int i=0;i<14;i++)
	{
		b(i)=pow(i,2);
	}
	for(int i=0;i<14;i++)
	{
		A(i,i)=2;
	}
	for(int i=1;i<14;i++)
	{
		A(i,i-1)=-1;
	}
	for(int i=0;i<13;i++)
	{
		A(i,i+1)=-1;
	}

	//cout<<A<<endl;
	//cout<<x0<<endl;
	//cout<<b<<endl;
	//boost::tie(x,ic)=Jocobi_Iteration(A,b,x0,tol,1);
	//boost::tie(x,ic)=Gauss_Siedel_Iteration(A,b,x0,tol,1);
	boost::tie(x,ic)=SOR_Iteration(A,b,x0,tol,1.98,0);
	cout<<"x="<<endl<<x.transpose()<<endl;
	cout<<"ic="<<ic;
}


