/*
 * Solver.hpp
 *
 *  Created on: Dec 6, 2014
 *      Author: wz
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include <cmath>


using namespace std;
using namespace Eigen;

VectorXd forward_subst(MatrixXd L,VectorXd b)
{
	int n=sqrt(L.size());
	VectorXd x(n);
    x(0)=b(0)/L(0,0);

    for(int j=1;j<n;j++)
    {
    	double sum=0.0;
    	for(int k=0;k<=j-1;k++)
    	{
    		 sum=sum+L(j,k)*x(k);
    	}
        x(j)=(b(j)-sum)/L(j,j);
    }
    return x;
}
VectorXd backward_subst(MatrixXd U,VectorXd b)
{
	int n=sqrt(U.size());
	VectorXd x(n);
	x(n-1)=b(n-1)/U(n-1,n-1);
	for(int j=(n-2);j>=0;j--)
	{
		double sum=0.0;
		for(int k=j+1;k<=n-1;k++)
		{
			sum=sum+U(j,k)*x(k);
		}
		x(j)=(b(j)-sum)/U(j,j);
	}
    return x;
}

void lu_no_pivoting(MatrixXd& L,MatrixXd& U,MatrixXd A)
{
	int n=sqrt(A.size());
	for(int i=0;i<=n-2;i++)
	{
		for(int k=i;k<=n-1;k++)
		{
			U(i,k)=A(i,k);
			L(k,i)=A(k,i)/U(i,i);
		}
		for(int j=i+1;j<=n-1;j++)
		{
			for(int k=i+1;k<=n-1;k++)
			{
				A(j,k)=A(j,k)-L(j,i)*U(i,k);
			}
		}
	}
	L(n-1,n-1)=1;
	U(n-1,n-1)=A(n-1,n-1);
}

VectorXd mult_inv_D(MatrixXd D,VectorXd x)
{
	return D.inverse()*x;
}
VectorXd mult_M_A(MatrixXd M,VectorXd x)
{
	return M*x;
}

VectorXd SOR_Iteration(MatrixXd A,VectorXd b,VectorXd x0,double tol,double w,int flag)
{
	IOFormat precision(12);
	if(w<=0||w>=2)
	{
		std::cout<<"The Omega should be in the right range.";
		exit(1);
	}
	int n=sqrt(A.size());
	VectorXd x(n);
	x=x0;
    VectorXd r0(n);
    r0=b-A*x0;
	VectorXd r(n);
	r=r0;
	double stop_iter_resid=tol*r0.norm();


	MatrixXd D(n,n);
	MatrixXd L(n,n);
	MatrixXd U(n,n);
	D.setZero();
	L.setZero();
	U.setZero();
	for(int i=0;i<n;i++)
	{
		D(i,i)=A(i,i);
	}
	for(int i=1;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			L(i,j)=A(i,j);
		}
	}
	for(int j=1;j<n;j++)
	{
		for(int i=0;i<j;i++)
		{
		    U(i,j)=A(i,j);
		}
	}

	VectorXd b_new(n);
	b_new=w*forward_subst(D+w*L,b);
	int ic=0;
	if(flag==0)
	{
		while(r.norm()>stop_iter_resid)
		{
			x=forward_subst(D+w*L,(1-w)*mult_M_A(D,x)-w*mult_M_A(U,x))+b_new;
			r=b-A*x;
			ic+=1;
		    if(ic<=3)
		    {
		    	cout<<"x="<<endl<<x.format(precision)<<endl;
		    }
		}
	}
	else if(flag==1)
	{
		VectorXd x_old(n);
		x_old=x0;
		VectorXd diff(n);
		diff.setOnes();
		while(diff.norm()>tol)
		{
			x=forward_subst(D+w*L,(1-w)*mult_M_A(D,x_old)-w*mult_M_A(U,x_old))+b_new;
			diff=x-x_old;
			x_old=x;
			ic+=1;
		}
	}
	return x;

}

VectorXd SOR_Iteration_American(MatrixXd A,VectorXd b,VectorXd x0,double tol,double w,int flag)
{
	IOFormat precision(12);
	if(w<=0||w>=2)
	{
		std::cout<<"The Omega should be in the right range.";
		exit(1);
	}
	int n=sqrt(A.size());
	VectorXd x(n);
	x=x0;

	if(flag==1)
	{
		VectorXd x_old(n);
		x_old=x0;
		VectorXd diff(n);
		diff.setOnes();
		while(diff.norm()>tol)
		{
			for(int j=0;j<n;j++){
				double sum1=0.0;
				double sum2=0.0;
				for(int k=0;k<=j-1;k++){
					sum1+=A(j,k)*x(k);
				}
				for(int k=j+1;k<n;k++){
					sum2+=A(j,k)*x_old(k);
				}
				x(j)=fmax(x0(j),(1-w)*x_old(j)-w/A(j,j)*(sum1+sum2)+w*b(j)/A(j,j));
			}
			diff=x-x_old;
			x_old=x;
		}
	}
	return x;

}





#endif /* SOLVER_HPP_ */
