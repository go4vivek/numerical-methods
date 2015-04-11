#ifndef LINEAR_SOLVER_HPP_
#define LINEAR_SOLVER_HPP_

#include <Eigen/Dense>
#include <boost/tuple/tuple.hpp>
#include "lu_decomposition.hpp"

using namespace Eigen;

VectorXd mult_inv_D(MatrixXd D,VectorXd x)
{
	return D.inverse()*x;
}
VectorXd mult_M_A(MatrixXd M,VectorXd x)
{
	return M*x;
}

boost::tuple<VectorXd,int> Jocobi_Iteration(MatrixXd A,VectorXd b,VectorXd x0,double tol,int flag)
{
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

	int ic=0;
	VectorXd b_new(n);
	b_new=mult_inv_D(D,b);
	// residual-based criterion

	if(flag==0)
	{
	    while(r.norm()>stop_iter_resid)
	    {
	    	//std::cout << x;
		    x=-mult_inv_D(D,L*x+U*x)+b_new;
		    r=b-A*x;
		    ic+=1;
	    }
	}

	// consecutive approximation criterion
	else if(flag==1)
	{
		VectorXd x_old(n);
		x_old=x0;
		VectorXd diff(n);
		diff.setOnes();

		while(diff.norm()>tol)
		{
			x=-mult_inv_D(D,L*x_old+U*x_old)+b_new;
			diff=x-x_old;
			x_old=x;
			ic+=1;
		}
	}

	return boost::make_tuple(x,ic);
}


boost::tuple<VectorXd,int> Gauss_Siedel_Iteration(MatrixXd A,VectorXd b,VectorXd x0,double tol,int flag)
{
	int n=sqrt(A.size());
	VectorXd x(n);
	x=x0;
	VectorXd r0(n);
	r0=b-A*x0;
	VectorXd r(n);
	r=r0;
	double stop_iter_resid=tol*r0.norm();


	MatrixXd M(n,n);
	MatrixXd N(n,n);
	M.setZero();
	N.setZero();
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<=i;j++)
		{
			M(i,j)=A(i,j);
		}
	}
	for(int j=1;j<n;j++)
	{
		for(int i=0;i<j;i++)
		{
			N(i,j)=A(i,j);
		}
	}

	VectorXd b_new(n);
	b_new=forward_subst(M,b);
	int ic=0;

	//  residual-based criterion
	if(flag==0)
	{
		while(r.norm()>stop_iter_resid)
		{
			//std::cout<<x;
			x=-forward_subst(M,mult_M_A(N,x))+b_new;
			r=b-A*x;
			ic+=1;
		}
	}
	// consecutive approximation criterion
	else if(flag==1)
	{
		VectorXd x_old(n);
		x_old=x0;
		VectorXd diff(n);
		diff.setOnes();
		while(diff.norm()>tol)
		{
			x=-forward_subst(M,mult_M_A(N,x_old))+b_new;
			diff=x-x_old;
			x_old=x;
			ic+=1;
		}
	}


	return boost::make_tuple(x,ic);
}


boost::tuple<VectorXd,int> SOR_Iteration(MatrixXd A,VectorXd b,VectorXd x0,double tol,double w,int flag)
{
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
			//cout << x;
			x=forward_subst(D+w*L,(1-w)*mult_M_A(D,x)-w*mult_M_A(U,x))+b_new;
			r=b-A*x;
			ic+=1;
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


	return boost::make_tuple(x,ic);

}

#endif /* LINEAR_SOLVER_HPP_ */