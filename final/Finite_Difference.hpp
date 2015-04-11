/*
 * Finite_Difference.hpp
 *
 *  Created on: Dec 6, 2014
 *      Author: wz
 */

#ifndef FINITE_DIFFERENCE_HPP_
#define FINITE_DIFFERENCE_HPP_
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include "Solver.hpp"
#include <boost/math/distributions/normal.hpp> // for normal_distribution
const double PI=3.14159265358979323846;
using namespace std;
using namespace Eigen;

class option{
public:
	double s;
	double k;
	double T;
	double sigma;
	double q;
	double r;

	option(double _s,double _k,double _T,double _sigma,double _q,double _r):s(_s),k(_k),T(_T),sigma(_sigma),q(_q),r(_r) {}
};

double std_normal_cdf(double t)
{
  	boost::math::normal s;
  	return cdf(s,t);
}

double vcall(option op)
{
	double d1 = (log(op.s/op.k)+(op.r-op.q+op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double d2 = (log(op.s/op.k)+(op.r-op.q-op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double value = op.s*exp(-op.q*op.T)*std_normal_cdf(d1)-op.k*exp(-op.r*op.T)*std_normal_cdf(d2);
	return value;
}

double vput(option op)
{
	double d1 = (log(op.s/op.k)+(op.r-op.q+op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double d2 = (log(op.s/op.k)+(op.r-op.q-op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double value = op.k*exp(-op.r*op.T)*std_normal_cdf(-d2)-op.s*exp(-op.q*op.T)*std_normal_cdf(-d1);
	return value;
}



class FD{
public:
	option op;
	double x_left;
	double x_right;
	double tor_final;
	double a;
	double b;
	int M;
	int N;

	double alpha_temp;
	double alpha;

	string type;

	double f(double x);         // boundary condition
	double g_left(double x);
	double g_right(double x);

	FD(option option,int m,double alpha_temp_,string type_);
	FD(option option,int m,int n,double alpha_,string type_);
	double pointwise_error1_European(MatrixXd u);
	double pointwise_error2_European(MatrixXd u);
	void print_European(MatrixXd u);
	double pointwise_error1_American(MatrixXd u,double V_exact);
	double pointwise_error2_American(MatrixXd u,double V_exact);
	void print_American(MatrixXd u);
	double RMS_error_European(MatrixXd u);
	void print_Greeks_European(MatrixXd u);
	void print_Greeks_American(MatrixXd u);

	double fd_value_American(MatrixXd u);
	double fd_value_European(MatrixXd u);

	double variance_reduction_American(MatrixXd uE,MatrixXd uA);

	MatrixXd Forward_Euler_European();
	MatrixXd Backward_Euler_LU_European();
	MatrixXd Backward_Euler_SOR_European();
	MatrixXd CN_LU_European();
	MatrixXd CN_SOR_European();
	MatrixXd Forward_Euler_American();
	MatrixXd CN_SOR_American();
	void early_exercise_domain(MatrixXd u);
};

double FD::f(double x){
	if(type=="European Put"){
		return op.k*exp(a*x)*fmax(1.0-exp(x),0);
	}
	if(type=="American Put"){
		return op.k*exp(a*x)*fmax(1.0-exp(x),0);
	}
	if(type=="American Call"){
		return op.k*exp(a*x)*fmax(exp(x)-1.0,0);
	}
	if(type=="European Call")
		return op.k*exp(a*x)*fmax(exp(x)-1.0,0);
}

double FD::g_left(double x){
	if(type=="European Put"){
		return op.k*exp(a*x_left+b*x)*(exp(-2*op.r*x/pow(op.sigma,2))-exp(x_left-2*op.q*x/pow(op.sigma,2)));
	}
	if(type=="American Put"){
		return op.k*exp(a*x_left+b*x)*(1-exp(x_left));
	}
	if(type=="American Call"){
		return 0;
	}
	if(type=="European Call")
		return 0;
}

double FD::g_right(double x){
	if(type=="European Put"){
		return 0;
	}
	if(type=="American Put"){
		return 0;
	}
	if(type=="American Call"){
		return op.k*exp(a*x_right+b*x)*(exp(x_right)-1);
	}
	if(type=="European Call")
		return op.k*exp(a*x_left+b*x)*(exp(x_left-2*op.q*x/pow(op.sigma,2))-exp(-2*op.r*x/pow(op.sigma,2)));
}

double calculate_N(option op,int M,double alpha_temp){
	double x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
	double x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	double tor_final=op.T*op.sigma*op.sigma/2.0;
	//double a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	//double b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	double delta_tor=tor_final/M;
	int N=floor((x_right-x_left)/sqrt(delta_tor/alpha_temp));
	cout<<"initial N is "<<N<<endl;
	return N;
}



double calculate_alpha(option op,int M,double alpha_temp){
	double x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
	double x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	double tor_final=op.T*op.sigma*op.sigma/2.0;
	//double a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	//double b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	double delta_tor=tor_final/M;
	int N=floor((x_right-x_left)/sqrt(delta_tor/alpha_temp));
	double delta_x=(x_right-x_left)/N;
	double alpha=delta_tor/(delta_x*delta_x);
	cout<<"initial alpha is "<<alpha<<endl;
	return alpha;
}
FD::FD(option option,int m,int n,double alpha_,string type_):op(option),M(m),N(n),alpha(alpha_),type(type_){
	x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
	x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	tor_final=op.T*op.sigma*op.sigma/2.0;
	a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	alpha_temp=0;

}
FD::FD(option option,int M_,double alpha_temp_,string type_):op(option),M(M_),alpha_temp(alpha_temp_),type(type_){
	x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
	x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	tor_final=op.T*op.sigma*op.sigma/2.0;
	a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	double delta_tor=tor_final/M;
	double delta_x;

	N=floor((x_right-x_left)/sqrt(delta_tor/alpha_temp));
	delta_x=(x_right-x_left)/N;
	alpha=delta_tor/(delta_x*delta_x);

}



MatrixXd FD::Forward_Euler_European(){
	// matrix used for recursive calculation
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	VectorXd b(N-1);
	VectorXd u_temp(N-1);
	u.setZero();
	A.setZero();
	b.setZero();
	u_temp.setZero();

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;
	// assignment according to boundary conditions
	// discrete x axis and assign values to Matrix U according to boundary conditions(t=0)
	double xtemp = 0;
	for (int i = 0; i <= N; i++){
		xtemp = x_left + i * dx;
		u(0,i) = f(xtemp); }

	// discrete t axis and assign values to Matrix U according to boundary conditions(x=0 and x=N)
	double t_temp = 0;
	for (int i = 0; i <= M; i++){
		t_temp = i * dt;
		u(i,0) = g_left(t_temp);
		u(i,N) = g_right(t_temp); }

	// initialize matrix A
	for (int i = 0; i < N-1; i++){
		A(i,i) = 1-2*alpha; }
	for (int i = 0; i < N-2; i++){
		A(i,i+1) = alpha;
		A(i+1,i) = alpha; }

	// initialize matrix u_temp, which is used for recursive calculation.
	// be careful about the range of the u_temp vector
	for (int i = 0; i <= N-2; i++){
		u_temp(i) = u(0,i+1); }

	// forward Euler method, begin from i=1,for i=0, it is the boundary condition
	for (int i = 1; i <= M; i++)
	{
		b(0) = alpha*u(i-1,0);
		b(N-2) = alpha*u(i-1,N);
		u_temp = A * u_temp + b;
		for (int j = 0; j <= N-2; j++)
		{ u(i,j+1) = u_temp(j) ; }
	}
	return u;
}

MatrixXd FD::Backward_Euler_LU_European(){
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	MatrixXd L(N-1,N-1);
	MatrixXd U(N-1,N-1);
	VectorXd b(N-1);
	VectorXd y(N-1);
	VectorXd u_temp(N-1);
	u.setZero();
	A.setZero();
	L.setZero();
	U.setZero();
	b.setZero();
	y.setZero();
	u_temp.setZero();

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;

	double xtemp = 0;
	for (int i = 0; i <= N; i++) { xtemp = x_left + i * dx; u(0,i) = f(xtemp); }

	double t_temp = 0;
	for (int i = 0; i <= M; i++) { t_temp = i * dt; u(i,0) = g_left(t_temp); u(i,N) = g_right(t_temp); }

	for (int i = 0; i < N-1; i++) { A(i,i) = 1+2*alpha; }
	for (int i = 0; i < N-2; i++) { A(i,i+1) = -alpha; A(i+1,i) = -alpha; }
	lu_no_pivoting(L,U,A);

	for (int i = 0; i <= N-2; i++) { u_temp(i) = u(0,i+1); }

	for (int i = 1; i <= M; i++)
	{
		b(0) = alpha*u(i,0); b(N-2) = alpha*u(i,N);
		y = forward_subst(L,(u_temp+b));
		u_temp = backward_subst(U,y);
		for (int j = 0; j <= N-2; j++) { u(i,j+1) = u_temp(j) ; }
	}
	return u;
}

MatrixXd FD::Backward_Euler_SOR_European(){
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	VectorXd b(N-1);
	VectorXd u_temp(N-1);
	VectorXd x0(N-1);
	u.setZero();
	A.setZero();
	b.setZero();
	u_temp.setZero();
	//x0.setOnes();
	int count = 0;

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;

	double xtemp = 0;
	for (int i = 0; i <= N; i++) { xtemp = x_left + i * dx; u(0,i) = f(xtemp); }

	double t_temp = 0;
	for (int i = 0; i <= M; i++) { t_temp = i * dt; u(i,0) = g_left(t_temp); u(i,N) = g_right(t_temp); }

	for (int i = 0; i < N-1; i++) { A(i,i) = 1+2*alpha; }
	for (int i = 0; i < N-2; i++) { A(i,i+1) = -alpha; A(i+1,i) = -alpha; }

	for (int i = 0; i <= N-2; i++) {
		u_temp(i) = u(0,i+1);
		x0(i)=u(0,i+1);
	}

	for (int i = 1; i <= M; i++)
	{
		b(0) = alpha*u(i,0); b(N-2) = alpha*u(i,N);
		u_temp=SOR_Iteration(A,u_temp+b,x0,0.000001,1.2,1);
		for (int j = 0; j <= N-2; j++) {
			u(i,j+1) = u_temp(j) ;
			x0(j) = u_temp(j);
		}
	}
	return u;
}

MatrixXd FD::CN_LU_European(){
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	MatrixXd L(N-1,N-1);
	MatrixXd U(N-1,N-1);
	MatrixXd B(N-1,N-1);
	VectorXd b(N-1);
	VectorXd y(N-1);
	VectorXd u_temp(N-1);
	u.setZero();
	A.setZero();
	L.setZero();
	U.setZero();
	B.setZero();
	b.setZero();
	y.setZero();
	u_temp.setZero();

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;

	double xtemp = 0;
	for (int i = 0; i <= N; i++) { xtemp = x_left + i * dx; u(0,i) = f(xtemp); }

	double t_temp = 0;
	for (int i = 0; i <= M; i++) { t_temp = i * dt; u(i,0) = g_left(t_temp); u(i,N) = g_right(t_temp); }

	for (int i = 0; i < N-1; i++) { A(i,i) = 1+alpha; }
	for (int i = 0; i < N-2; i++) { A(i,i+1) = -alpha/2; A(i+1,i) = -alpha/2; }
	lu_no_pivoting(L,U,A);

	for (int i = 0; i < N-1; i++) { B(i,i) = 1-alpha; }
	for (int i = 0; i < N-2; i++) { B(i,i+1) = alpha/2; B(i+1,i) = alpha/2; }

	for (int i = 0; i <= N-2; i++) { u_temp(i) = u(0,i+1); }

	for (int i = 1; i <= M; i++)
	{
		b = B * u_temp;
		b(0) += alpha/2*(u(i,0)+u(i-1,0));
		b(N-2) += alpha/2*(u(i,N)+u(i-1,N));
		y = forward_subst(L,b);
		u_temp = backward_subst(U,y);
		for (int j = 0; j <= N-2; j++) { u(i,j+1) = u_temp(j) ; }
	}
	return u;
}

MatrixXd FD::CN_SOR_European(){
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	MatrixXd B(N-1,N-1);
	VectorXd b(N-1);
	VectorXd u_temp(N-1);
	VectorXd x0(N-1);
	u.setZero();
	A.setZero();
	B.setZero();
	b.setZero();
	u_temp.setZero();
	//x0.setOnes();
	int count = 0;

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;

	double xtemp = 0;
	for (int i = 0; i <= N; i++) { xtemp = x_left + i * dx; u(0,i) = f(xtemp); }

	double t_temp = 0;
	for (int i = 0; i <= M; i++) { t_temp = i * dt; u(i,0) = g_left(t_temp); u(i,N) = g_right(t_temp); }

	for (int i = 0; i < N-1; i++) { A(i,i) = 1+alpha; }
	for (int i = 0; i < N-2; i++) { A(i,i+1) = -alpha/2; A(i+1,i) = -alpha/2; }

	for (int i = 0; i < N-1; i++) { B(i,i) = 1-alpha; }
	for (int i = 0; i < N-2; i++) { B(i,i+1) = alpha/2; B(i+1,i) = alpha/2; }

	for (int i = 0; i <= N-2; i++) {
		u_temp(i) = u(0,i+1);        // boundary condition
		x0(i) = u(0,i+1);            // first initial guess
	}



	for (int i = 1; i <= M; i++)
	{
		b = B * u_temp;
		b(0) += alpha/2*(u(i,0)+u(i-1,0)); b(N-2) += alpha/2*(u(i,N)+u(i-1,N));
		u_temp=SOR_Iteration(A,b,x0,0.00000001,1.2,1);
		for (int j = 0; j <= N-2; j++) {
			u(i,j+1) = u_temp(j) ;
			x0(j)=u_temp(j);        // update the initial guess
		}
	}
	return u;
}
void FD::print_European(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
    cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	cout<<"S[i+1]"<<si1<<endl;

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	cout<<"the approximate value of option(type II) "<<V_app2<<endl;
	cout<<endl;

}

double FD::pointwise_error1_European(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	//cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	//cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	//cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	//cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
	//cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	//cout<<"S[i+1]"<<si1<<endl;

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	//cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	//cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	//cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	//cout<<"the approximate value of option(type II) "<<V_app2<<endl;

	double V_exact=vput(op);         // put option
	double error1_European=abs(V_app1-V_exact);
	double error2_European=abs(V_app2-V_exact);

	return error1_European;
}

double FD::pointwise_error2_European(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	//cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	//cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	//cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	//cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
	//cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	//cout<<"S[i+1]"<<si1<<endl;

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	//cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	//cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	//cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	//cout<<"the approximate value of option(type II) "<<V_app2<<endl;

	double V_exact=vput(op);        // put option
	double error1_European=abs(V_app1-V_exact);
	double error2_European=abs(V_app2-V_exact);

	return error2_European;
}

double FD::RMS_error_European(MatrixXd u){
	double delta_x=(x_right-x_left)/N;
	double s0=op.s;
	int count=0;
	double sum=0.0;

	for(int k=0;k<=N;k++){
		double xk=x_left+k*delta_x;
		double sk=op.k*exp(xk);
		double V_app=exp(-a*xk-b*tor_final)*u(M,k);
		op.s=sk;
		double V_exact=vput(op);        // put option
		if(V_exact>0.00001*s0){
			sum+=(V_app-V_exact)*(V_app-V_exact)/(V_exact*V_exact);
			count++;
		}
	}
	double error=sqrt(sum/count);
	return error;
}

void FD::print_Greeks_European(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	double xi=x_left+i*delta_x;
	double xi_minus1=x_left+(i-1)*delta_x;
	double xi_plus1=x_left+(i+1)*delta_x;
	double xi_plus2=x_left+(i+2)*delta_x;

	double si=op.k*exp(xi);
	double si_minus1=op.k*exp(xi_minus1);
	double si_plus1=op.k*exp(xi_plus1);
	double si_plus2=op.k*exp(xi_plus2);

	double Vi=exp(-a*xi-b*tor_final)*u(M,i);
	double Vi_minus1=exp(-a*xi_minus1-b*tor_final)*u(M,i-1);
	double Vi_plus1=exp(-a*xi_plus1-b*tor_final)*u(M,i+1);
	double Vi_plus2=exp(-a*xi_plus2-b*tor_final)*u(M,i+2);

	double delta=(Vi_plus1-Vi)/(si_plus1-si);
	cout<<"delta(European) "<<delta<<endl;
	double gamma=((Vi_plus2-Vi_plus1)/(si_plus2-si_plus1)-(Vi-Vi_minus1)/(si-si_minus1))/((si_plus2+si_plus1)*0.5-(si+si_minus1)*0.5);
	cout<<"gamma(European) "<<gamma<<endl;
	double delta_tor=tor_final/M;
	double delta_t=2*delta_tor/(op.sigma*op.sigma);

	double Vi_delta_t=exp(-a*xi-b*(tor_final-delta_tor))*u(M-1,i);
	double Vi_plus1_delta_t=exp(-a*xi_plus1-b*(tor_final-delta_tor))*u(M-1,i+1);

	double V_app_delta_t=((si_plus1-op.s)*Vi_delta_t+(op.s-si)*Vi_plus1_delta_t)/(si_plus1-si);
	double V_app_0=((si_plus1-op.s)*Vi+(op.s-si)*Vi_plus1)/(si_plus1-si);

	double theta=-(V_app_0-V_app_delta_t)/delta_t;   // the negative sign here. be careful
	cout<<"theta(European) "<<theta<<endl;
	cout<<endl;
}

MatrixXd FD::Forward_Euler_American(){
	// matrix used for recursive calculation
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	VectorXd B(N-1);
	VectorXd u_temp(N-1);
	MatrixXd early_ex_premium(M,N-1);   // for American options
	u.setZero();
	A.setZero();
	B.setZero();
	u_temp.setZero();

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;
	// assignment according to boundary conditions
	// discrete x axis and assign values to Matrix U according to boundary conditions(t=0)
	double xtemp = 0;
	for (int i = 0; i <= N; i++){
		xtemp = x_left + i * dx;
		u(0,i) = f(xtemp); }

	// discrete t axis and assign values to Matrix U according to boundary conditions(x=0 and x=N)
	double t_temp = 0;
	for (int i = 0; i <= M; i++){
		t_temp = i * dt;
		u(i,0) = g_left(t_temp);
		u(i,N) = g_right(t_temp); }

	// initialize matrix A
	for (int i = 0; i < N-1; i++){
		A(i,i) = 1-2*alpha; }
	for (int i = 0; i < N-2; i++){
		A(i,i+1) = alpha;
		A(i+1,i) = alpha; }

	for(int m=0;m<M;m++){
		for(int i=0;i<N-1;i++){
			double x=x_left+(i+1)*dx;
			double tor=(m+1)*dt;
			early_ex_premium(m,i)=op.k*exp(a*x+b*tor)*fmax(1-exp(x),0);
		}
	}

	// initialize matrix u_temp, which is used for recursive calculation.
	// be careful about the range of the u_temp vector
	for (int i = 0; i <= N-2; i++){
		u_temp(i) = u(0,i+1); }

	// forward Euler method, begin from i=1,for i=0, it is the boundary condition
	for (int i = 1; i <= M; i++)
	{
		B(0) = alpha*u(i-1,0);
		B(N-2) = alpha*u(i-1,N);
		u_temp = A * u_temp + B;
		for (int j = 0; j <= N-2; j++)
		{ u(i,j+1) = fmax(u_temp(j),early_ex_premium(i-1,j));
		  u_temp(j) = u(i,j+1);}

	}

	return u;
}
MatrixXd FD::CN_SOR_American(){
	MatrixXd u(M+1,N+1);
	MatrixXd A(N-1,N-1);
	MatrixXd B(N-1,N-1);
	VectorXd _b(N-1);
	VectorXd u_temp(N-1);
	VectorXd x0(N-1);
	MatrixXd early_ex_premium(M,N-1); //for American option
	u.setZero();
	A.setZero();
	B.setZero();
	_b.setZero();
	u_temp.setZero();
	//x0.setOnes();
	int count = 0;

	double dx = (x_right-x_left)/N;
	double dt = tor_final/M;

	double xtemp = 0;
	for (int i = 0; i <= N; i++) { xtemp = x_left + i * dx; u(0,i) = f(xtemp); }

	double t_temp = 0;
	for (int i = 0; i <= M; i++) { t_temp = i * dt; u(i,0) = g_left(t_temp); u(i,N) = g_right(t_temp); }

	for (int i = 0; i < N-1; i++) { A(i,i) = 1+alpha; }
	for (int i = 0; i < N-2; i++) { A(i,i+1) = -alpha/2; A(i+1,i) = -alpha/2; }

	for (int i = 0; i < N-1; i++) { B(i,i) = 1-alpha; }
	for (int i = 0; i < N-2; i++) { B(i,i+1) = alpha/2; B(i+1,i) = alpha/2; }

	for (int i = 0; i <= N-2; i++) { u_temp(i) = u(0,i+1); }
	// American option

	for(int m=0;m<M;m++){
		for(int i=0;i<N-1;i++){
			double x=x_left+(i+1)*dx;
			double tor=(m+1)*dt;
			early_ex_premium(m,i)=op.k*exp(a*x+b*tor)*fmax(1-exp(x),0);
		}
	}

	for (int i = 1; i <= M; i++)
	{
		_b = B * u_temp;
		_b(0) += alpha/2*(u(i,0)+u(i-1,0)); _b(N-2) += alpha/2*(u(i,N)+u(i-1,N));
		x0=early_ex_premium.row(i-1);
		u_temp=SOR_Iteration_American(A,_b,x0,0.00000001,1.2,1);
		for (int j = 0; j <= N-2; j++) {
			u(i,j+1)=u_temp(j);
		}
	}
	return u;
}

double FD::pointwise_error1_American(MatrixXd u,double V_exact){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	//cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	//cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	//cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	//cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
	//cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	//cout<<"S[i+1]"<<si1<<endl;

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	//cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	//cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	//cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	//cout<<"the approximate value of option(type II) "<<V_app2<<endl;

	double error1=abs(V_app1-V_exact);
	return error1;
}
double FD::pointwise_error2_American(MatrixXd u,double V_exact){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	//cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	//cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	//cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	//cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
	//cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	//cout<<"S[i+1]"<<si1<<endl;

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	//cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	//cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	//cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	//cout<<"the approximate value of option(type II) "<<V_app2<<endl;

	double error2=abs(V_app2-V_exact);
	return error2;
}

void FD::print_American(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	cout<<"x_compute "<<x_compute<<endl;
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	cout<<"i "<<i<<endl;
	double xi=x_left+i*delta_x;
	cout<<"x[i] "<<xi<<endl;
	double xi1=x_left+(i+1)*delta_x;
	cout<<"x[i+1] "<<xi1<<endl;

	double si=op.k*exp(xi);
	cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	cout<<"S[i+1]"<<si1<<endl;

	cout << "haha" << u(M,i) << "\t" << u(M,i+1) << endl;
	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	cout<<"the approximate value of option(type II) "<<V_app2<<endl;
}

void FD::print_Greeks_American(MatrixXd u){
	std::cout<<setprecision(16);
	double x_compute=log(op.s/op.k);
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	double xi=x_left+i*delta_x;
	double xi_minus1=x_left+(i-1)*delta_x;
	double xi_plus1=x_left+(i+1)*delta_x;
	double xi_plus2=x_left+(i+2)*delta_x;

	double si=op.k*exp(xi);
	double si_minus1=op.k*exp(xi_minus1);
	double si_plus1=op.k*exp(xi_plus1);
	double si_plus2=op.k*exp(xi_plus2);

	double Vi=exp(-a*xi-b*tor_final)*u(M,i);
	double Vi_minus1=exp(-a*xi_minus1-b*tor_final)*u(M,i-1);
	double Vi_plus1=exp(-a*xi_plus1-b*tor_final)*u(M,i+1);
	double Vi_plus2=exp(-a*xi_plus2-b*tor_final)*u(M,i+2);

	double delta=(Vi_plus1-Vi)/(si_plus1-si);
	cout<<"delta(American) "<<delta<<endl;
	double gamma=((Vi_plus2-Vi_plus1)/(si_plus2-si_plus1)-(Vi-Vi_minus1)/(si-si_minus1))/((si_plus2+si_plus1)*0.5-(si+si_minus1)*0.5);
	cout<<"gamma(American) "<<gamma<<endl;
	double delta_tor=tor_final/M;
	double delta_t=2*delta_tor/(op.sigma*op.sigma);

	double Vi_delta_t=exp(-a*xi-b*(tor_final-delta_tor))*u(M-1,i);
	double Vi_plus1_delta_t=exp(-a*xi_plus1-b*(tor_final-delta_tor))*u(M-1,i+1);

	double V_app_delta_t=((si_plus1-op.s)*Vi_delta_t+(op.s-si)*Vi_plus1_delta_t)/(si_plus1-si);
	double V_app_0=((si_plus1-op.s)*Vi+(op.s-si)*Vi_plus1)/(si_plus1-si);

	double theta=-(V_app_0-V_app_delta_t)/delta_t;   // the negative sign here. be careful
	cout<<"theta(American) "<<theta<<endl;
}

void FD::early_exercise_domain(MatrixXd u){

	double delta_x = (x_right-x_left)/N;
	double delta_tor = tor_final/M;
	MatrixXd early_ex_premium(M,N-1);
	vector<int> mm;
	vector<int> nn;
	vector<double> s;
	vector<double> t;
	for(int m=0;m<M;m++){
		for(int i=0;i<N-1;i++){
			double x=x_left+(i+1)*delta_x;
			double tor=(m+1)*delta_tor;
			early_ex_premium(m,i)=op.k*exp(a*x+b*tor)*fmax(1-exp(x),0);
		}
	}

	for(int m=1;m<=M;m++){
		for(int n=1;n<=N-2;n++){
			if((u(m,n)==early_ex_premium(m-1,n-1))&&(u(m,n+1)>early_ex_premium(m-1,n))){
				mm.push_back(m);
				nn.push_back(n);
				break;
			}
		}
	}

	for(int i=0;i<nn.size();i++){
		double xn=x_left+nn[i]*delta_x;
		double xn1=x_left+(nn[i]+1)*delta_x;
		double sn=op.k*exp(xn);
		double sn1=op.k*exp(xn1);
		s.push_back(0.5*(sn+sn1));
		double tt=op.T-2*mm[i]*delta_tor/(op.sigma*op.sigma);
		t.push_back(tt);
	}

	cout<<"early exercise domain:"<<endl;
	for(int i=0;i<s.size();i++){
		cout<<setprecision(16);
		cout<<t[i]<<"       "<<s[i]<<endl;
	}


}

double FD::fd_value_American(MatrixXd u){
	double x_compute=log(op.s/op.k);
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	double xi=x_left+i*delta_x;
	double xi1=x_left+(i+1)*delta_x;

	double si=op.k*exp(xi);
	double si1=op.k*exp(xi1);

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);

	return V_app1;
}

double FD::fd_value_European(MatrixXd u){
	double x_compute=log(op.s/op.k);
	double delta_x=(x_right-x_left)/N;
	int i=floor((x_compute-x_left)/delta_x);
	double xi=x_left+i*delta_x;
	double xi1=x_left+(i+1)*delta_x;

	double si=op.k*exp(xi);
	double si1=op.k*exp(xi1);

	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	return V_app1;
}

double FD::variance_reduction_American(MatrixXd uE,MatrixXd uA){
	double v_bs;
	double v_red;
	if(type=="American Call") v_bs=vcall(op);
	if(type=="American Put") v_bs=vput(op);
	v_red=fd_value_American(uA)+(v_bs-fd_value_European(uE));
	cout << fd_value_American(uA) << "\t" << fd_value_European(uE) << "\t" << v_red << endl;
	return v_red;
}

void implied_vol(){
	double sigma_old=0.1;
	double sigma_new=0.4;
	double sigma_oldest;

	double tol=0.0001;
	MatrixXd u_new;
	MatrixXd u_old;
	double value_new;
	double value_old;
	double value_oldest;
	double T=0.083333333333333;

	while(fabs(sigma_new-sigma_old)>tol){
		option op_old(38,40,T,sigma_old,0.01,0.04);  // change of sigma here
		FD fd1(op_old,16,0.4,"American Put");
		u_old=fd1.Forward_Euler_American();
		value_old=fd1.fd_value_American(u_old);

		option op_new(38,40,T,sigma_new,0.01,0.04);  // change of sigma here
		FD fd2(op_new,16,0.4,"American Put");
		u_new=fd2.Forward_Euler_American();
		value_new=fd2.fd_value_American(u_new);

		sigma_oldest=sigma_old;
		value_oldest=value_old;
		sigma_old=sigma_new;
		value_old=value_new;
		sigma_new=sigma_old-(value_old-2.45)*(sigma_old-sigma_oldest)/(value_old-value_oldest);
		cout<<sigma_new<<endl;
	}
}




#endif /* FINITE_DIFFERENCE_HPP_ */
