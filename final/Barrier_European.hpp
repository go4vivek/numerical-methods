/*
 * Boundary_European.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: wz
 */

#ifndef BARRIER_EUROPEAN_HPP_
#define BARRIER_EUROPEAN_HPP_

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

class barrier_option{
public:
	double s;
	double k;
	double T;
	double sigma;
	double q;
	double r;
	double B;
    barrier_option(double _s,double _k,double _T,double _sigma,double _q,double _r,double _B):s(_s),k(_k),T(_T),sigma(_sigma),q(_q),r(_r),B(_B) {}
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

double vbarriercall(barrier_option op){
	double barrier_strike=op.B*op.B/op.s;
	option op1(op.s,op.k,op.T,op.sigma,op.q,op.r);
	option op2(barrier_strike,op.k,op.T,op.sigma,op.q,op.r);
	double aa=(op.r-op.q)/(op.sigma*op.sigma)-0.5;
	double value=vcall(op1)-pow(op.B/op.s,2*aa)*vcall(op2);
	return value;
}

double vput(option op)
{
	double d1 = (log(op.s/op.k)+(op.r-op.q+op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double d2 = (log(op.s/op.k)+(op.r-op.q-op.sigma*op.sigma/2)*op.T)/(op.sigma*sqrt(op.T));
	double value = op.k*exp(-op.r*op.T)*std_normal_cdf(-d2)-op.s*exp(-op.q*op.T)*std_normal_cdf(-d1);
	return value;
}



class barrier_FD{
public:
	barrier_option op;
	double x_left;
	double x_right;
	double tor_final;
	double a;
	double b;
	int M;
	int N;
	double alpha_temp;
	double alpha;

	string barrier_type;   //"Down and Out" "Up and Out"
	string type;   // "Call" "Put"

	double f(double x);         // boundary condition
	double g_left(double x);
	double g_right(double x);

	barrier_FD(barrier_option,int m,int n,double alpha,string barrier_type,string type_);
	barrier_FD(barrier_option,int m,double alpha_temp_,string barrier_type_,string type_);
	double pointwise_error1_European(MatrixXd u);
	double pointwise_error2_European(MatrixXd u);
	void print_European(MatrixXd u);
	void print_Greeks_European(MatrixXd u);

	MatrixXd Forward_Euler_European();
	MatrixXd Backward_Euler_LU_European();
	MatrixXd Backward_Euler_SOR_European();
	MatrixXd CN_LU_European();
	MatrixXd CN_SOR_European();
};

double barrier_FD::f(double x){
	if(type=="Call"){
		return op.k*exp(a*x)*fmax(exp(x)-1,0);
	}
}

double barrier_FD::g_left(double x){
	if(type=="Call"){
		return 0;
	}
}

double barrier_FD::g_right(double x){
	if(type=="Call"){
		return op.k*exp(a*x_right+b*x)*(exp(x_right-2*op.q*x/pow(op.sigma,2))-exp(-2*op.r*x/pow(op.sigma,2)));
	}
}

double calculate_N(barrier_option op,int M,double alpha_temp,string barrier_type){
	double x_left,x_right;
	if(barrier_type=="Down and Out"){
		x_left=log(op.B/op.k);
		x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	}
	if(barrier_type=="Up and Out"){
		x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
		x_right=log(op.B/op.k);
	}
	if(barrier_type=="Exam") {
		x_left = log(40.0/op.k);
		x_right = log(60.0/op.k);
	}
	double tor_final=op.T*op.sigma*op.sigma/2.0;
	//double a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	//double b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	double delta_tor=tor_final/M;
	int N=floor((x_right-x_left)/sqrt(delta_tor/alpha_temp));
	cout<<"initial N is "<<N<<endl;
	return N;
}

double calculate_alpha(barrier_option op,int M,double alpha_temp,string barrier_type){
	double x_left,x_right;
	if(barrier_type=="Down and Out"){
		x_left=log(op.B/op.k);
		x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	}
	if(barrier_type=="Up and Out"){
		x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
		x_right=log(op.B/op.k);
	}
	if(barrier_type=="Exam") {
		x_left = log(40.0/op.k);
		x_right = log(60.0/op.k);
	}
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




barrier_FD::barrier_FD(barrier_option option,int m,int n,double alpha_,string barrier_type_,string type_):op(option),M(m),N(n),alpha(alpha_),barrier_type(barrier_type_),type(type_){
	if(barrier_type=="Down and Out"){
		x_left=log(op.B/op.k);
		x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	}
	if(barrier_type=="Up and Out"){
		x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
		x_right=log(op.B/op.k);
	}
	if(barrier_type=="Exam") {
		x_left = log(40.0/op.k);
		x_right = log(60.0/op.k);
	}
	tor_final=op.T*op.sigma*op.sigma/2.0;
	a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);
	alpha_temp=0;

}

barrier_FD::barrier_FD(barrier_option option,int M_,double alpha_temp_,string barrier_type_,string type_):op(option),M(M_),alpha_temp(alpha_temp_),barrier_type(barrier_type_),type(type_){
	if(barrier_type=="Down and Out"){
		x_left=log(op.B/op.k);
		x_right=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T+3*op.sigma*sqrt(op.T);
	}
	if(barrier_type=="Up and Out"){
		x_left=log(op.s/op.k)+(op.r-op.q-0.5*op.sigma*op.sigma)*op.T-3*op.sigma*sqrt(op.T);
		x_right=log(op.B/op.k);
	}
	if(barrier_type=="Exam") {
		x_left = log(40.0/op.k);
		x_right = log(60.0/op.k);
	}
	tor_final=op.T*op.sigma*op.sigma/2.0;
	a=(op.r-op.q)/pow(op.sigma,2)-0.5;
	b=pow((op.r-op.q)/pow(op.sigma,2)+0.5,2)+2*op.q/pow(op.sigma,2);

	double delta_tor=tor_final/M;
	N=floor((x_right-x_left)/sqrt(delta_tor/alpha_temp));
	double delta_x=(x_right-x_left)/N;
	alpha=delta_tor/(delta_x*delta_x);

}

MatrixXd barrier_FD::Forward_Euler_European(){
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

MatrixXd barrier_FD::Backward_Euler_LU_European(){
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

MatrixXd barrier_FD::Backward_Euler_SOR_European(){
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
			x0(j)=u_temp(j);
		}
	}
	return u;
}

MatrixXd barrier_FD::CN_LU_European(){
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

MatrixXd barrier_FD::CN_SOR_European(){
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
		u_temp(i) = u(0,i+1);
		x0(i)=u(0,i+1);
	}

	for (int i = 1; i <= M; i++)
	{
		b = B * u_temp;
		b(0) += alpha/2*(u(i,0)+u(i-1,0)); b(N-2) += alpha/2*(u(i,N)+u(i-1,N));
		//u_temp = SOR_Iteration(A,b,x0,1.2,0.000001,count);
		u_temp=SOR_Iteration(A,b,x0,0.000001,1.2,1);
		for (int j = 0; j <= N-2; j++) {
			u(i,j+1) = u_temp(j) ;
			x0(j)=u_temp(j);
		}
	}
	return u;
}
void barrier_FD::print_European(MatrixXd u){
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
	//cout << "hehe\t" << alpha << "\t" << N << "\t" << x_left << "\t" << x_right << "\t" << x_compute << "\t" << tor_final << endl;
	//cout<<"U(M,i) "<<u(M,i);
	//cout<<endl;
	//cout<<"U(M,i+1) "<<u(M,i+1);
	//cout<<endl;
	double si=op.k*exp(xi);
    //cout<<"S[i] "<<si<<endl;
	double si1=op.k*exp(xi1);
	//cout<<"S[i+1] "<<si1<<endl;

	//for (int k=0; k<=N; k++) {
	//	cout << exp(-a*xi-b*tor_final)*u(M,k) << "\t";
	//}
	double vi=exp(-a*xi-b*tor_final)*u(M,i);
	//cout<<"V[i] "<<vi<<endl;
	double vi1=exp(-a*xi1-b*tor_final)*u(M,i+1);
	//cout<<"V[i+1] "<<vi1<<endl;

	double V_app1=((si1-op.s)*vi+(op.s-si)*vi1)/(si1-si);
	//cout<<"the approximate value of option(type I) "<<V_app1<<endl;

	double u_inter=((xi1-x_compute)*u(M,i)+(x_compute-xi)*u(M,i+1))/(xi1-xi);
	double V_app2=exp(-a*x_compute-b*tor_final)*u_inter;
	//cout<<"the approximate value of option(type II) "<<V_app2<<endl;

	double V_exact=vbarriercall(op);         // barrier call
	double error1_European=abs(V_app1-V_exact);
	//cout<<"error_pointwise(type I) "<<error1_European<<endl;
	//cout<<endl;

	cout << u(M,i) << "\t" << u(M,i+1) << "\t" << V_app1 << "\t";
}

double barrier_FD::pointwise_error1_European(MatrixXd u){
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

	double V_exact=vbarriercall(op);         // barrier call
	double error1_European=abs(V_app1-V_exact);
	double error2_European=abs(V_app2-V_exact);

	return error1_European;
}

double barrier_FD::pointwise_error2_European(MatrixXd u){
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

	double V_exact=vbarriercall(op);      // barrier call
	double error1_European=abs(V_app1-V_exact);
	double error2_European=abs(V_app2-V_exact);

	return error2_European;
}



void barrier_FD::print_Greeks_European(MatrixXd u){
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
	//cout<<"delta(European) "<<delta<<endl;
	double gamma=((Vi_plus2-Vi_plus1)/(si_plus2-si_plus1)-(Vi-Vi_minus1)/(si-si_minus1))/((si_plus2+si_plus1)*0.5-(si+si_minus1)*0.5);
	//cout<<"gamma(European) "<<gamma<<endl;
	double delta_tor=tor_final/M;
	double delta_t=2*delta_tor/(op.sigma*op.sigma);

	double Vi_delta_t=exp(-a*xi-b*(tor_final-delta_tor))*u(M-1,i);
	double Vi_plus1_delta_t=exp(-a*xi_plus1-b*(tor_final-delta_tor))*u(M-1,i+1);

	double V_app_delta_t=((si_plus1-op.s)*Vi_delta_t+(op.s-si)*Vi_plus1_delta_t)/(si_plus1-si);
	double V_app_0=((si_plus1-op.s)*Vi+(op.s-si)*Vi_plus1)/(si_plus1-si);

	double theta=-(V_app_0-V_app_delta_t)/delta_t;   // the negative sign here. be careful
	//cout<<"theta(European) "<<theta<<endl;
	//cout<<endl;

	cout << delta << "\t" << gamma << "\t" << theta << endl;
}









#endif /* BARRIER_EUROPEAN_HPP_ */
