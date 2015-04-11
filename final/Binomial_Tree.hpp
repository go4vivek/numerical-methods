/*
 * Binomial_Tree.hpp
 *
 *  Created on: Dec 15, 2014
 *      Author: wz
 */

#ifndef BINOMIAL_TREE_HPP_
#define BINOMIAL_TREE_HPP_
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/math/distributions/normal.hpp>
const double PI=3.14159265358979323846;


using namespace std;


class Option{
public:
	double T;
	double sigma;
	double r;
	int K;
	double q;
	double S;
	string type;        // Put or Call
	string style;       // European or American
	Option(double _T, double _sigma, double _r, int _K, double _S, double _q, string _type, string _style="European");
	double payoff(double spot);
	double Black_Scholes(double spot, double t);
	double Black_Scholes();
	double delta();
	double gamma();
	double theta();
	boost::tuple <double, double, double ,double> Binomial_Tree(int N);
	boost::tuple <double, double, double ,double> Average_Binomial_Tree(int N);
	boost::tuple <double, double, double ,double> Binomial_BS(int N);
	boost::tuple <double, double, double ,double> BBSR(int N);
	boost::tuple <double, double, double ,double> Trinomial_Tree(int N);
	boost::tuple <double, double, double ,double> Trinomial_BS(int N);
	boost::tuple <double, double, double ,double> TBSR(int N);
};

class Barrier: public Option{
public:
	double B;
	string direction;
	Barrier(double _S, double _K, double _B, double _r, double _sigma, double _q, double _T, string _type, string _direction);
	double Barrier_Binomial_Tree(int N);
	double Barrier_Trinomial_Tree(int N);
};

Barrier::Barrier(double _S, double _K, double _B, double _r, double _sigma, double _q, double _T, string _type, string _direction):
		Option(_T, _sigma, _r, _K, _S, _q, _type)
{
	B = _B;
	direction = _direction;
}



double Barrier::Barrier_Binomial_Tree(int N){

	double v10, v11, s10, s11, v20, v21, v22, s20, s21, s22;
	double v, delta, gamma, theta;
	double * V = new double[N+1];
	double det_T, u, d, p_u, p_d;
	det_T = T / N;
	u = exp(sigma * sqrt(det_T));
	d = 1 / u;
	p_u = (exp((r-q)*det_T) - d) / (u-d);
	p_d = 1 - p_u;

	s10=S*u;
	s11=S*d;
	s20=S*u*u;
	s21=S;
	s22=S*d*d;

	for (int i = 0; i <= N; ++i)
	{

		V[i] = payoff(S * pow(u, N-i) * pow(d, i));
		if(direction=="dao" && S * pow(u, N-i) * pow(d, i)<=B) V[i] = 0;
		else if(direction=="uao" && S * pow(u, N-i) * pow(d, i)>=B) V[i] = 0;
	}
	for (int i = N-1; i >=0; --i)
	{
		for (int j=0; j<=i; ++j)
				{
					V[j] = exp(-det_T*r)*(p_u*V[j]+p_d*V[j+1]);
					if(direction=="dao" && S * pow(u, i-j) * pow(d, j)<=B) V[j] = 0;
					else if(direction=="uao" && S * pow(u, i-j) * pow(d, j)>=B) V[j] = 0;
					//if (i == N-1) cout << V[j] << endl;
				}

		//if (i == 1) {v10 = V[0]; v11 = V[1];}
		//if (i == 2) {v20 = V[0]; v21 = V[1]; v22 = V[2];}
	}
	v = V[0];
	//delta = (v10 - v11) / (s10 - s11);
	//gamma = ((v20 - v21) / (s20 - s21) - (v21 - v22) / (s21 - s22)) / ((s20 - s22) / 2);
	//theta = (v21 - v) / (2 * det_T);
	delete []V;
	//return boost::make_tuple(v, delta, gamma, theta);
	return v;

}

double Barrier::Barrier_Trinomial_Tree(int N){
	double v10, v11, v12, s10, s12, v20, v24, v22, s20, s24, s22;
	double v, delta, gamma, theta;
	double * V = new double[2*N+1];
	double det_T, u, d, p_u, p_d, p_m;
	det_T = T / N;
	u = exp(sigma * sqrt(3 * det_T));
	d = 1 / u;
	p_u = double(1) / 6 + (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_d = double(1) / 6 - (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_m = double(2) / 3;
	s10 = S * u;
	s12 = S * d;
	s20 = S * u * u;
	s22 = S;
	s24 = S * d * d;
	for (int i = 0; i <= 2 * N; ++i)
	{
		V[i] = payoff(S * pow(u, N-i));
		if(direction=="dao" && S * pow(u, N-i)<=B) V[i] = 0;
		else if(direction=="uao" && S * pow(u, N-i)>=B) V[i] = 0;
	}
	for (int i = N-1; i >=0; --i)
	{
		for (int j=0; j<= 2*i; ++j)
				{
					V[j] = exp(-det_T*r)*(p_u*V[j]+p_m*V[j+1]+p_d*V[j+2]);
					if(direction=="dao" && S * pow(u, i-j)<=B) V[j] = 0;
					else if(direction=="uao" && S * pow(u, i-j)>=B) V[j] = 0;
					//if (i == N-1) cout << V[j] << endl;
				}
		//if (i == 1) {v10 = V[0]; v11 = V[1]; v12 = V[2];}
		//if (i == 2) {v20 = V[0]; v22 = V[2]; v24 = V[4];}	/* code */
	}
	v = V[0];
	//delta = (v10 - v12) / (s10 - s12);
	//gamma = ((v20 - v22) / (s20 - s22) - (v22 - v24) / (s22 - s24)) / (s10 - s12);
	//theta = (v11 - v) / det_T;
	delete []V;
	//return boost::make_tuple(v, delta, gamma, theta);
	return v;
}

Option::Option(double _T, double _sigma, double _r, int _K, double _S, double _q, string _type, string _style){
	T = _T;
	sigma = _sigma;
	r = _r;
	K = _K;
	S = _S;
	q = _q;
	type = _type;
	style = _style;
}

double Option::payoff(double spot){
	if(type == "Put") return fmax(K - spot, double(0));
	else return fmax(spot-K, 0.0);
}


double Option::Black_Scholes(double spot, double t ){
	double v, d1, d2;
	boost::math::normal norm;
	d1 = (log(spot/K) + (r-q+sigma*sigma/2)*t) / (sigma * sqrt(t));
	d2 = d1 - sigma * sqrt(t);
	if(type == "Put") v = exp(-r*t)*K*cdf(norm, -d2) - spot * exp(-q*t)*cdf(norm, -d1) ;
	else v = spot * exp(-q*t)*cdf(norm, d1) - exp(-r*t)*K*cdf(norm, d2) ;
	return v;
}

double Option::Black_Scholes(){
	double v, d1, d2;
	boost::math::normal norm;
	d1 = (log(S/K) + (r-q+sigma*sigma/2)*T) / (sigma * sqrt(T));
	d2 = d1 - sigma * sqrt(T);
	if(type == "Put") v = exp(-r*T)*K*cdf(norm, -d2) - S * exp(-q*T)*cdf(norm, -d1);
	else v = S * exp(-q*T)*cdf(norm, d1) - exp(-r*T)*K*cdf(norm, d2) ;
	return v;
}

double Option::delta(){
	boost::math::normal norm;
	if(type=="Call"){
		double d1 = (log(S/K)+(r-q+sigma*sigma/2)*T)/(sigma*sqrt(T));
		return exp(-q*T)*cdf(norm,d1);
	}
	if(type=="Put"){
		double d1 = (log(S/K)+(r-q+sigma*sigma/2)*T)/(sigma*sqrt(T));
		return -exp(-q*T)*cdf(norm,-d1);
	}
}

double Option::gamma(){
	double d1 = (log(S/K)+(r-q+sigma*sigma/2)*T)/(sigma*sqrt(T));
	return (exp(-q*T-d1*d1/2))/(S*sigma*sqrt(2*T*PI));
}

double Option::theta(){
	boost::math::normal norm;
	if(type=="Call"){
		double d1 = (log(S/K)+(r-q+sigma*sigma/2)*T)/(sigma*sqrt(T));
		double d2 = (log(S/K)+(r-q-sigma*sigma/2)*T)/(sigma*sqrt(T));
		double theta = -S*sigma*exp(-q*T-d1*d1/2)/(2*sqrt(2*T*PI))+q*S*exp(-q*T)*cdf(norm,d1)-r*K*exp(-r*T)*cdf(norm,d2);
		return theta;
	}
	if(type=="Put"){
		double d1 = (log(S/K)+(r-q+sigma*sigma/2)*T)/(sigma*sqrt(T));
		double d2 = (log(S/K)+(r-q-sigma*sigma/2)*T)/(sigma*sqrt(T));
		double theta = -S*sigma*exp(-q*T-d1*d1/2)/(2*sqrt(2*T*PI))-q*S*exp(-q*T)*cdf(norm,-d1)+r*K*exp(-r*T)*cdf(norm,-d2);
		return theta;
	}
}


boost::tuple <double, double, double ,double> Option::Binomial_Tree(int N){

	double v10, v11, s10, s11, v20, v21, v22, s20, s21, s22;
	double v, delta, gamma, theta;
	double * V = new double[N+1];
	double det_T, u, d, p_u, p_d;
	det_T = T / N;
	u = exp(sigma * sqrt(det_T));
	d = 1.0/ u;
	p_u = (exp((r-q)*det_T) - d) / (u-d);
	p_d = 1 - p_u;

	for (int i = 0; i <= N; ++i)
	{
		V[i] = payoff(S * pow(u, N-i) * pow(d, i));
	}
	for (int i = N-1; i >=0; --i)
	{
		for (int j=0; j<=i; ++j)
				{
					V[j] = exp(-det_T*r)*(p_u*V[j]+p_d*V[j+1]);
					if(style == "American") V[j] = fmax(V[j], payoff(S * pow(u, i-j) * pow(d, j)));
					//if (i == N-1) cout << V[j] << endl;
				}
		if (i == 1) {v10 = V[0]; v11 = V[1];}
		if (i == 2) {v20 = V[0]; v21 = V[1]; v22 = V[2];}
	}
	v = V[0];
	s10=S*u;
	s11=S*d;
	s20=S*u*u;
	s21=S;
	s22=S*d*d;

	delta = (v10 - v11) / (s10 - s11);
	gamma = ((v20 - v21) / (s20 - s21) - (v21 - v22) / (s21 - s22)) / ((s20 - s22) / 2);
	theta = (v21 - v) / (2 * det_T);
	delete []V;
	return boost::make_tuple(v, delta, gamma, theta);
}

boost::tuple <double, double, double ,double> Option::Trinomial_Tree(int N){
	double v10, v11, v12, s10, s12, v20, v24, v22, s20, s24, s22;
	double v, delta, gamma, theta;
	double * V = new double[2*N+1];
	double det_T, u, d, p_u, p_d, p_m;
	det_T = T / N;
	u = exp(sigma * sqrt(3 * det_T));
	d = 1 / u;
	p_u = double(1) / 6 + (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_d = double(1) / 6 - (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_m = double(2) / 3;
	s10 = S * u;
	s12 = S * d;
	s20 = S * u * u;
	s22 = S;
	s24 = S * d * d;

	for (int i = 0; i <= 2 * N; ++i)
	{
		V[i] = payoff(S * pow(u, N-i));
	}
	for (int i = N-1; i >=0; --i)
	{
		for (int j=0; j<= 2*i; ++j)
				{
					V[j] = exp(-det_T*r)*(p_u*V[j]+p_m*V[j+1]+p_d*V[j+2]);
					if(style == "American") V[j] = fmax(V[j], payoff(S * pow(u, i-j)));
					//if (i == N-1) cout << V[j] << endl;
				}
		if (i == 1) {v10 = V[0]; v11 = V[1]; v12 = V[2];}
		if (i == 2) {v20 = V[0]; v22 = V[2]; v24 = V[4];}	/* code */
	}
	v = V[0];
	delta = (v10 - v12) / (s10 - s12);
	gamma = ((v20 - v22) / (s20 - s22) - (v22 - v24) / (s22 - s24)) / (s10 - s12);
	theta = (v11 - v) / det_T;
	delete []V;
	return boost::make_tuple(v, delta, gamma, theta);
}

boost::tuple <double, double, double ,double> Option::Average_Binomial_Tree(int N){
	double v1, v2, delta1, delta2, gamma1, gamma2, theta1, theta2;
	boost::tie(v1, delta1, gamma1, theta1) = Binomial_Tree(N);
	boost::tie(v2, delta2, gamma2, theta2) = Binomial_Tree(N+1);
	double v,delta,gamma,theta;
	v=(v1+v2)/2.0;
	delta=(delta1+delta2)/2.0;
	gamma=(gamma1+gamma2)/2.0;
	theta=(theta1+theta2)/2.0;
	return boost::make_tuple(v,delta,gamma,theta);
}

boost::tuple <double, double, double ,double> Option::Binomial_BS(int N){
	double v10, v11, s10, s11, v20, v21, v22, s20, s21, s22;
	double v, delta, gamma, theta;
	double * V = new double[N];
	double det_T, u, d, p_u, p_d;
	det_T = T / N;
	u = exp(sigma * sqrt(det_T));
	d = 1 / u;
	p_u = (exp((r-q)*det_T) - d) / (u-d);
	p_d = 1 - p_u;
	s10 = S * u;
	s11 = S * d;
	s20 = S * u * u;
	s21 = S;
	s22 = S * d * d;
	for (int i = 0; i <= N-1; ++i)
	{
		V[i] = Black_Scholes(S*pow(u, N-1-i)*pow(d,i), det_T);
		if(style == "American") V[i] = fmax(payoff(S*pow(u, N-1-i)*pow(d,i)), V[i]);
		//cout << V[i] << endl; // this part needs to be done
	}
	for (int i = N-2; i >=0; --i)
	{
		for (int j=0; j<=i; ++j)
				{
					V[j] = exp(-det_T*r)*(p_u*V[j]+p_d*V[j+1]);
					if(style == "American") V[j]=fmax(V[j], payoff(S * pow(u, i-j) * pow(d, j)));
				}
		if (i == 1) {v10 = V[0]; v11 = V[1];}
		if (i == 2) {v20 = V[0]; v21 = V[1]; v22 = V[2];}	/* code */
	}
	v = V[0];
	delta = (v10 - v11) / (s10 - s11);
	gamma = ((v20 - v21) / (s20 - s21) - (v21 - v22) / (s21 - s22)) / ((s20 - s22) / 2);
	theta = (v21 - v) / (2 * det_T);
	delete []V;
	return boost::make_tuple(v, delta, gamma, theta);
}

boost::tuple <double, double, double ,double> Option::Trinomial_BS(int N){
	double v10, v11, v12, s10, s12, v20, v24, v22, s20, s24, s22;
	double v, delta, gamma, theta;
	double * V = new double[2*N+1];
	double det_T, u, d, p_u, p_d, p_m;
	det_T = T / N;
	u = exp(sigma * sqrt(3 * det_T));
	d = 1 / u;
	p_u = double(1) / 6 + (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_d = double(1) / 6 - (r-q-sigma*sigma/2) * sqrt(det_T/(12*sigma*sigma));
	p_m = double(2) / 3;
	s10 = S * u;
	s12 = S * d;
	s20 = S * u * u;
	s22 = S;
	s24 = S * d * d;
	for (int i = 0; i <= 2 * N - 2; ++i)
	{
		V[i] = Black_Scholes(S*pow(u, N-1-i), det_T);
		if(style == "American") V[i] = fmax(payoff(S*pow(u, N-1-i)), V[i]);
		//V[i] = payoff(S * pow(u, N-i));
	}
	for (int i = N-2; i >=0; --i)
	{
		for (int j=0; j<= 2*i; ++j)
				{


					V[j] = exp(-det_T*r)*(p_u*V[j]+p_m*V[j+1]+p_d*V[j+2]);
					if(style == "American") V[j] = fmax(V[j], payoff(S * pow(u, i-j)));
					//if (i == N-1) cout << V[j] << endl;
				}
		if (i == 1) {v10 = V[0]; v11 = V[1]; v12 = V[2];}
		if (i == 2) {v20 = V[0]; v22 = V[2]; v24 = V[4];}	/* code */
	}
	v = V[0];
	delta = (v10 - v12) / (s10 - s12);
	gamma = ((v20 - v22) / (s20 - s22) - (v22 - v24) / (s22 - s24)) / (s10 - s12);
	theta = (v11 - v) / det_T;
	delete []V;
	return boost::make_tuple(v, delta, gamma, theta);
}

boost::tuple <double, double, double ,double> Option::BBSR(int N){
	double v1, v2, delta1, delta2, gamma1, gamma2, theta1, theta2;
	boost::tie(v1, delta1, gamma1, theta1) = Binomial_BS(N);
	boost::tie(v2, delta2, gamma2, theta2) = Binomial_BS(N/2);
	return boost::make_tuple(2*v1-v2, 2*delta1-delta2, 2*gamma1-gamma2, 2*theta1-theta2);
}

boost::tuple <double, double, double ,double> Option::TBSR(int N){
	double v1, v2, delta1, delta2, gamma1, gamma2, theta1, theta2;
	boost::tie(v1, delta1, gamma1, theta1) = Trinomial_BS(N);
	boost::tie(v2, delta2, gamma2, theta2) = Trinomial_BS(N/2);
	return boost::make_tuple(2*v1-v2, 2*delta1-delta2, 2*gamma1-gamma2, 2*theta1-theta2);
}

#endif /* BINOMIAL_TREE_HPP_ */
