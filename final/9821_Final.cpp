//============================================================================
// Name        : 9821_Final.cpp
// Author      : wz
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
//#include "Finite_Difference.hpp"
#include "Barrier_European.hpp"
using namespace std;

int main() {

	/*
	//European and American options
	IOFormat precision(16);
	cout<<setprecision(16);

	double N_initial,alpha_initial;          //initial N and alpha

	MatrixXd u;
	option op(52,50,11.0/12,0.3,0.015,0.02);

	N_initial=calculate_N(op,4,0.5);          // initial N and alpha
	alpha_initial=calculate_alpha(op,4,0.5);   // initial N and alpha

	double error;
	double v;
	for(int i=1;i<=1;i++){                    // in the for loop, every time M*4, then N*2 to keep the alpha constant
		int M=pow(4,i);
		//cout<<M<<endl;
		int N=pow(2,i-1);
		//cout<<N*N_initial<<endl;
		FD fd(op,M,N*N_initial,alpha_initial,"American Call");    // here it is!!!! alpha is constant here
		//cout << fd.alpha << "\t" << fd.N << "\t" << fd.x_left << "\t" << fd.x_right << "\t" << fd.tor_final << "\t" << endl;
		u=fd.CN_SOR_American();
		cout<<u.format(precision)<<endl;
		//v=fd.fd_value_American(u);
		//cout<<v<<endl;
		//fd.print_American(u);
		//fd.print_Greeks_American(u);
		//fd.early_exercise_domain(u);
		//error=fd.pointwise_error1_American(u,4.08381705117638);
		//cout<<error<<endl;
	}
	*/

	//cout<<u.format(precision)<<endl;

	//FD fd(op,4,0.45,"American Put");

	//u=fd.CN_SOR_American();
	//cout<<u.format(precision)<<endl;
	//fd.print_American(u);
    //fd.print_Greeks_American(u);

	//fd.early_exercise_domain(u);


	/*
	//variance reduction
	IOFormat precision(16);
	cout<<setprecision(16);
	MatrixXd uE;
	MatrixXd uA;
	option op(52,50,11.0/12,0.3,0.015,0.02);

	double N_initial=calculate_N(op,4,4);          // initial N and alpha
	double alpha_initial=calculate_alpha(op,4,4);   // initial N and alpha

	for (int i=1; i<=4; i++) {
		int M = pow(4,i);
		int N = pow(2,i-1);
		FD fdE(op,M,N*N_initial,4,"European Call");
		FD fdA(op,M,N*N_initial,4,"American Call");

		uE=fdE.CN_SOR_European();
		uA=fdA.CN_SOR_American();

		double v_red = fdA.variance_reduction_American(uE,uA);
		//cout << v_red << endl;
	}
	*/



	
    // barrier option
	IOFormat precision(16);
	double N_initial,alpha_initial;

	MatrixXd u;
	barrier_option op(50,48,3.0/12,0.3,0.005,0.02,35);
	N_initial=calculate_N(op,4,0.5,"Exam");          // initial N and alpha
	alpha_initial=calculate_alpha(op,4,0.5,"Exam");   // initial N and alpha

	for(int i=1;i<=1;i++){
		int M=pow(4,i);
		int N=pow(2,i-1);
		barrier_FD fd(op,M,N*N_initial,alpha_initial,"Exam","Call");
		u=fd.Forward_Euler_European();
		//cout<<u.format(precision)<<endl;
		//cout<<endl;
		//cout<<"M= "<<M<<endl;
	    fd.print_European(u);
	    fd.print_Greeks_European(u);
	}
	


}
