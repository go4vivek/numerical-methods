//============================================================================
// Name        : 9821_final_tree.cpp
// Author      : wz
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Binomial_Tree.hpp"
using namespace std;

int main() {
	cout<<setprecision(16);
	
	Option op(1,0.3,0.03,40, 41,0.01,"Put","European");
	double v,delta,gamma,theta;

	//cout << "\t" << op.Black_Scholes() << "\t" << op.delta() << "\t" << op.gamma() << "\t" << op.theta() << endl;
	/*
	boost::tie(v,delta,gamma,theta)=op.Average_Binomial_Tree(10000);
	cout<<"\t"<<v<<"\t"<<delta<<"\t"<<gamma<<"\t"<<theta<<endl; //<< "\t" << op.delta() << "\t" << op.gamma() << "\t" << op.theta() << endl;
	*/

	
	int N[]={10,20,40,80,160,320,640,1280};
	int M[]={11,21,41,81,161,321,641,1281};
	int R[]={20,40,80,160,320,640,1280};
	for(int i=0;i<8;i++){
	    boost::tie(v,delta,gamma,theta)=op.Trinomial_Tree(N[i]);
	    cout << "\t" << v<< "\t" << delta << "\t" << gamma << "\t" << theta << endl;
	}

	/*
	Barrier op(50, 48, 45, 0.02,0.3, 0.01, 0.6666666667, "Call", "dao");
	double v;
	for(int i=10;i<=1000;i++){
		v=op.Barrier_Trinomial_Tree(i);
		cout<<v<<endl;
	}
	*/


}
