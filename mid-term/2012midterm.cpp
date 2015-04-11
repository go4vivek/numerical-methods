#include "linear_solvers.h"
#include <iomanip>
//#include "selfio.h"
#include "portfolio.h"
#include "eigen.h"
#include <fstream> 

IOFormat CSVFormat(FullPrecision, 0, "\t", "\n");

int main(int argc, char const *argv[])
{
	MatrixXd data(35,9);
	MatrixXd U, X, beta, b, con;
	ifstream ifile("data.csv");
	ifile >> data;
	//cout << data << endl;
	price_mat re(data);
	//cout << re.re << endl;
	//cout << re.cov << endl;
	//U = cholesky(re.cov);
	//cout << U <<endl;

	b = re.re.col(0);
	X.resize(b.rows(), 8);
	//con.resize(b.rows(), 1);
	//con.setOnes();
	X << re.re.rightCols(8);
	beta = least_squares(X, b);
	//cout << beta << endl;
	cout << std::fixed << std::setprecision(50) << (X*beta - b).norm()<< endl;
	ofstream ofile("out.csv");
    ofile << (beta).format(CSVFormat);
	return 0;
}