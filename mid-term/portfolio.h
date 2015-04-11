#include <Eigen/Dense>
#include <string>

using namespace std;
using namespace Eigen;

class price_mat{
public:
	MatrixXd price;
	MatrixXd re;
	MatrixXd cov;
	price_mat(MatrixXd _price, string method);

};

price_mat::price_mat(MatrixXd _price, string method="log"){
	VectorXd v;
	price = _price;
	int rows = _price.rows() - 1;
	int cols = _price.cols();
	re.resize(rows, cols);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if(method == "log") re(i, j) = log(_price(i, j) / _price(i+1, j));
			else re(i, j) = _price(i, j) / _price(i+1, j) - 1;
		}
	}
	cov.resize(cols, cols);
	MatrixXd mean = re.colwise().mean();
	//cout << mean << endl;
	v.resize(rows);
	v.setOnes();
	cov = (re - v * mean).transpose() * (re - v * mean) / (rows - 1);
}