#include "matrix.hpp"
#include "vector.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Diamond::Matrix<double> A(3, 3, 0);
	A[0][0] = A[1][1] = A[2][2] = 5;
	cout << A << endl;
	Diamond::Matrix<double> B(3, 3, 0);
	cout << B << endl;
	Diamond::Matrix<double> C = (A + B) / 5;
	cout << C << endl;
	Diamond::Vector<double> v(3);
	v[0] = 1;
	v[1] = 2;
	v[2] = 3;
	cout << v << endl;
	cout << A * v << endl;
	cout << Transpose(v) * v << endl;
	cout << v * Transpose(v) << endl;
	cout << 100.0 * Transpose(v) << endl;
	cout << 100.0 * v << endl;
}