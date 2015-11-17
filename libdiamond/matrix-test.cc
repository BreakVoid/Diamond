#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include "lu-decom.hpp"

#include <iostream>

using namespace std;
using namespace Diamond;

Vector<double> GenerateFx(const Vector<double> &x)
{
	Vector<double> res(2);
	res[0] = x[0] + 2 * x[1] - 2;
	res[1] = x[0] * x[0] + 4 * x[1] * x[1] - 4;
	return res;
}

Matrix<double> GenerateJfx(const Vector<double> &x)
{
	Matrix<double> res(2, 2);
	res[0][0] = 1; res[0][1] = 2;
	res[1][0] = 2 * x[0]; res[1][1] = 8 * x[1];
	return res;
}

int main(int argc, char const *argv[])
{
	cout << "Solve Nonlinear equation system:" << endl;
	cout << 
		"	 x_1    +  2x_2    - 2 = 0\n"
		"	(x_1)^2 + 4(x_2)^2 - 4 = 0\n";
	Vector<double> x(2);
	x[0] = 1; x[1] = 2;
	while (true) {
		Vector<double> f = GenerateFx(x);
		Matrix<double> J = GenerateJfx(x);
		Vector<double> s = LU::SolveLinearEquationSystem(J, f * double(-1));
		auto newX = x + s;
		if (EqualZero(newX[0] - x[0]) && EqualZero(newX[1] - x[1])) {
			x = newX;
			break;
		}
		x = newX;
	}
	cout << "Final result:" << x << endl;
	return 0;
}