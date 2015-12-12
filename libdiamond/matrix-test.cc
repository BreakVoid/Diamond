#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include "lu-decom.hpp"
#include "optimization.hpp"

#include <iostream>

using namespace std;
using namespace Diamond;

Vector<long double> GenGradient(const Vector<long double> &x)
{
	Vector<long double> result(2);
	result[0] = x[0];
	result[1] = 5.0 * x[1];
	return result;
}

long double F(const Vector<long double> &x)
{
	if (x.Size() != 2) {
		throw std::invalid_argument("the number of arguments is incorrect.");
	}
	return 0.5 * x[0] * x[0] + 2.5 * x[1] * x[1];
}

int main(int argc, char const *argv[])
{
	cout << "Conjugate Gradient Method test" << endl;
	Vector<long double> x0(2);
	x0[0] = 5;
	x0[1] = 1;
	cout << Optimization::ConjugateGradientMethod<long double>(2, F, GenGradient, x0);
	return 0;
}
