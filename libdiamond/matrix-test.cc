#include <iostream>
#include "matrix.hpp"
#include "util.hpp"
#include "vector.hpp"
#include "optimization.hpp"
#include "lu-decom.hpp"

using namespace std;
using namespace Diamond;

int main(int argc, char const *argv[])
{
	Vector<double> a = GenerateRandomVector<double>(5);
	Vector<double> b = GenerateRandomVector<double>(5);
	cout << a << endl;
	cout << b << endl;
	cout << a + b << endl;
	VectorT<double> aT = GenerateRandomVector<double>(5);
	VectorT<double> bT = GenerateRandomVector<double>(5);
	cout << aT << endl;
	cout << bT << endl;
	cout << aT + bT << endl;

	cout << a << endl;
	cout << Transpose(a) << endl;
	cout << Transpose(Transpose(a)) << endl;
	cout << Transpose(a) * a << endl;
	cout << a * Transpose(a) << endl;
	return 0;
}