#include "../../libdiamond/polynomial.hpp"

#include <iostream>

using namespace std;
using namespace Diamond;

int main(int argc, char const *argv[])
{
	Polynomial<double> p1({1.0, 1.0});
	Polynomial<double> p2({-1.0, 1.0});
	auto p3 = pow(p1, 20);
	cout << "--- Test pow ---" << endl;
	cout << p3 << endl;
	cout << "--- Test derivate ---" << endl;
	cout << Derivate(p3) << endl;
	cout << Derivate(p3, 2) << endl;
	cout << Derivate(p3, 4) << endl;
	return 0;
}