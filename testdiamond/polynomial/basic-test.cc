#include "../../libdiamond/polynomial.hpp"

#include <iostream>

using namespace std;
using namespace Diamond;

int main(int argc, char const *argv[])
{
	Polynomial<double> p1({1.0, 1.0});
	Polynomial<double> p2({-1.0, 1.0});
	cout << "--- Testing single polynomial ---" << endl;
	cout << p1 << endl;
	cout << p2 << endl;
	cout << "--- Testing addition ---" << endl;
	cout << p1 + p2 << endl;
	cout << 5 + p1 << endl;
	cout << p1 + 20 << endl;
	cout << "--- Testing substraction ---" << endl;
	cout << p1 - p2 << endl;
	cout << 12 - p1 << endl;
	cout << p2 - 16 << endl;
	cout << "--- Testing multiplication ---" << endl;
	cout << p1 * p2 << endl;
	cout << 5 * p1 << endl;
	cout << p2 * 8 << endl;
	return 0;
}