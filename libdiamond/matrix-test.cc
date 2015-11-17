#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom.hpp"
#include "lu-decom.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Matrix<double> a(2, 2);
	a[0][0] = 1; a[0][1] = 2;
	a[1][0] = -1.67; a[1][1] = 11.3;
	Vector<double> b(2);
	b[0] = 0; b[1] = -4.72;
	const auto s = LU::SolveLinearEquationSystem(a, b);
	cout << "Solution :" << s << endl;
}