#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom-double.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Matrix<double> A(3, 3, 1);
	A[0][2] = 0;
	A[2][0] = 0;
	auto qrRes = QR_double::QR_DecompositionPivoting(A);
	cout << qrRes.first.first << endl;
	cout << qrRes.first.second << endl;

	cout << A << endl;
	
	cout << qrRes.first.first * qrRes.first.second * Transpose(qrRes.second) << endl;
}