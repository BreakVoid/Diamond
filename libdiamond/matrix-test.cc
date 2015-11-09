#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom-double.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	double data[4][4] = {
		{2.9766, 0.3945, 0.4198, 1.1159},
		{0.3945, 2.7328, -0.3097, 0.1129},
		{0.4198, -0.3097, 2.5675, 0.6079},
		{1.1159, 0.1129, 0.6079, 1.7231}
	};
	Matrix<double> A(4, 4);
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			A[i][j] = data[i][j];
		}
	}

	std::vector<std::pair<std::pair<Matrix<double>, Matrix<double>>, Matrix<double>>> qrRes;
	std::vector<Matrix<double>> arrA;
	arrA.push_back(A);
	qrRes.push_back(QR_double::QR_DecompositionPivoting(A));
	for (int i = 1; i <= 10; ++i) {
		qrRes.push_back(QR_double::QR_DecompositionPivoting(arrA[i - 1]));
		arrA.push_back(qrRes[i].first.second * (Transpose(qrRes[i].second) * qrRes[i].first.first));
		cout << arrA.back() << endl;
	}
}