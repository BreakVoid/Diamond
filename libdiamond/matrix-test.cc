#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Matrix<double> A(GenerateRandomSymmetricMatrix<double>(4, -10, 10));
	std::vector<std::pair<std::pair<Matrix<double>, Matrix<double>>, Matrix<double>>> qrRes;
	std::vector<Matrix<double>> arrA;
	arrA.push_back(A);
	cout << A << endl;
	qrRes.push_back(QR::QR_DecompositionPivoting(A));
	for (int i = 1; i <= 1000; ++i) {
		qrRes.push_back(QR::QR_DecompositionPivoting(arrA[i - 1]));
		arrA.push_back(qrRes[i].first.second * Transpose(qrRes[i].second) * qrRes[i].first.first);
	}
	cout << arrA.back() << endl;
}