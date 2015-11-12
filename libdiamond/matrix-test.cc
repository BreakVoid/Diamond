#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Matrix<double> A(GenerateRandomMatrix<double>(4, 4, -10, 10));

	std::vector<std::pair<Matrix<double>, Matrix<double>>> qrRes;
	std::vector<Matrix<double>> arrA;
	arrA.push_back(A);
	qrRes.push_back(QR::QR_Decomposition(A));
	for (int i = 1; i <= 1000; ++i) {
		qrRes.push_back(QR::QR_Decomposition(arrA[i - 1]));
		arrA.push_back(qrRes[i].second * qrRes[i].first);
		cout << arrA.back() << endl;
	}
}