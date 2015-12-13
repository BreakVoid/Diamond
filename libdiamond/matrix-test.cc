#include <iostream>
#include "matrix.hpp"
#include "util.hpp"
#include "vector.hpp"
#include "optimization.hpp"
#include "lu-decom.hpp"

using namespace std;
using namespace Diamond;

template<typename _Td>
Matrix<_Td> GeneratePositiveDefiniteMatrix(const size_t &n, const _Td &minValue = static_cast<_Td>(0), const _Td &maxValue = static_cast<_Td>(1))
{
	Matrix<_Td> res(n, n);
	std::uniform_real_distribution<_Td> random(minValue, maxValue);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			res[i][j] = random(engine);
		}
	}
	return Transpose(res) * res;
}

const size_t nSize = 5;

int main(int argc, char const *argv[])
{
	Matrix<long double> A = GeneratePositiveDefiniteMatrix(nSize, -5.0L, 5.0L);
	Vector<long double> b = GenerateRandomVector(nSize, -50.0L, 50.0L);

	cout << LU::SolveLinearEquationSystem(A, b) << endl;
	cout << Optimization::ConjugateGradientMethod<long double>(nSize, [&A, &b](const Vector<long double> &x) -> long double {
		return 0.5L * Transpose(x) * A * x - Transpose(x) * b;
	}, [&A, &b](const Vector<long double> &x) -> Vector<long double> {
		return A * x - b;
	}) << endl;
	return 0;
}