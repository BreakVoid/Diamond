#include <iostream>
#include "matrix.hpp"
#include "util.hpp"
#include "vector.hpp"
#include "optimization.hpp"
#include "lu-decom.hpp"
#include <sstream>
#include <fstream>

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

int savedSolutionCnt = 0;

void SaveSolution(const Matrix<long double> &A, const Vector<long double> &b, const Vector<long double> &LU_Sol, const Vector<long double> &CGM_Sol)
{
	ostringstream oss;
	oss << setw(3) << setfill('0') << savedSolutionCnt;
	string filename = "matrix-and-solution-" + oss.str();
	ofstream ofs(filename);
	ofs << "The positive definite matrix is " << A << "\n";
	ofs << "The vector b is " << b << "\n";
	ofs << "The solution by using LU decomposition is " << LU_Sol << "\n";
	ofs << "The solution by using Conjugate Gradient Methos is " << CGM_Sol << "\n";
	ofs.close();
	savedSolutionCnt++;
}

int main(int argc, char const *argv[])
{
	while (true) {
		Matrix<long double> A = GeneratePositiveDefiniteMatrix(nSize, -5.0L, 5.0L);
		Vector<long double> b = GenerateRandomVector(nSize, -50.0L, 50.0L);
		auto x = LU::SolveLinearEquationSystem(A, b);
		auto y = Optimization::ConjugateGradientMethod<long double>(nSize, [&A, &b](const Vector<long double> &x) -> long double {
			return 0.5L * Transpose(x) * A * x - Transpose(x) * b;
		}, [&A, &b](const Vector<long double> &x) -> Vector<long double> {
			return A * x - b;
		});
		cout << x << endl;
		cout << y << endl;
		cout << "Do you want to save this case?(y/n): ";
		string s;
		cin >> s;
		if (s == "y") {
			SaveSolution(A, b, x, y);
		}
	}
	return 0;
}