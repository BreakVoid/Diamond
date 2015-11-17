#include "matrix.hpp"
#include "vector.hpp"
#include "qr-decom.hpp"
#include "lu-decom.hpp"
#include <iostream>

using namespace std;
using namespace Diamond;

int main()
{
	Matrix<double> A(GenerateRandomMatrix<double>(4, 4, -10, 10));
	cout << A << endl;
	auto luRes = LU::LU_Decomposition(A);
	cout << "L = " << luRes.first << endl;
	cout << "U = " << luRes.second << endl;
	cout << "L * U = " << luRes.first * luRes.second << endl;
	auto luPRes = LU::LU_DecompositionPivoting(A);
	cout << "L = " << luPRes.first.first << endl;
	cout << "U = " << luPRes.first.second << endl;
	cout << "P = " << luPRes.second << endl;
	cout << "LUP' = " << luPRes.first.first * luPRes.first.second * Transpose(luPRes.second) << endl;
}