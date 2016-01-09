#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include "lu-decom.hpp"
#include <map>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <fstream>

using namespace std;
using namespace Diamond;

template<typename _Td>
inline _Td Sqr(const _Td &x)
{
	return x * x;
}

template<typename _Td>
inline _Td Cub(const _Td &x)
{
	return x * x * x;
}

vector<long double> pos;

long double t(const int &i)
{
	if (i >= 0 && i < pos.size()) {
		return pos[i];
	} else if (i < 0) {
		return (pos[1] - pos[0]) * i;
	} else {
		return (pos[pos.size() - 1] - pos[pos.size() - 2]) * (i - pos.size() + 1) + pos.back();
	}
}

long double B(int k, int i, long double x)
{
	long double t1 = t(i);
	long double t2 = t(i + 1);
	if (k == 0) {
		return (t1 <= x && x < t2);
	} else {
		long double t11 = t(i + k);
		long double t22 = t(i + k + 1);
		long double k1 = (x - t1) / (t11 - t1);
		long double k2 = (x - t2) / (t22 - t2);
		return B(k - 1, i, x) * k1 + B(k - 1, i + 1, x) * (1 - k2);
	}
}

long double B2(int k, int i, long double x)
{
	return (B(k, i, x + 1e-4) + B(k, i, x - 1e-4) - B(k, i, x) * 2) / 1e-8;
}

vector<string> filename({"3-parts", "4-parts", "5-parts", "6-parts"});

void FiniteElementMethod()
{
	for (const int &size : {3, 4, 5, 6}) {
		pos.clear();
		for (int i = 0; i < size; ++i) {
			pos.push_back((long double)i / (size - 1));
		}
		for (int i = 0; i < 3; ++i) {
			pos.push_back(pos.back() + (long double)1 / (size - 1));
		}
		Matrix<long double> A(size + 2, size + 2);
		Vector<long double> b(size + 2);
		for (int i = 0; i < size + 2; ++i) {
			if (i == size) {
				A[i][0] = B(3, -3, t(0));
				A[i][1] = B(3, -2, t(0));
				A[i][2] = B(3, -1, t(0));
				b[i] = 0;
			} else if (i == size + 1) {
				A[i][size - 1] = B(3, size - 4, t(size - 1));
				A[i][size] = B(3, size - 3, t(size - 1));
				A[i][size + 1] = B(3, size - 2, t(size - 1));
				b[i] = 1;
			} else {
				A[i][i] = B2(3, i - 3, t(i)) - 3 * B(3, i - 3, t(i));
				A[i][i + 1] = B2(3, i - 2, t(i)) - 3 * B(3, i - 2, t(i));
				A[i][i + 2] = B2(3, i - 1, t(i)) - 3 * B(3, i - 1, t(i));
				b[i] = pow(t(i), 2);
			}
		}
		cerr << A << endl;
		cerr << Transpose(b) << endl;
		auto k = LU::SolveLinearEquationSystem(A, b);
		cerr << Transpose(k) << endl;
		cerr << Transpose(A * k) << endl;
		system("pause");
		ofstream fout(filename.front());
		fout.precision(10);
		filename.erase(filename.begin());
		for (int step = 1; step <= 1000; ++step) {
			long double x = step * 1e-3;
			long double y = 0;
			for (int i = 0; i < size + 2; ++i) {
				y += B(3, i - 3, x) * k[i];
			}
			cout << x << "\t" << y << endl;
			fout << x << "\t" << y << endl;
		}
		fout.close();
		cout << "-------------------------" << endl;
	}
}

int main(int argc, char const *argv[])
{
	FiniteElementMethod();
	return 0;
}