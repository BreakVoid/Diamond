#include "../../libdiamond/interpolation.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;
using namespace Diamond;

int main(int argc, char const *argv[])
{
	ifstream dataFile("data");
	int n;
	vector<pair<double, double>> data;
	dataFile >> n;
	for (int i = 0; i < n; ++i) {
		double t, y;
		dataFile >> t >> y;
		data.push_back(make_pair(t, y));
	}
	auto bspline = Interpolation::Bspline<double>(data);
	// for (double t = -1; t <= 1.0; t += 1e-3) {
	// 	cout << t << '\t' << bspline.y(t) << endl;
	// }
	return 0;
}