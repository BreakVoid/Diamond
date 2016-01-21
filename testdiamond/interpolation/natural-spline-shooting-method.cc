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
	for (size_t i = 0; i < n; ++i) {
		cerr << data[i].first << " " << data[i].second << endl;
	}
	auto interpolationResult = Interpolation::NaturalSplineInterpolationEco(data);
	for (const auto &seg : interpolationResult) {
		for (const auto &para : seg) {
			cout << para << "\t";
		}
		cout << endl;
	}
	return 0;
}