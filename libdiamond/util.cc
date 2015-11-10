#include "util.hpp"

namespace Diamond {

extern const double EPS_DOUBLE = 1e-8;
extern const long double EPS_LONG_DOUBLE = 1e-10;

bool EqualZero(const double &x)
{
	return -EPS_DOUBLE < x && x < EPS_DOUBLE;
}

bool EqualZero(const long double &x)
{
	return -EPS_LONG_DOUBLE < x && x < EPS_LONG_DOUBLE;
}

}