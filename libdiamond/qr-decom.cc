#include "qr-decom.hpp"

namespace Diamond {
namespace QR {

int Sgn(const double &x)
{
	if (x < -EPS_DOUBLE) {
		return -1;
	} else {
		return 1;
	}
}

int Sgn(const long double &x)
{
	if (x < -EPS_LONG_DOUBLE) {
		return -1;
	} else {
		return 1;
	}
}


}
}