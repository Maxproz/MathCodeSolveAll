#include "ExponentialFunction.h"




double FindCompoundIntrest(const double & P, const double & r, const double & t)
{
	const double Exponent = r*t;
	const double Base = Eulere;

	double TempResult = std::pow(Base, Exponent);
	TempResult = TempResult * P;

	return TempResult;
}

double FindCompoundIntrest(
	const double& P,
	const double& r,
	const double& n,
	const double& t)
{
	const double Exponent = n*t;
	const double InsideParanthesis = 1.0 + (r / n);

	double TempResult = std::pow(InsideParanthesis, Exponent);
	TempResult = TempResult * P;

	return TempResult;

}
