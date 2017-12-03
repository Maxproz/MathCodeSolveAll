#include "TrigonometricFunction.h"










//double FixAngleBetweenZeroAndTwoRad(const double & InAngle)
//{
//	double FixedAngle = InAngle;
//
//	while (FixedAngle > (2 * M_PI))
//	{
//		FixedAngle = FixedAngle - (2 * M_PI);
//	}
//	while (FixedAngle < 0)
//	{
//		FixedAngle = FixedAngle + (2 * M_PI);
//	}
//
//	return FixedAngle;
//}

//double TrigonometricFunction::SinOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
//{
//	double FixedAngle = InAngle;
//
//	while (FixedAngle > (2 * M_PI))
//	{
//		FixedAngle = FixedAngle - (2 * M_PI);
//	}
//	while (FixedAngle < 0)
//	{
//		FixedAngle = FixedAngle + (2 * M_PI);
//	}
//
//	auto Map = InUnitCircle.GetRadianSinCosMap();
//
//	double SinOfAngle = Map[FixedAngle].second;
//
//	return SinOfAngle;
//}
//
//double TrigonometricFunction::CosOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
//{
//	double FixedAngle = InAngle;
//
//	while (FixedAngle > (2 * M_PI))
//	{
//		FixedAngle = FixedAngle - (2 * M_PI);
//	}
//	while (FixedAngle < 0)
//	{
//		FixedAngle = FixedAngle + (2 * M_PI);
//	}
//
//	auto Map = InUnitCircle.GetRadianSinCosMap();
//
//	double CosOfAngle = Map[FixedAngle].first;
//
//	return CosOfAngle;
//}
//
//double TrigonometricFunction::TanOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
//{
//	double FixedAngle = InAngle;
//
//	while (FixedAngle > (2.0 * M_PI))
//	{
//		FixedAngle = FixedAngle - (2.0 * M_PI);
//	}
//	while (FixedAngle < 0)
//	{
//		FixedAngle = FixedAngle + (2.0 * M_PI);
//	}
//
//	auto Map = InUnitCircle.GetRadianSinCosMap();
//
//
//	double CosOfAngle = Map[FixedAngle].first;
//	double SinOfAngle = Map[FixedAngle].second;
//
//
//	// TODO: used this as debug code (remove later once you memorize numeric limits stuff)
//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
//	//	<< SevenM_PIOverFour << std::endl;
//
//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//	//	<< FixedAngle << std::endl;
//
//	return SinOfAngle / CosOfAngle;
//}

void AutoSetDerivativeFunction(TrigometricFunction<MPSIN, 1>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	Derivative<TrigometricFunction<MPSIN, 1>, TrigometricFunction<MPCOS, 1>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}
