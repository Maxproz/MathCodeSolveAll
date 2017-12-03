#include "TrigonometricFunction.h"

#include "Derivative.h"


void AutoSetCOSDerivativeFunction(MPCOS<1>& InFunc)
{
	Derivative<MPCOS<1>, MPNEGSIN<1>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}

void AutoSetSinDerivativeFunction(MPSIN<1>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	Derivative<MPSIN<1>, MPCOS<1>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}

void AutoSetTANDerivativeFunction(MPTAN<1>& InFunc)
{
	Derivative<MPTAN<1>, MPSEC<2> > Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}

void AutoSetCOTDerivativeFunction(MPCOT<1>& InFunc)
{
	Derivative<MPCOT<1>, MPNEGCSC<2>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}

void AutoSetSECDerivativeFunction(MPSEC<1>& InFunc)
{
	Derivative<MPSEC<1>, MPSECTAN<1>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}

void AutoSetSECDerivativeFunction(MPSEC<2>& InFunc)
{
	//Derivative<MPSEC<2>, TODO: Find out and finish function > Derivative(InFunc);
	//InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());



}


void AutoSetCSCDerivativeFunction(MPCSC<1>& InFunc)
{
	Derivative<MPCSC<1>, MPNEGCSCCOT<1>> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());
}




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

