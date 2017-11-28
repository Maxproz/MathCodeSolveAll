#pragma once


#ifndef TRIGOMETRICFUNCTION_H
#define TRIGOMETRICFUNCTION_H

#include "TranscendentalFunction.h"
#include "Circle.h"

#include <vector>
#include <utility>

using std::vector;
using std::pair;


// Formally, a function f is periodic if there exists a number p such that f(x + p) = f(x) 
// for all x. The smallest possible value of p is the period.
// The reciprocal of period is frequency.
enum class Period
{
	OnePI,
	TwoPI
};


// PI/180 == 180/PI == 1
// Convert 225 degrees to radians = (225 * pi) / 180 = 5pi/4
// Convert 5pi/3 to degrees = 5pi/3 * 180/pi == 900 / 3 == 300

// A trigonometric function relates the ratios of two sides of a right triangle. 
// They are sinx,cosx,tanx,cotx,secx,andcscx.

// going to define this function using the unit circle
/*
Let P=(x,y) be a point on the unit circle centered at the origin O.
Let θ be an angle with an initial side along the positive x-axis
and a terminal side given by the line segment OP.
The trigonometric functions are then defined as

sin(theta) = y
cos(theta) = x
tan(theta) = y/x

*/


// Trig Identity Link below
// http://www.mathwords.com/t/trig_identities.htm

//  Acute Angle
// An angle that has measure less than 90°.
// SOHCAHTOA Defining only works on acute angles 

// The sine, cosine, secant, and cosecant functions have a period of 2π.
// Since the tangent and cotangent functions repeat on an interval of length π

// TODO: When im more experienced I need to figure out a way to implement the squeeze theorm for limits 
// and implement it into this function 
// https://cnx.org/contents/i4nRcikn@2.66:-xC--8XH@6/The-Limit-Laws


// Trigonometric functions are continuous over their entire domains.

class MPSIN
{
	
};

class MPCOS
{

};

class MPTAN
{

};

// TODO: Still got to add necessary checks/functionality for how to handle the input when a/b/c/d variables are not all filled out


// NOTE: These functions currently expect the full form to be imputed into the function... all 5 variables + sin/cos/tan..
template <typename T>
class TrigometricFunction
{
public:
	

};

template <>
class TrigometricFunction<MPSIN>
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

public:

	explicit TrigometricFunction(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	//double operator()(const double& x)
	//{
	//	double First = (x + m_c) * m_b;
	//	double InFunc = std::sin(First);
	//	double StretchOrShrink = InFunc * m_a;
	//	double VertShift = StretchOrShrink + m_d;

	//	return VertShift;
	//}

	// try this out for a bit? until I better understand the inverse crap
	double operator()(const double& x, const bool& bIsInverseFunction = false)
	{
		double X = x;

		if (bIsInverseFunction == true)
		{
			 X = std::pow(x, -1);
		}
		else
		{
			
		}

		double First = (X + m_c) * m_b;
		double InFunc = std::sin(First);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};

template <>
class TrigometricFunction<MPCOS>
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

public:

	explicit TrigometricFunction(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	// try this out for a bit? until I better understand the inverse crap
	double operator()(const double& x, const bool& bIsInverseFunction = false)
	{
		double X = x;

		if (bIsInverseFunction == true)
		{
			X = std::pow(x, -1);
		}
		else
		{

		}

		double First = (X + m_c) * m_b;
		double InFunc = std::cos(First);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};


template <>
class TrigometricFunction<MPTAN>
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

public:

	explicit TrigometricFunction(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	// try this out for a bit? until I better understand the inverse crap
	double operator()(const double& x, const bool& bIsInverseFunction = false)
	{
		double X = x;

		if (bIsInverseFunction == true)
		{
			X = std::pow(x, -1);
		}
		else
		{

		}

		double First = (X + m_c) * m_b;
		double InFunc = std::tan(First);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};
inline double FixAngleBetweenZeroAndTwoRad(const double & InAngle)
{
	double FixedAngle = InAngle;

	while (FixedAngle > (2 * M_PI))
	{
		FixedAngle = FixedAngle - (2 * M_PI);
	}
	while (FixedAngle < 0)
	{
		FixedAngle = FixedAngle + (2 * M_PI);
	}

	return FixedAngle;
}



//
//class TrigonometricFunction : TranscendentalFunction
//{
//private:
//	Period m_Period;
//
//	double m_PeriodB;
//
//	double m_HorizontalShifta;
//
//	double m_VerticleStretchA;
//
//	double m_VerticleShiftC;
//
//
//	//bool m_SinFunctionInUse = false;
//	//bool m_CosFunctionInUse = false;
//
//
//	//std::vector<std::pair<double, double>> m_SinFunctionPoints
//	//{
//	//	std::pair<double,double>(0,0),
//	//	std::pair<double,double>(M_PIOverTwo, 1),
//	//	std::pair<double,double>(M_PI, 0),
//	//	std::pair<double,double>(ThreeM_PIOverTwo, -1),
//	//	std::pair<double,double>(TwoM_PI, 0)
//	//};
//
//	//std::vector<std::pair<double, double>> m_CosFunctionKeyPoints
//	//{
//	//	std::pair<double,double>(0,1),
//	//	std::pair<double,double>(M_PIOverTwo, 0),
//	//	std::pair<double,double>(M_PI, -1),
//	//	std::pair<double,double>(ThreeM_PIOverTwo, 0),
//	//	std::pair<double,double>(TwoM_PI, 1)
//	//};
//
//	double FixAngleBetweenZeroAndTwoRad(const double& InAngle);
//
//
//public:
//
//	explicit TrigonometricFunction(const double& A, const double& B, const double& a, const double& C)
//		//	: m_VerticleStretchA(A), m_PeriodB(B), m_HorizontalShifta(a), m_VerticleShiftC(C)
//	{
//		// TODO: This is currently setup to automatically evaluate a sin function
//		// Add option for cos/tan/create empty etc..
//
//
//		// why when i used these two variables in the formula it didnt work???
//		m_PeriodB = TwoM_PI / std::abs(B);
//		m_VerticleStretchA = std::abs(A);
//
//
//		m_HorizontalShifta = a;
//		m_VerticleShiftC = C;
//
//		//m_SinFunctionInUse = true;
//
//		//for (const auto& Point : m_SinFunctionPoints)
//		//{
//		//	std::cout << "x: " << Point.first << " f(x): " << Point.second << std::endl;
//		//}
//
//		//for (int i = 0; i < m_SinFunctionPoints.size(); ++i)
//		//{
//		//	// This didnt work the one below did?
//		//	//m_SinFunctionPoints[i].second = ((std::sin((m_SinFunctionPoints[i].first - m_HorizontalShifta)
//		//	//		* m_PeriodB) * m_VerticleStretchA) + m_VerticleShiftC);
//
//		//	m_SinFunctionPoints[i].second =
//		//		((std::sin((m_SinFunctionPoints[i].first - a) * B) * A) + C);
//		//}
//
//		//PrintCurrentSinFunctionKeyPoints();
//	}
//
//	double SinOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);
//	double CosOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);
//	double TanOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);
//
//	//inline void PrintCurrentSinFunctionKeyPoints() const
//	//{
//	//	for (const auto& Point : m_SinFunctionPoints)
//	//		std::cout << "x: " << Point.first << " f(x): " << Point.second << std::endl;
//	//}
//
//	inline double EvaluateInverseSinComposition(const double& InSinNum)
//	{
//		UnitCircle TempUnitCircle;
//		auto Map = TempUnitCircle.GetRadianSinCosMap();
//
//		//std::cout << InSinNum << std::endl;
//
//		// Inverse sin domain [1, 1]
//		// Inverse sin range [-90,+90] in deg
//
//		//const double NegPIOverTwo = M_PIOverTwo * (-1);	
//
//		double MaxRange = ThreeM_PIOverTwo;
//
//		double OutResult(0);
//		double OutRange(0);
//
//		for (auto Num : Map)
//		{
//			if (Num.second.second == InSinNum)
//			{
//				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
//				//	<< Num.second.second << std::endl;
//
//				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//				//	<< InSinNum << std::endl;
//
//				OutResult = Num.first;
//				break;
//			}
//		}
//		bool FirstRange = (OutResult >= ThreeM_PIOverTwo && OutResult <= 2 * M_PI);
//		bool SecondRange = (OutResult >= 0 && OutResult <= M_PIOverTwo);
//
//		bool InValidRange = (FirstRange || SecondRange);
//
//		if (InValidRange)
//		{
//
//			std::cout << OutResult << std::endl;
//			std::cout << " or \n";
//			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)
//
//			return OutResult;
//		}
//		else
//		{
//			std::cout << OutResult << std::endl;
//
//			return OutResult;
//		}
//	}
//
//	inline double EvaluateInverseCosComposition(const double& InCosNum)
//	{
//		UnitCircle TempUnitCircle;
//		auto Map = TempUnitCircle.GetRadianSinCosMap();
//
//		//std::cout << InSinNum << std::endl;
//
//		// Inverse cos domain [-1, 1]
//		// Inverse cos range [0,PI] in rad
//
//		double OutResult(0);
//		double OutRange(0);
//
//		for (auto Num : Map)
//		{
//			if (Num.second.first == InCosNum)
//			{
//				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
//				//	<< Num.second.second << std::endl;
//
//				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//				//	<< InSinNum << std::endl;
//
//				OutResult = Num.first;
//				break;
//			}
//		}
//		bool FirstRange = (OutResult <= M_PI && OutResult >= 0);
//		bool InValidRange = (FirstRange); // || SecondRange);
//
//		if (InValidRange)
//		{
//
//			std::cout << OutResult << std::endl;
//			std::cout << " or \n";
//			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)
//
//			return OutResult;
//		}
//		else
//		{
//			std::cout << OutResult << std::endl;
//
//			return OutResult;
//		}
//	}
//
//	inline double EvaluateInverseTanComposition(const double& InTanNum)
//	{
//		UnitCircle TempUnitCircle;
//		auto Map = TempUnitCircle.GetRadianTanMap();
//
//
//		// Inverse tan domain R
//		// Inverse tan range [-pi/2,pi/2] in rad
//
//		double OutResult(0);
//		double OutRange(0);
//
//		for (auto Num : Map)
//		{
//			if (Num.second == InTanNum)
//			{
//				OutResult = Num.first;
//				break;
//			}
//		}
//		bool Valid = (OutResult <= ThreeM_PIOverTwo && OutResult >= M_PIOverTwo);
//		bool InValidRange = (Valid);
//
//		if (InValidRange)
//		{
//
//			std::cout << OutResult << std::endl;
//			std::cout << " or \n";
//			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)
//
//			return OutResult;
//		}
//		else
//		{
//			std::cout << OutResult << std::endl;
//
//			return OutResult;
//		}
//	}
//
//	// TODO: How do i do this?
//	double operator()(const double x) const
//	{
//
//	}
//
//};

inline double EvaluateCommonAngleSin(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, Point>  Map = TestCircle.GetRadianSinCosMap();

	double TestSinRes;

	for (const auto& Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestSinRes = Val.second.second;
		}
	}
	return TestSinRes;
}

inline double EvaluateCommonAngleCos(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, Point>  Map = TestCircle.GetRadianSinCosMap();

	double TestCosRes;

	for (const auto& Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestCosRes = Val.second.first;
		}
	}
	return TestCosRes;
}

inline double EvaluateCommonAngleTan(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, double>  Map = TestCircle.GetRadianTanMap();

	double TestTanRes(0);

	for (auto Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestTanRes = Val.second;

		}
	}
	return TestTanRes;
}


inline double COSH(const double& x)
{
	return (std::exp(x) + std::exp(-x)) / 2;
}

inline double SINH(const double& x)
{
	return (std::exp(x) - std::exp(-x)) / 2;
}

inline double TANH(const double& x)
{
	return SINH(x) / COSH(x);
}

inline double CSCH(const double& x)
{
	return 1.0 / SINH(x);
}
inline double SECH(const double& x)
{
	return 1.0 / COSH(x);
}

inline double COTH(const double& x)
{
	return COSH(x) / SINH(x);
}


// NOTE:
/*
Inverse Hyperbolic Functions
From the graphs of the hyperbolic functions,
we see that all of them are one-to-one except coshx and sechx.
If we restrict the domains of these two functions to the interval [0,∞),
then all the hyperbolic functions are one-to-one,
and we can define the inverse hyperbolic functions.
*/


// HYPERBLOIC Identitys link below
// https://image.slidesharecdn.com/lesson3derivativeofhyperbolicfunctions-150718183502-lva1-app6892/95/lesson-3-derivative-of-hyperbolic-functions-5-638.jpg?cb=1437244604
// sinh^-1 -- arcsinhx -- ln(x+sqrt(x^2+1))
inline double FindInverseHyperbolicSin(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (x + (std::sqrt(std::pow(x, 2) + 1)));
	double LastStep = std::log(InsideParathesis);

	return LastStep;

}

// cosh^-1
inline double FindInverseHyperbolicCos(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (x + (std::sqrt(std::pow(x, 2) - 1)));
	double LastStep = std::log(InsideParathesis);

	return LastStep;

}

inline double FindInverseHyperbolicTan(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = ((1 + x) / (1 - x));
	double NextStep = std::log(InsideParathesis);
	double LastStep = NextStep * 0.5;

	return LastStep;
}

inline double FindInverseHyperbolicCot(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = ((x + 1) / (x - 1));
	double NextStep = std::log(InsideParathesis);
	double LastStep = NextStep * 0.5;

	return LastStep;
}

inline double FindInverseHyperbolicSec(const double& x)
{
	// TODO: should i block the invalid input or what?
	// undefined for x <= 0, x > 1

	// setup the formula and return result
	double InsideParathesis = (1 + (std::sqrt(1 - std::pow(x, 2)))) / x;
	double LastStep = std::log(InsideParathesis);

	return LastStep;
}

inline double FindInverseHyperbolicCsc(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (1 / x) + ((std::sqrt(1 + std::pow(x, 2))) / std::abs(x));
	double LastStep = std::log(InsideParathesis);

	return LastStep;
}



#endif

