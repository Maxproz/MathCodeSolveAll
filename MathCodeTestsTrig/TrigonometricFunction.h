#pragma once


#ifndef TRIGOMETRICFUNCTION_H
#define TRIGOMETRICFUNCTION_H

#include "TranscendentalFunction.h"
#include "Circle.h"



#include <vector>
#include <utility>
#include <tuple>
#include <cmath>

using std::vector;
using std::pair;
using std::tuple;

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


//// Trigonometric functions are continuous over their entire domains.
//class TrigometricFunc
//{
//
//};
//
//
//class MPSIN : public TrigometricFunc {};
//class MPCOS : public TrigometricFunc {};
//class MPTAN : public TrigometricFunc {};
//class MPCOT : public TrigometricFunc {};
//class MPSEC : public TrigometricFunc {};
//class MPCSC : public TrigometricFunc {};
//
//// Might change how I do this later but adding classes that you get from deriviatives that are multiple trig functions
//// example f(x) = sec(x)tan(x)
//class MPSECTAN : public TrigometricFunc {};
//class MPNEGCSCCOT : public TrigometricFunc {};
//
//// Added the classes below to track the negative functions you get from derivatives
//class MPNEGSIN : public TrigometricFunc {};
//class MPNEGCOS : public TrigometricFunc {};
//class MPNEGCSC : public TrigometricFunc {};



// TODO: Still got to add necessary checks/functionality for how to handle the input when a/b/c/d variables are not all filled out
// TODO: Add this functionality listed in the following 2 lines so its tracked in the sin/cos functions
// f(x) = sinx is increasing, f′(x) = cosx > 0 
// f(x) = sinx is decreasing, f′(x) = cosx < 0.
// TODO: Add this tracking in variables in sin function listed in the 2 lines below.
// Notice that at the points where f(x) = sinx has a horizontal tangent, 
// its derivative f′(x) = cosx takes on the value zero. 
// TODO IMPORTANT: How I setup the trig functions to handle NEGCOS NEGSIN is going to cause me a lot of confusion later on if I don't handle
// how they will work with the current m_a variable and how the signs of them will interact,
// as of right now I am just using a hardcoded constant in the operator() function to represent it
// UPDATE: after looking at it again i am pretty sure it won't interact with the A variable so I should be ok to move forward,
// and do additional testing as I go on

// TODO: Make PrintFunction()'s for these classes to make them easier to debug.
// TODO: Figure out a way to make reading the type of function into string from the MPTrigFunction classes easier.

// TODO IMPORTANT: It appears I can make derivative functions for the higher ordered powers by using the method I saw on Line 539: - void AutoSetSECDerivativeFunction(MPSEC<2>& InFunc);
// I would just have to figure out what the return type would be and then make a corresponding derivative function in the derivative header.




//// NOTE: These functions currently expect the full form to be imputed into the function... all 5 variables + sin/cos/tan..
//template <typename TrigFunc, int Power>
//class TrigometricFunction
//{
//public:
//
//
//
//};
template<int POWER = 1> class MPSIN;
template<int POWER = 1> class MPCOS;
template<int POWER = 1> class MPTAN;
template<int POWER = 1> class MPNEGSIN;
template<int POWER = 1> class MPCOT;
template<int POWER = 1> class MPSEC;
template<int POWER = 1> class MPSECTAN;
template<int POWER = 1> class MPNEGCSC;
template<int POWER = 1> class MPNEGCSCCOT;
template<int POWER = 1> class MPCSC;
template<int POWER = 1> class MPNEGCOS;




template <int POWER>
class MPNEGSIN
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	//int m_Power = Power;
	int m_Power = POWER;

public:


	explicit MPNEGSIN() = default;
	MPNEGSIN(const MPNEGSIN&) = default;

	explicit MPNEGSIN(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}


	// try this out for a bit? until I better understand the inverse crap
	double operator()(const double& x, const bool& bIsInverseFunction = false)
	{
		double X = x;

		if (bIsInverseFunction == true)
		{
			X = std::pow(x, (-1));
		}
		else
		{

		}

		double First = (X + m_c) * m_b;
		double InFunc = (std::pow(std::sin(First), m_Power)) * (-1);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};


void AutoSetCOSDerivativeFunction(MPCOS<1>& InFunc);

template<int POWER>
class MPCOS
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	//int m_Power = Power;
	int m_Power = POWER;

	MPNEGSIN<1> m_DerivativeFunction = MPNEGSIN<1>(1, 1, 0, 0);

	//template <typename T>
	//inline T GetHigherOrderedDerivativeFunction(const unsigned int& DesiredDerivativeNum)
	//{


	//}

public:

	explicit MPCOS() = default;
	MPCOS(const MPCOS&) = default;
	//explicit TrigometricFunction(TrigometricFunction<MPCOS, 1>&) = default;
	//TrigometricFunction<MPCOS,1>& operator=(const TrigometricFunction<MPCOS,1>&) = default;

	 explicit MPCOS(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{
		 AutoSetCOSDerivativeFunction(*this);

	}

	 inline void SetDerivativeFunction(const MPNEGSIN<1>& InFunc)
	 {
		 //if (Power != 1)
		 //{
		 //	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
		 //}

		 m_DerivativeFunction = InFunc;
	 }


	 inline MPNEGSIN<1> GetDerivativeFunction() const
	 {
		 //if (Power != 1)
		 //{
		 //	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
		 //}

		 return m_DerivativeFunction;
	 }


	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
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

	std::string GetFunctionString() const;

};




void AutoSetSinDerivativeFunction(MPSIN<1>& InFunc);


// Now accepts functions of the form sin^2(x)
// If a power is not specified on construction it's assumed to be in the form f(x) = sin(x) with also equal to - sin^1(x)
// https://calculus.subwiki.org/wiki/Sine_function
template <int POWER>
class MPSIN
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;
	
	//int m_Power = Power;
	int m_Power = POWER;
	
	// NOTE IMPORTANT: Need to think of a way to get the derivative from higher ordered Trig functions using the Power Variable
	MPCOS<1> m_DerivativeFunction = MPCOS<1>(1,1,0,0);



public:
	
	// NOTE IMPORTANT: Need to think of a way to get the derivative from higher ordered Trig functions using the Power Variable




	explicit MPSIN() = default;

	MPSIN(const MPSIN&) = default;

	explicit MPSIN(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

		AutoSetSinDerivativeFunction(*this);
	}

	void SetDerivativeFunction(const MPCOS<1>& InFunc);
	MPCOS<1> GetDerivativeFunction() const;

	inline tuple<double, double, double, double> GetABCD() const { return tuple<double, double, double, double>(m_a, m_b, m_c, m_d); }


	
	void PrintHigherOrderedDerivativeFunctionType(const unsigned int& DesiredDerivativeNum);
	

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
			 X = std::pow(x, (-1));
		}
		else
		{
			
		}

		double First = (X + m_c) * m_b;
		double InFunc = std::pow(std::sin(First), m_Power);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

	std::string GetFunctionString() const;

};


template <int POWER>
class MPNEGCOS
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	//int m_Power = Power;
	int m_Power = POWER;

public:

	MPNEGCOS() = default;

	explicit MPNEGCOS(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{


	}

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
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
		double InFunc = (std::pow(std::cos(First), Power)) * (-1);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};


void AutoSetTANDerivativeFunction(MPTAN<1>& InFunc);

template <int POWER>
class MPTAN
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;
	
	//int m_Power = Power;
	int m_Power = POWER;
	
	MPSEC<2> m_DerivativeFunction = MPSEC<2>(1, 1, 0, 0);


public:

	MPTAN() = default;

	explicit MPTAN(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

		AutoSetTANDerivativeFunction(*this);

	}

	void SetDerivativeFunction(const MPSEC<2>& InFunc);
	MPSEC<2> GetDerivativeFunction() const;


	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
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

void AutoSetCOTDerivativeFunction(MPCOT<1>& InFunc);

template <int POWER>
class MPCOT
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	int m_Power = POWER;
	
	MPNEGCSC<2> m_DerivativeFunction = MPNEGCSC<2>(1, 1, 0, 0);

public:

	MPCOT() = default;

	explicit MPCOT(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{
		AutoSetCOTDerivativeFunction(*this);
	}

	void SetDerivativeFunction(const MPNEGCSC<2>& InFunc);
	MPNEGCSC<2> GetDerivativeFunction() const;

	inline tuple<double, double, double, double> GetABCD() const {	return tuple<double, double, double, double>(m_a, m_b, m_c, m_d); }

	std::string GetEquationForATangentLineAtInputAngle(const double& x);

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
		double InFunc = (std::pow(std::cos(First), m_Power) / std::pow(std::sin(First), m_Power));
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};

void AutoSetSECDerivativeFunction(MPSEC<1>& InFunc);
void AutoSetSECDerivativeFunction(MPSEC<2>& InFunc);

template <int POWER>
class MPSEC
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	int m_Power = POWER;
	
	MPSECTAN<1> m_DerivativeFunction = MPSECTAN<1>(1, 1, 0, 0);

public:

	MPSEC() = default;

	explicit MPSEC(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{
		AutoSetSECDerivativeFunction(*this);
	}

	void SetDerivativeFunction(const MPSECTAN<1>& InFunc);
	MPSECTAN<1> GetDerivativeFunction() const;


	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
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
		double InFunc = 1.0 / std::pow(std::cos(First), m_Power);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};

void AutoSetCSCDerivativeFunction(MPCSC<1>& InFunc);

template <int POWER>
class MPCSC
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	int m_Power = POWER;
	
	MPNEGCSCCOT<1> m_DerivativeFunction = MPNEGCSCCOT<1>(1, 1, 0, 0);

public:

	MPCSC() = default;

	explicit MPCSC(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

		AutoSetCSCDerivativeFunction(*this);

	}

	void SetDerivativeFunction(const MPNEGCSCCOT<1>& InFunc);
	MPNEGCSCCOT<1> GetDerivativeFunction() const;

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
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
		double InFunc = 1.0 / std::pow(std::sin(First), m_Power);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};



template <int POWER>
class MPNEGCSC
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	//int m_Power = Power;
	int m_Power = POWER;

public:

	MPNEGCSC() = default;

	explicit MPNEGCSC(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}

	// TODO IMPORTANT: Hmmm getting an inkling the way I am doing inverses here will fuck me in these function.
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
		double InFunc = (1.0 / std::pow(std::sin(First), m_Power)) * (-1);
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};

template <int POWER>
class MPSECTAN
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	int m_Power = POWER;

public:

	MPSECTAN() = default;

	explicit MPSECTAN(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}

	// TODO IMPORTANT: Hmmm getting an inkling the way I am doing inverses here will fuck me in these function.
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

		// the denominator power in this function is twice the Power variable
		double CosDenomPower = m_Power * 2;
		double InFunc = (std::pow(std::sin(First), m_Power)) / (std::pow(std::cos(First), CosDenomPower));
		double StretchOrShrink = InFunc * m_a;
		double VertShift = StretchOrShrink + m_d;

		return VertShift;
	}

};

template <int POWER>
class MPNEGCSCCOT
{
private:
	double m_a;
	double m_b;
	double m_c;
	double m_d;

	int m_Power = POWER;

public:

	MPNEGCSCCOT() = default;

	explicit MPNEGCSCCOT(const double& a, const double& b, const double& c, const double& d)
		: m_a{ a }, m_b{ b }, m_c{ c }, m_d{ d }
	{

	}

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}

	// TODO IMPORTANT: Hmmm getting an inkling the way I am doing inverses here will fuck me in these function.
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

		// the denominator power in this function is twice the Power variable
		double SinDenomPower = m_Power * 2;
		// Calculate it using the alternate form and then times it by negative one to get the result you need
		double InFunc = ((std::pow(std::cos(First), m_Power)) / (std::pow(std::sin(First), SinDenomPower))) * (-1);
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


template<int POWER>
inline void MPTAN<POWER>::SetDerivativeFunction(const MPSEC<2>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	m_DerivativeFunction = InFunc;

}

template<int POWER>
inline MPSEC<2> MPTAN<POWER>::GetDerivativeFunction() const
{
	return m_DerivativeFunction;
}

template<int POWER>
inline void MPCOT<POWER>::SetDerivativeFunction(const MPNEGCSC<2>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	m_DerivativeFunction = InFunc;
}

template<int POWER>
inline MPNEGCSC<2> MPCOT<POWER>::GetDerivativeFunction() const
{
	return m_DerivativeFunction;
}


template<int POWER>
inline std::string MPCOT<POWER>::GetEquationForATangentLineAtInputAngle(const double & x)
{
	double y = operator()(x);
	MPNEGCSC<2> DerivativeFunc = GetDerivativeFunction();
	double Slope = DerivativeFunc(x);
	cout << "TestSlope Value: " << endl;
	// Using the point-slope equation of the line, we obtain

	// Y - y = Slope(X - x)
	// format into std::string
	// return
	double SlopeX = Slope;


	// TODO: Something tells me I need a better system for determining the signs for these variables.

	double RHSPI = ((Slope * x) * (-1));
	double RHSY = y;
	
	//FlipSign(RHSY);


	std::string OutString("y = ");
	OutString.append(std::to_string(SlopeX));
	OutString.append("x + ");
	OutString.append(std::to_string(RHSY));
	OutString.append(" + ");
	OutString.append(std::to_string(RHSPI));

	return OutString;
}



template<int POWER>
inline void MPSEC<POWER>::SetDerivativeFunction(const MPSECTAN<1>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	m_DerivativeFunction = InFunc;
}

template<int POWER>
inline MPSECTAN<1> MPSEC<POWER>::GetDerivativeFunction() const
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	return m_DerivativeFunction;
}

template<int POWER>
inline void MPCSC<POWER>::SetDerivativeFunction(const MPNEGCSCCOT<1>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	m_DerivativeFunction = InFunc;
}

template<int POWER>
inline MPNEGCSCCOT<1> MPCSC<POWER>::GetDerivativeFunction() const
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	return m_DerivativeFunction;
}

template<int POWER>
inline void MPSIN<POWER>::SetDerivativeFunction(const MPCOS<1>& InFunc)
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	m_DerivativeFunction = InFunc;
}

template<int POWER>
inline MPCOS<1> MPSIN<POWER>::GetDerivativeFunction() const
{
	//if (Power != 1)
	//{
	//	throw std::logic_error("You have not setup taking derivatives with higher level sin function powers");
	//}

	return m_DerivativeFunction;
}

template<int POWER>
inline void MPSIN<POWER>::PrintHigherOrderedDerivativeFunctionType(const unsigned int & DesiredDerivativeNum)
{
	int Remainder = (DesiredDerivativeNum % 4);

	//double a = m_a * m_b;
	//double b = m_b;
	//double c = m_c;
	//double d = 0;

	if (Remainder == 0)
	{
		cout << "SIN";
		//return MPCOS(a, b, c, d);
	}
	else if (Remainder == 1)
	{
		cout << "COS";
		//return MPNEGSIN(a, b, c, d);
	}
	else if (Remainder == 2)
	{
		cout << "-SIN";
		//return MPNEGCOS(a, b, c, d);
	}
	else if (Remainder == 3)
	{
		cout << "-COS";
		//return MPSIN(a, b, c, d);
	}
	else if (Remainder == 4)
	{
		cout << "SIN";
		//return MPSIN(a, b, c, d);
	}
	

}

template<int POWER>
inline std::string MPSIN<POWER>::GetFunctionString() const
{
	std::string OutString;
	OutString.append(std::to_string(m_a));
	OutString.append("sin");
	OutString.append(std::to_string(m_b));
	OutString.append("(x - ");
	OutString.append(std::to_string(m_c));
	OutString.append(") + ");
	OutString.append(std::to_string(m_d));

	return OutString;
}

template<int POWER>
inline std::string MPCOS<POWER>::GetFunctionString() const
{
	std::string OutString;
	OutString.append(std::to_string(m_a));
	OutString.append("cos");
	OutString.append(std::to_string(m_b));
	OutString.append("(x - ");
	OutString.append(std::to_string(m_c));
	OutString.append(") + ");
	OutString.append(std::to_string(m_d));

	return OutString;
}


#endif

