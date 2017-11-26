#pragma once



#include <functional>


/* old includes before refactor */
#include "Circle.h"
/* newest includes */

#include "TranscendentalFunction.h"
#include "PolynomialFunction.h"
#include "CubicFunction.h"
#include "QuadraticFunction.h"
#include "LinearFunction.h"
#include "RootFunction.h"
#include "MathConstants.h" // for euluere
#include "RationalFunction.h"
#include "PiecewiseFunction.h"
#include "MiscMathEquations.h"
#include "Limit.h"

typedef pair<double, double> Point;


// TODO: Move the important constants into one header
// TODO: Calculate the domain and ranges of each individual function after splitting up files (and other data)
// TODO: Set the domain and range variables in the constructors of all the different function types.

// Basic Important Notes on Functions
/* 1. (algebraic functions can only have powers that are rational numbers.)
*/


void PrintPointSlopeForm(const double& Slope, const Point& InPoint);

double GetSlope(const Point& FirstPoint, const Point& SecondPoint);

void PrintSlopeInterceptForm(const Point& Point, const double& Slope);



// Start of calculus chapter 
// TODO: Sort/organize this stuff later

// Slope of a secant line formula
// through Points (a,f(a)) and (x,f(x))
inline double GetSlopeOfSecantLineTwoPoints(const Point& FirstPoint, const Point& SecondPoint)
{
	double a = FirstPoint.first;
	double FOfa = FirstPoint.second;

	double x = SecondPoint.first;
	double FOfx = SecondPoint.second;

	// use formula return result
	double Numerator = FOfx - FOfa;
	double Denominator = x - a;

	double SlopeOfSecantLine = Numerator / Denominator;

	return SlopeOfSecantLine;
}


// TODO: This functions variables need to be better understood
//  average velocity of an object over a time period to be the change in its position divided by the length of the time period.
inline double GetAverageVelocity(const double& a, const double& t,
	std::function<double(const double&)> s)
{
	// x == t
	/*
	Let s(t)s(t) be the position of an object moving along a coordinate axis at time t. 
	The average velocity of the object over a time interval [a,t] where a<t
	(or [t,a] if t<a)t<a) 
	is
	
	s(t) - s(a)
	/
	t - a

	*/
	if (a < t)
	{
		std::cout << "a < t using iterval [a, t]\n";

		// use formula return result
		double Numerator = s(t) - s(a);
		double Denominator = t - a;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}

	if (a > t)
	{

		std::cout << "t < a using iterval [t, a]\n";
		// use formula return result
		double Numerator = s(a) - s(t);
		double Denominator = a - t;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}
}

// Interval [0,3] under estimate means
// evaluate the function for 0,1,2 and add the results
// underestimate = false means evaluate for 1 2 3
inline double GetAreaUnderCurve(const int& IntervalStart, const int& IntervalEnd,
	std::function<double(const double&)> func,
	const bool& bIsUnderEstimate = true)
{
	double OutResult{ 0.0 };

	if (bIsUnderEstimate)
	{
		for (int i = IntervalStart; i < IntervalEnd; i++)
		{
			OutResult = OutResult + func(i);
		}
	}

	// TODO: is this bottom part correct?
	if (bIsUnderEstimate == false)
	{
		for (int i = IntervalStart + 1; i <= IntervalEnd; ++i)
		{
			OutResult = OutResult + func(i);
		}
	}

	return OutResult;
}

//template <typename FirstFunc, typename SecondFunc>
//bool DetermineContinunityAtAPoint(RationalFunction<FirstFunc,SecondFunc>& InFunc, const int& InPoint);

bool DetermineContinunityAtAPoint(PiecewiseFunction<QuadraticFunction, LinearFunction>& InFunc, const int& InPoint);


// TODO: How I did this is super sloppy and specific, should alter it somehow.
template <typename FirstFunc, typename SecondFunc, int ThirdFunctionConstant>
inline bool DetermineContinunityAtAPoint(const PiecewiseFunctionThreeFunctions<FirstFunc, SecondFunc, ThirdFunctionConstant>& InFunc, const int& InPoint)
{
	//QuadraticFunction FirstFuntion = InFunc.GetFirstFunction();
	//LinearFunction SecondFunction = InFunc.GetSecondFunction();

	// Step 1: check to see if f(a) is defined
	double TestOne = InFunc(InPoint);

	if (std::isnan(TestOne))
	{
		// failed 
		std::cout << "TestOneFailed: The function is not continuous at " << InPoint << "\n";
		return false;
	}
	else
	{
		//  If f(a) is defined, continue to step 2.

		// Step 2: Compute Limit from both sides
		// If Limit does not exist (that is, it is not a real number),
		// then the function is not continuous at a and the problem is solved.



		Limit TestTwoLimit(InFunc, InPoint);
		double TestTwo = TestTwoLimit.GetLimitResult();

		// if Limit Exists go to step 3
		// TODO: Need a rational function check to return bool to see if it exists
		// TODO: Need more example data in order to fix

		if (TestOne != TestTwo)
		{
			std::cout << "TestThreeFailed: The function is not continuous at " << InPoint << "\n";
			std::cout << "Because f(" << InPoint << ") " << " = " << TestOne << " != " << TestTwo << " = Limit fof(x)_x->m_a \n";
			return false;
		}
		else
		{
			std::cout << "TestThreePassed: The function is continuous at " << InPoint << "\n";
			return true;
		}

	}
}


/* Types of discontinunities 
Intuitively, a removable discontinuity is a discontinuity for which there is a hole in the graph,
a jump discontinuity is a noninfinite discontinuity for which the sections of the function do not meet up,
and an infinite discontinuity is a discontinuity located at a vertical asymptote.
*/

// If a function is not continuous then you can classify what type its discontinunity is.

// 1. Removeable if the limit exists

// 2. Jump if limit exists from both sides but they are not equal.

// 3. Infinity if limit is +- infinity on both sides

inline bool ApplyIntermediateValueTherom(const CubicFunction& InCubicFunc, const int& ClosedIntervalStart, const int& ClosedIntervalEnd)
{
	auto FullFunctionForm = InCubicFunc.GetABCD();

	//auto a = std::get<0>(FullFunctionForm);
	//auto b = std::get<1>(FullFunctionForm);
	//auto c = std::get<2>(FullFunctionForm);
	//auto d = std::get<3>(FullFunctionForm);

	double StartIntervalRes = InCubicFunc(ClosedIntervalStart);
	double EndIntervalRes = InCubicFunc(ClosedIntervalEnd);

	bool OppositeSigns = (StartIntervalRes > 0 && EndIntervalRes < 0) || (StartIntervalRes < 0 && EndIntervalRes > 0);

	if (OppositeSigns)
	{
		// There exists at least one zero in the interval
		return true;
	}
}


//inline void ProveLinearFunctionLimitEpsilonDelta(const Limit& LinearFunctionLimit)
//{
//	LinearFunction LinearFunc = LinearFunctionLimit.GetLinearFunctionIfExists();
//
//	double a = LinearFunc.GetA();
//	double b = LinearFunc.GetB();
//
//	// Let  epsilon > 0
//
//	double Epsilon = 1;
//	double Delta = 1;
//	// This means we must prove that whatever follows is true no matter what positive value of ε is chosen.
//
//	auto L = LinearFunctionLimit.GetLimitResult();
//
//	// Factor the form
//	// (ax + b) = L
//	if (L > 0)
//	{
//		b = b - L;
//	}
//	else if (L < 0)
//	{
//		b = b + L;
//	}
//	else
//	{
//		// shouldnt reach here
//		throw std::logic_error("reached invalid location Prove epsilon delta linear func");
//	}
//
//	// current form is |ax + b| < epsilon
//	// check if needs factoring
//	
//	// Should only need this variable sometimes?
//	double Divisor = 0;
//	
//	if (a == b)
//	{
//		Divisor = a;
//
//		a = a / Divisor;
//		b = b / Divisor;
//		// Epsilion = Epsilion / Divisior;
//
//		// current form
//		// std::abs(LinearFunction NewForm(a, b)) <  Epsilion / Divisior;
//		
//		// Thus, it would seem that delta = Epsilion / Divisior is appropriate.
//		// Delta is std::abs(LinearFunction NewForm(a, b)) 
//		// Delta = Epsilion / Divisior;
//
//
//	}
//
//	// Now Assume 0 < | x − 1 | < Delta
//	// When Delta has been chosen
//
//	// Then 0<|x−1|< delta ,  then |(2x+1)−3|< epsilion.
//
//	// |2| |x-1|
//	// 2 *|x-1|
//	// < 2 * Delta ~~~~~~ here’s where we use the assumption that 0<|x−1|< delta
//	// Let our choice of delta = epsilion / divisior
//	// = 2*epsilion/divisor  // if divisior == 2  = epsilion
//
//
//}
