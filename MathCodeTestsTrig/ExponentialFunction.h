#pragma once


#ifndef EXPONENTIALFUNCTION_H
#define EXPONENTIALFUNCTION_H


#include "TranscendentalFunction.h"
#include "FunctionEnums.h"
#include "MathConstants.h" // for euluere


#include <exception>
#include <cmath>

using std::exception;
using std::pow;


// Exponent Rules
// http://dailprice.weebly.com/uploads/1/3/2/7/13279202/2545731_orig.png

//  An exponential function is a function of the form f(x)= a*b^x, where a > 0
// where the base 0 < b < 1 or b > 1
class ExponentialFunction : TranscendentalFunction
{
private:
	double m_a;
	double m_b;

	// If not growth its decay
	bool m_bIsExponentialGrowthFunc = true;

	Domain m_Domain;
	Range m_Range;

	EndBehavior m_AsXGoesPosEndBehavior;
	EndBehavior m_AsXGoesNegEndBehavior;

public:
	// Exponential functions have constant bases and variable exponents
	explicit ExponentialFunction(const double& a, const double& b)
		: m_a(a), m_b(b)
	{
		// a must be positive
		const bool AIsNotGreaterThanZero = (!(a > 0));

		if (AIsNotGreaterThanZero)
			throw std::exception("a needs to be greater than 0");

		// exponential growth or exponential decay?
		if (b > 1)
		{
			// growth
			m_bIsExponentialGrowthFunc = true;
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_Range = Range::ExclusiveZeroToPosInfinity;
			// Increasing on (neg inf to pos inf)
			// b^x -> pos inf AS x -> pos inf 
			// b^x -> 0 AS x -> neg inf
			m_AsXGoesPosEndBehavior = EndBehavior::AsXGoesToPosInfinityFOfXGoesToPosInfinity;
			m_AsXGoesNegEndBehavior = EndBehavior::AsXGoesToNegInfinityFOfXGoesToZero;
		}
		else if (b > 0 && b < 1)
		{
			// decay
			m_bIsExponentialGrowthFunc = false;
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_Range = Range::ExclusiveZeroToPosInfinity;
			// decreasing on (neg inf to pos inf)
			// b^x -> 0 AS x -> pos inf 
			// b^x -> pos inf AS x -> neg inf
			m_AsXGoesPosEndBehavior = EndBehavior::AsXGoesToPosInfinityFOfXGoesToZero;
			m_AsXGoesNegEndBehavior = EndBehavior::AsXGoesToNegInfinityFOfXGoesToPosInfinity;
		}
		else if (b == 1)
		{
			throw std::exception("b cannot be == to 1");
		}
		else
		{
			// b < 0
			throw std::exception("b needs to be greater than 0");
		}

	}

};

// TODO: maybe move these 2 formula functions to a different file later.
//http://www.mathwords.com/c/continuously_compounded_interest.htm
// Continous Compound Intrest
// P = principal - Starting Amount 
// r = rate of intrest per year // Example: .032 which is 3.2% as a decimal
// t = number of years
// return A = final amount 
double FindCompoundIntrest(
	const double& P,
	const double& r,
	const double& t);

// http://www.softschools.com/formulas/math/compound_interest_formula/138/
// Compound Intrest Formula (Exponential Func)
// P = principal
// r = intrest rate as a decimal // Example: .032 which is 3.2% as a decimal
// n = number of times compounded per year
// t = number of years
// return A = the future value a particular investment will have.
double FindCompoundIntrest(
	const double& P,
	const double& r,
	const double& n,
	const double& t);


#endif