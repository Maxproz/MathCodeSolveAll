#pragma once


#ifndef CUBICFUNCTION_H
#define CUBICFUNCTION_H

#include "PolynomialFunction.h"
#include "FunctionEnums.h"

#include <utility>
#include <cmath>
#include <tuple>

using std::pair;
using std::tuple;
using std::pow;

// The function of the coefficient aa in the general equation is to make the graph 
// "wider" or "skinnier", or to reflect it (if negative):
// 
//  f(x)=ax^3+bx^2+cx+d

// its natural domain and range including all real numbers x and y respectively. 
// TODO: Firstly, a cubic equation always has its center point on the vertical line x = -b/3*a
// TODO: when the quadratic term bx^2is not present in the cubic, the inflection point is on the y axis.
// TODO: We also know that a cubic function can have at most three x-intercepts.
// TODO: The constant d in the equation is the y -intercept of the graph. at  x = d
// TODO: Fifthly, we know that if the coefficient of x^3 is negative, then the cubic, notwithstanding any local extrema, 
// is falling with a negative slope from top left to bottom right. If the coefficient of x^3 is positive, 
// the curve is generally rising with positive slope from bottom left to top right.



class CubicFunction : public PolynomialFunction
{
private:
	double m_a = 0;
	double m_b = 0;
	double m_c = 0;
	double m_d = 0;

	pair<double, double> m_JustAAndDCubic = pair<double, double>(0, 0);
	tuple<double, double, double> m_ACDForm;


	bool m_bIsInAAndDForm = false;
	bool m_bIsInACDForm = false;

	inline tuple<double, double, double> GetACD() const
	{
		return tuple<double, double, double>(m_a, m_c, m_d);
	}

public:

	CubicFunction() = default;

	explicit CubicFunction(const double& a, const double& b, const double& c, const double& d)
		: m_a(a), m_b(b), m_c(c), m_d(d)
	{
		m_PolyFunctionType = PolynomialFunctionType::CUBIC;

		// TODO: might have to change this later?
		// Linear functions have all real number domain/ranges
		m_Domain = Domain::NegInfinityToPosInfinity;
		m_Range = Range::NegInfinityToPosInfinity;

		if (b == 0 && c == 0)
		{
			m_JustAAndDCubic = pair<double, double>(a, d);
			m_bIsInAAndDForm = true;
		}

		if (b == 0)
		{
			m_ACDForm = GetACD();
			m_bIsInACDForm = true;
		}

	}

	inline pair<double, double> GetAAndDCubicFuncForm() const { return m_JustAAndDCubic; }

	inline tuple<double, double, double, double> GetABCD() const
	{
		return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}



	bool GetIsFuncInAAndDForm() const { return m_bIsInAAndDForm; }
	bool GetIsFuncInACDForm() const { return m_bIsInACDForm; }


	tuple<double, double, double> GetACDForm() const { return m_ACDForm; }


	// Runs the full function form version
	double operator()(const double& x) const
	{
		double TermOne, TermTwo, TermThree, TermFour;
		TermOne = (m_a * (std::pow(x, 3)));
		TermTwo = (m_b * (std::pow(x, 2)));
		TermThree = x*m_c;
		TermFour = m_d;

		return TermOne + TermTwo + TermThree + TermFour;
	}


};



#endif





