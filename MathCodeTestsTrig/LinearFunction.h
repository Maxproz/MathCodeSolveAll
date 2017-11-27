#pragma once

#ifndef LINEARFUNCTION_H
#define LINEARFUNCTION_H

#include "PolynomialFunction.h"

#include <utility>
#include <exception>
#include <iostream>

using std::pair;
using std::exception;
using std::cout;
using std::endl;


// f(x)=ax+b
class LinearFunction : public PolynomialFunction
{
private:

	double m_a;
	double m_b;

	LineBehavior m_LineBehavior;


	float m_Slope;

	pair<double, double> m_YIntercept;
	pair<double, double> m_XIntercept;

	bool m_bIsBOnlyForm = false;

public:

	// if a > b etc.. TODO:

	// TODO: finish filling out the specific information gathering functions

	LinearFunction() = default;

	LinearFunction(const LinearFunction&) = default;

	//LinearFunction(LinearFunction&&) = default;

	//LinearFunction& operator =(const LinearFunction &) = default;

	explicit LinearFunction(double a, double b) : m_a(a), m_b(b)
	{
		m_PolyFunctionType = PolynomialFunctionType::LINEAR;

		// TODO: might have to change this later?
		// Linear functions have all real number domain/ranges
		m_Domain = Domain::NegInfinityToPosInfinity;
		m_Range = Range::NegInfinityToPosInfinity;

		m_Slope = a;

		if (m_Slope == 0)
		{
			m_Degree = 0;
			m_bIsBOnlyForm = true;
		}
		else if (m_Slope != 0)
		{
			m_Degree = 1;
		}
		else
		{
			throw std::exception("Something went wrong with assigning the degree");
		}

		m_XIntercept.first = (-1 * b) / m_Slope;
		m_XIntercept.second = 0;

		m_YIntercept.first = 0;
		m_YIntercept.second = b;

		if (a > 0)
		{
			// Line rises as x increases
			m_LineBehavior = LineBehavior::Increasing;


		}
		else if (a < 0)
		{
			// line falls as x increases
			m_LineBehavior = LineBehavior::Decreasing;
		}
		else
		{
			// a == 0
			// horizontal line
			m_LineBehavior = LineBehavior::Horizontal;
		}

	}

	double operator()(double x) const { return ((m_a*x) + m_b); }

	double GetA() const { return m_a; }
	double GetB() const { return m_b; }

	// TODO: add functioanlity for verticle lines... 
	//ax + by = c,
	//	ax + by = c,
	//where a, ba, b are both not zero, to denote the standard form of a line.

	inline void PrintLinearFunctionInfo() const { cout << " f(x) = " << m_a << "x^2 + " << m_b << endl; }

	inline bool IsBOnlyForm() const { return m_bIsBOnlyForm; }


};




#endif





