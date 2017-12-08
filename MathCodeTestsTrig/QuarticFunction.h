#pragma once

#ifndef QUARTICFUNCTION_H
#define QUARTICFUNCTION_H




// In algebra, a quartic function is a function of the form
// f(x) = ax^4+bx^3+cx^2+dx+e

#include "PolynomialFunction.h"

#include <cmath>
#include <exception>
#include <tuple>

using std::pow;
using std::exception;
using std::tuple;


class QuarticFunction : public PolynomialFunction
{

private:
	double m_a;
	double m_b; 
	double m_c;
	double m_d; 
	double m_e;


	virtual void FindCriticalPoints() override;

	virtual void SetDefaultDomainInterval() override;
	virtual void SetDefaultRangeInterval() override;
	virtual void SetIncreasingDecreasingIntervals() override;

protected:
	

public:

	explicit QuarticFunction(const double& a, const double& b, const double& c, const double& d, const double& e)
		: m_a(a), m_b(b), m_c(c), m_d(d), m_e(e)
	{
		m_Degree = 4;
		m_PolyFunctionType = PolynomialFunctionType::QUARTIC;

		if (m_a == 0)
		{
			throw exception("a cannot == 0");
		}


	}


	QuarticFunction operator+(QuarticFunction const& rhs) const;

	// Runs the full function form version
	double operator()(const double& x) const
	{
		double TermOne, TermTwo, TermThree, TermFour, TermFive;
		
		TermOne = (m_a * (std::pow(x, 4)));
		TermTwo = (m_b * (std::pow(x, 3)));
		TermThree = (m_c * (std::pow(x, 2)));
		TermFour = m_d * x;
		TermFive = m_e;

		return TermOne + TermTwo + TermThree + TermFour + TermFive;
	}

	inline tuple<double, double, double, double, double> GetABCDE() const
	{
		return tuple<double, double, double, double, double>(m_a, m_b, m_c, m_d, m_e);
	}

	void PrintFunction() const;

};



#endif