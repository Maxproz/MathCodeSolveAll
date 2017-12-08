#pragma once


#ifndef CUBICFUNCTION_H
#define CUBICFUNCTION_H

#include "PolynomialFunction.h"
#include "FunctionEnums.h"
#include "QuadraticFunction.h"



#include <utility>
#include <cmath>
#include <tuple>
#include <vector>

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

// TODO: https://math.stackexchange.com/questions/516662/find-the-equations-of-the-tangent-lines-that-are-parallel-to-the-line-y-4x
// Find the equations of the tangent line(s) to (*this) cubic
// that are parallel to the (given linearfunction)
// Set the Derivative of the cubic equal to the derivative of the given linear function this
// makes it so that the tangent line to that point will be parallel to your line now solve for x.
// Plug the x values into the function to get the corresponding y values
// Now that you have cordinate pairs, use those pairs in the point slope form with the linear function
// derivative found above to solve for the functions.


class CubicFunction;

void AutoSetCubicDerivativeFunction(CubicFunction& InFunc);

class CubicFunction : public PolynomialFunction
{
private:
	double m_a = 0;
	double m_b = 0;
	double m_c = 0;
	double m_d = 0;

	bool m_bIsInAAndDForm = false;
	bool m_bIsInACDForm = false;


	pair<double, double> m_JustAAndDCubic = pair<double, double>(0, 0);
	tuple<double, double, double> m_ACDForm;
	

	// The function has horizontal tanget lines at these x values
	std::vector<double> m_HorizontalTangentLines;
	
	std::vector<std::pair<double, double>> m_ConcaveUpIntervals;
	std::vector<std::pair<double, double>> m_ConcaveDownIntervals;

	std::vector<Point> m_InflectionPoints;

	QuadraticFunction m_DerivativeFunction;// = QuadraticFunction(1, 0, 0);

	inline tuple<double, double, double> GetACD() const { return tuple<double, double, double>(m_a, m_c, m_d); } // pretty sure this is not needed

	inline void AutoSetHorizontalTangetLines()
	{
		std::vector<double> DerivZerosVec = m_DerivativeFunction.GetAllZerosVec();
		m_HorizontalTangentLines = (DerivZerosVec);
	}

	void AutoSetCubicInflectionPoints();
	void AutoSetConcaveUpAndDownIntervals();

	virtual void FindCriticalPoints() override;
	
	virtual void SetDefaultDomainInterval() override;
	virtual void SetDefaultRangeInterval() override;

	virtual void SetIncreasingDecreasingIntervals() override;




public:

	void FindGlobalExtremums();

	void PrintConcaveIntervalData() const;


	std::vector<Point> GetInflectionPoints() const { return m_InflectionPoints; }
	std::vector<std::pair<double, double>> GetConcaveUpIntervals() const { return m_ConcaveUpIntervals; }
	std::vector<std::pair<double, double>> GetConcaveDownIntervals() const { return m_ConcaveDownIntervals; }
	void AddAdditionalConcaveUpInterval(const std::pair<double, double>& Interval) { m_ConcaveUpIntervals.push_back(Interval); }
	void AddAdditionalConcaveDownInterval(const std::pair<double, double>& Interval) { m_ConcaveDownIntervals.push_back(Interval); }

	CubicFunction() = default;
	//CubicFunction(const CubicFunction&) = default;
	//CubicFunction(CubicFunction &&data) = default;
	//CubicFunction& operator=(CubicFunction&&data) = default;
	////QuadraticFunction(QuadraticFunction&&) = default;

	explicit CubicFunction(const double& a, const double& b, const double& c, const double& d)
		: m_a(a), m_b(b), m_c(c), m_d(d)
	{
		m_PolyFunctionType = PolynomialFunctionType::CUBIC;

		SetDefaultDomainInterval();
		SetDefaultRangeInterval();


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

		SetDegree(3);
		SetLeadingCoefficent(m_a);
		SetPolynomialEndBehaviours();
		
		AutoSetCubicDerivativeFunction(*this);

		AutoSetCubicInflectionPoints();
		// Here I will Use the Derivative function to automatically set the zeros for horizontal tangent liens
		AutoSetHorizontalTangetLines();


		m_AbsoluteMaximum = NAN;
		m_AbsoluteMinimum = NAN;

		FindCriticalPoints();
	}
	
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


	inline pair<double, double> GetAAndDCubicFuncForm() const { return m_JustAAndDCubic; }
	inline tuple<double, double, double, double> GetABCD() const{	return tuple<double, double, double, double>(m_a, m_b, m_c, m_d);}
	inline tuple<double, double, double> GetACDForm() const { return m_ACDForm; }
	
	inline QuadraticFunction GetDerivativeFunction() const { return m_DerivativeFunction; }
	inline void SetDerivativeFunction(QuadraticFunction& InFunc) { m_DerivativeFunction = InFunc; }

	inline bool GetIsFuncInAAndDForm() const { return m_bIsInAAndDForm; }
	inline bool GetIsFuncInACDForm() const { return m_bIsInACDForm; }


	std::vector<double> GetZerosOfDerivative() const { return m_HorizontalTangentLines; }

	bool IsCubicConcaveUpOverInterval(const double& ClosedIntervalStart,
		const double& ClosedIntervalEnd);
	bool IsCubicConcaveDownOverInterval(const double& ClosedIntervalStart,
		const double& ClosedIntervalEnd);



	
	
	std::string GetFunctionString() const;
	void PrintFunction() const;
	void PrintHorizontalTangetLineXValues() const;
};




#endif





