#pragma once

#ifndef LINEARFUNCTION_H
#define LINEARFUNCTION_H

#include "PolynomialFunction.h"
#include "MiscMathEquations.h"
#include "MathConstants.h"



#include <utility>
#include <exception>
#include <iostream>
#include <vector>

using std::pair;
using std::exception;
using std::cout;
using std::endl;
using std::vector;


class QuarticFunction;
class QuadraticFunction;
class CubicFunction;


// f(x)=ax+b
class LinearFunction : public PolynomialFunction
{
private:

	friend class ConstantFunction;


	double m_a;
	double m_b;

	LineBehavior m_LineBehavior;


	float m_Slope;

	pair<double, double> m_YIntercept;
	pair<double, double> m_XIntercept;

	bool m_bIsBOnlyForm = false;

	//std::vector<double> m_AllRealZeros;
	std::vector<double> m_AllZeros;


	// TODO: Validate this function
	inline void FindAndStoreAllRealZeros()
	{
		// form y = ax + b

		// set linear func == to 0 solve
		double LocalA = m_a;
		double LocalB = m_b;

		//bool bIsANegative = m_a < 0;
		bool bIsBNegative = m_b < 0.0;

		double LeftHandSide = LocalA; // * x
		
		double RightHandSide{ 0.0 };


		if (bIsBNegative)
		{
			RightHandSide = RightHandSide + (LocalB * (-1));
		}
		else
		{
			RightHandSide = RightHandSide + (LocalB * (-1));
		}


		double RealZero = RightHandSide / LocalA;

		// TODO: will this give me issues later?
		const bool bIsANumber = (!(std::isnan(RealZero)));
		if (bIsANumber)
		{
			//m_AllRealZeros.push_back(RealZero);
			m_AllZeros.push_back(RealZero);
		}
	}

public:

	//inline std::vector<double> GetRealNumberZerosVec() const { return m_AllRealZeros; }
	inline std::vector<double> GetAllZerosVec() const { return m_AllZeros; }

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


		FindAndStoreAllRealZeros();

	}

	// multiplication
	QuarticFunction operator*(CubicFunction const& rhs) const;
	QuadraticFunction operator*(LinearFunction const& rhs) const;
	
	//template <typename T>

	inline LinearFunction operator-(const LinearFunction& rhs) const
	{
		double OutA(0);
		double OutB(0);

		OutA = m_a - rhs.m_a;
		OutB = m_b - rhs.m_b;

		//if (m_a == 0)
		//{
		//	return ConstantFunction(OutB);
		//}
		//else
		//{
	
		//}
		return LinearFunction(OutA, OutB);
	}


	// Use this operator as squared
	QuadraticFunction GetSquaredFunction() const;

	/*QuarticFunction& operator*=(CubicFunction const& rhs);*/

	double operator()(double x) const { return ((m_a*x) + m_b); }

	double GetA() const { return m_a; }
	double GetB() const { return m_b; }

	// TODO: add functioanlity for verticle lines... 
	//ax + by = c,
	//	ax + by = c,
	//where a, ba, b are both not zero, to denote the standard form of a line.

	inline void PrintFunction() const 
	{
		//cout << "f(x) = " << m_a << "x";
		if (m_a != 0)
			cout << m_a << "x";

		if (m_b == 0)
		{
			//cout << endl;
			return;
		}
		else
		{
			char PlusOrMinus;
			if (m_b < 0)
			{
				PlusOrMinus = ' ';
			}
			else
			{
				PlusOrMinus = '+';
			}

			cout << PlusOrMinus << m_b;
		}

	}

	inline bool IsBOnlyForm() const { return m_bIsBOnlyForm; }

	void PrintAllZeros() const;
};

inline void ConvertFromPointSlopeFromToLinearForm(const Point& InPoint, const double& InSlope, LinearFunction& OutTangentLine)
{
	const double Localx = InPoint.first * (-1);
	const double Localy = InPoint.second * (-1);

	// Convert this form into a working LinearFunction result you can return.
	// Y - y = InSlope(X - x)

	double FarRHS = (InSlope * (Localx));
	double ChangedSideY = Localy;
	FlipSign(ChangedSideY);
	FarRHS = FarRHS + ChangedSideY;

	OutTangentLine = LinearFunction(InSlope, FarRHS);
}


#endif





