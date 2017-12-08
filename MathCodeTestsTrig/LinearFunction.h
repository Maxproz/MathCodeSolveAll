#pragma once

#ifndef LINEARFUNCTION_H
#define LINEARFUNCTION_H

#include "PolynomialFunction.h"
#include "MiscMathEquations.h"
#include "MathConstants.h"
#include "FunctionEnums.h"

#include "ConstantFunction.h"


#include <utility>
#include <exception>
#include <iostream>
#include <vector>
#include <tuple>

using std::pair;
using std::exception;
using std::cout;
using std::endl;
using std::vector;


class QuarticFunction;
class QuadraticFunction;
class CubicFunction;
class LinearFunction;


void AutoSetDerivativeFunction(LinearFunction& InFunc);

// f(x)=ax+b
class LinearFunction : public PolynomialFunction
{
private:
	friend class ConstantFunction;

	double m_a;
	double m_b;

	float m_Slope;

	pair<double, double> m_YIntercept;
	pair<double, double> m_XIntercept;

	bool m_bIsConstantFunction = false;

	//std::vector<double> m_AllRealZeros;
	std::vector<double> m_AllZeros;

	LineBehavior m_LineBehavior;
	
	// Domain Ranges specific to LinearFunctions 	//  set of real numbers
	// for linear functions (NEGINFINITY, POSINFINITY) 


	// TODO: Validate this function
	void FindAndStoreAllRealZeros();


	void AutoSetXAndYIntercepts(const double& b);
	void AutoSetSlopeAndFunctionForm(const double& a);
	void AutoSetFunctionsLineBehaviour(const double& a);
	void AutoCheckSetDegree(const double& a);


	//void AutoSetDomainInterval();
	//void AutoSetRangeInterval();
	
	ConstantFunction m_DerivativeFunction = ConstantFunction(1);


	virtual void FindCriticalPoints() override;

	virtual void SetDefaultDomainInterval() override;
	virtual void SetDefaultRangeInterval() override;
	virtual void SetIncreasingDecreasingIntervals() override;
public:

	//void CheckIsContinuousFunction();

	//inline std::vector<double> GetRealNumberZerosVec() const { return m_AllRealZeros; }
	inline std::vector<double> GetAllZerosVec() const { return m_AllZeros; }

	// if a > b etc.. TODO:

	// TODO: finish filling out the specific information gathering functions

	LinearFunction() = default;
	LinearFunction(const LinearFunction&) = default;
	//LinearFunction(LinearFunction&&) = default;
	//LinearFunction& operator =(const LinearFunction &) = default;

	explicit LinearFunction(const double& a, const double& b) : m_a(a), m_b(b)
	{
		m_PolyFunctionType = PolynomialFunctionType::LINEAR;


		SetDefaultDomainInterval();
		SetDefaultRangeInterval();


		SetLeadingCoefficent(m_a);
		AutoCheckSetDegree(a);
		SetPolynomialEndBehaviours(); // This function is only safe to run after the degree and leading coefficent have been set
	

		AutoSetSlopeAndFunctionForm(a);
		AutoSetXAndYIntercepts(b);
		AutoSetFunctionsLineBehaviour(a);

		// Normal polynomials are continuous 
		SetIsContinuousFunction(true);

		AutoSetDerivativeFunction(*this);

		FindAndStoreAllRealZeros();

	}

	// multiplication
	QuarticFunction operator*(CubicFunction const& rhs) const;
	QuadraticFunction operator*(LinearFunction const& rhs) const;
	
	// subtraction
	LinearFunction operator-(const LinearFunction& rhs) const;


	// Use this operator as squared
	QuadraticFunction GetSquaredFunction() const;


	double operator()(const double& x) const { return ((m_a*x) + m_b); }

	double GetA() const { return m_a; }
	double GetB() const { return m_b; }

	// UPDATE: what....
	// TODO: add functioanlity for verticle lines...  
	//ax + by = c,
	//	ax + by = c,
	//where a, ba, b are both not zero, to denote the standard form of a line.

	void PrintFunction() const;
	void PrintAllZeros() const;

	// Should I change the way I do these bools?
	inline bool IsConstantFunction() const { return m_bIsConstantFunction; }

	std::string GetFunctionString() const;


	ConstantFunction GetDerivativeFunction() const;
	void SetDerivativeFunction(ConstantFunction& InFunc);


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





