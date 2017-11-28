#pragma once


#ifndef DERIVATIVE_H
#define DERIVATIVE_H


#include "Limit.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "ConstantFunction.h"

template <typename InFunction, typename OutFunction>
class Derivative
{
private:
	InFunction m_InFunction;
	OutFunction m_OutFunction;




	// Evaluate a quadratic function derivative get a linear function
	inline LinearFunction EvaluateFunctionDerivative(QuadraticFunction& InFunction)
	{
		// When getting derivative here the c variable doesnt matter
		auto AB = InFunction.GetAB();

		double a = std::get<0>(AB);
		double b = std::get<1>(AB);

		// need to change for higher powered functions
		a = 2 * a;

		// loses x variable
		b = b;

		LinearFunction OutFunc(a, b);

		return OutFunc;
	}

	inline ConstantFunction EvaluateFunctionDerivative(LinearFunction& InFunction)
	{
		// When getting derivative here the c variable doesnt matter
		auto a = InFunction.GetA();
		auto b = InFunction.GetB();

		// a loses x variable
		a = a;

		// b is dropped from derivative function
		b = 0;

		ConstantFunction OutFunc(a);

		return OutFunc;
	}


public:

	explicit Derivative(InFunction& InFunc)
	{
		m_InFunction = std::move(InFunc);

		m_OutFunction = EvaluateFunctionDerivative(InFunc);
	}

	inline OutFunction GetDerivativeFunction() const { return m_OutFunction; }

	//double GetDerivative(const double& x)
	//{
	//	const unsigned int h = 0;

	//	//Limit<Function> LocalLimit(InFunc, h);
	//	double Numerator = m_Function(x + h) - m_Function(x);
	//	double Denominator = h;

	//	return Numerator / Denominator;

	//}

};



#endif