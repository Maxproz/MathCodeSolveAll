#pragma once

#ifndef CONSTANTFUNCTION_H
#define CONSTANTFUNCTION_H


#include "PolynomialFunction.h"


class LinearFunction;
class QuadraticFunction;


class ConstantFunction : public PolynomialFunction
{
private: 
	double m_b{ 0 };


	virtual void FindCriticalPoints() override;

public:

	ConstantFunction() = default;

	ConstantFunction(const ConstantFunction&) = default;

	explicit ConstantFunction(const double& b)
		: m_b{ b }
	{
		m_PolyFunctionType = PolynomialFunctionType::CONSTANT;


		SetDegree(0);
	}
	
	double operator()(const double& x)
	{
		return m_b;
	}


	void PrintConstantFunctionInfo() const;

	inline double GetB() const { return m_b; }
	
	QuadraticFunction operator*(const QuadraticFunction& rhs) const;
	LinearFunction operator*(const LinearFunction& rhs) const;

};


#endif