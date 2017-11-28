#pragma once

#ifndef CONSTANTFUNCTION_H
#define CONSTANTFUNCTION_H


#include "PolynomialFunction.h"


class ConstantFunction : public PolynomialFunction
{
private: 
	double m_b{ 0 };

public:

	ConstantFunction() = default;

	ConstantFunction(const ConstantFunction&) = default;

	explicit ConstantFunction(const double& b)
		: m_b{ b }
	{
		m_PolyFunctionType = PolynomialFunctionType::CONSTANT;


		m_Degree = 0;
	}
	
	double operator()(const double& x)
	{
		return m_b;
	}
	
};


#endif