#pragma once

#ifndef CONSTANTFUNCTION_H
#define CONSTANTFUNCTION_H


#include "PolynomialFunction.h"

#include <iostream>

using std::cout;
using std::endl;


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


	inline void PrintConstantFunctionInfo() const
	{
		cout << "f(x) = " << m_b << endl;
	}

	inline double GetB() const { return m_b; }
	
	inline QuadraticFunction operator*(const QuadraticFunction& rhs) const 
	{
		
		double OutQuadA(0);
		double OutQuadB(0);
		double OutQuadC(0);

		OutQuadA = m_b * rhs.m_a;

		OutQuadB = m_b * rhs.m_b;

		OutQuadC = m_b * rhs.m_c;

		return QuadraticFunction(OutQuadA, OutQuadB, OutQuadC);
	}

	inline LinearFunction operator*(const LinearFunction& rhs) const
	{

		double OutLinearA(0);
		double OutLinearB(0);
		
		OutLinearA = m_b * rhs.m_a;
		OutLinearB = m_b * rhs.m_b;

		return LinearFunction(OutLinearA, OutLinearB);
	}

};


#endif