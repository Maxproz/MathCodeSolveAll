#pragma once


#ifndef POLYNOMIALFUNCTION_H
#define POLYNOMIALFUNCTION_H

#include "FunctionEnums.h"

enum class PolynomialFunctionType
{
	LINEAR,
	RATIONAL,
	CUBIC,
	QUADRATIC,
	ROOT, // not polynomial
};

class PolynomialFunction
{
private:


protected:
	int m_Degree;
	bool m_bIsEvenFunction;
	EndBehavior m_EndBehavior;

	Domain m_Domain;
	Range m_Range;

	PolynomialFunctionType m_PolyFunctionType;

public:

	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }

};


#endif