#pragma once


#ifndef POLYNOMIALFUNCTION_H
#define POLYNOMIALFUNCTION_H

#include "FunctionEnums.h"
//#include "PowerFunction.h"
//
//#include <list>
//
//using std::list;

enum class PolynomialFunctionType
{
	LINEAR,
	RATIONAL,
	CUBIC,
	QUADRATIC,
	CONSTANT,
	POWER,
	QUARTIC,
	//ROOT, // not polynomial
};

class PolynomialFunction
{
private:

	//std::list<double> m_Coefficents;

protected:
	int m_Degree;
	bool m_bIsEvenFunction;
	EndBehavior m_EndBehavior;

	Domain m_Domain;
	Range m_Range;

	PolynomialFunctionType m_PolyFunctionType;

	

public:
	PolynomialFunction() = default;

	PolynomialFunction(const PolynomialFunction&) = default;

	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }

	//explicit PolynomialFunction()
	//{
	//	int i = 0;
	//	while (HighestDegree - i != 0)
	//	{
	//		

	//		i = i + 1;
	//	}
	//}


};


#endif