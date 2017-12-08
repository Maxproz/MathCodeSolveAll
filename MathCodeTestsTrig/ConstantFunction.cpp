#include "ConstantFunction.h"

#include "LinearFunction.h"
#include "QuadraticFunction.h"

#include <iostream>

using std::cout;
using std::endl;


void ConstantFunction::FindCriticalPoints()
{

}

void ConstantFunction::SetDefaultDomainInterval()
{
}

void ConstantFunction::SetDefaultRangeInterval()
{
}

void ConstantFunction::SetIncreasingDecreasingIntervals()
{

}

QuadraticFunction ConstantFunction::operator*(const QuadraticFunction & rhs) const
{
	double OutQuadA(0);
	double OutQuadB(0);
	double OutQuadC(0);

	OutQuadA = m_b * rhs.m_a;

	OutQuadB = m_b * rhs.m_b;

	OutQuadC = m_b * rhs.m_c;

	return QuadraticFunction(OutQuadA, OutQuadB, OutQuadC);
}

LinearFunction ConstantFunction::operator*(const LinearFunction & rhs) const
{
	double OutLinearA(0);
	double OutLinearB(0);

	OutLinearA = m_b * rhs.m_a;
	OutLinearB = m_b * rhs.m_b;

	return LinearFunction(OutLinearA, OutLinearB);
}


void ConstantFunction::PrintConstantFunctionInfo() const
{
	std::cout << "f(x) = " << m_b << std::endl;
}