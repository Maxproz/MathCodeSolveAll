#include "QuarticFunction.h"
#include "MiscMathEquations.h"


#include <iostream>


using std::cout;
using std::endl;

void QuarticFunction::FindCriticalPoints()
{

}

void QuarticFunction::SetDefaultDomainInterval()
{
}

void QuarticFunction::SetDefaultRangeInterval()
{
}

void QuarticFunction::SetIncreasingDecreasingIntervals()
{

}

QuarticFunction QuarticFunction::operator+(QuarticFunction const & rhs) const
{
	double NewA = m_a + rhs.m_a;
	double NewB = m_b + rhs.m_b;
	double NewC = m_c + rhs.m_c;
	double NewD = m_d + rhs.m_d;
	double NewE = m_e + rhs.m_e;

	return QuarticFunction(NewA, NewB, NewC, NewD, NewE);
}

void QuarticFunction::PrintFunction() const
{
	// TODO: Edit this function to filter out the positive signs if the variable is negative.

	cout << "f(x) = " << m_a << "x^4";
	
	if (m_b != 0)
		cout << " + " << m_b << "x^3";
	
	if (m_c != 0)
		cout << " + " << m_c << "x^2";

	if (m_d != 0)
		cout << " + " << m_d << "x";

	if (m_e != 0)
		cout << " + " << m_e;


	//if (m_b == 0)
	//{
	//	// do nothing
	//}
	//else
	//{
	//	// Handle b
	//	bool bIsBPos = IsPositive<double>(m_b);
	//	char BVarPlusOrMinus(' ');

	//	if (bIsBPos)
	//	{
	//		BVarPlusOrMinus = '+';
	//	}
	//	else
	//	{
	//		BVarPlusOrMinus = ' ';
	//	}

	//	cout << BVarPlusOrMinus << m_b;
	//}

	//if (m_c == 0)
	//{
	//	// do nothing new line return
	//	//cout << endl;
	//	return;
	//}
	//else
	//{
	//	// Handle c
	//	bool bIsCPos = IsPositive<double>(m_c);
	//	char CVarPlusOrMinus(' ');

	//	if (bIsCPos)
	//	{
	//		CVarPlusOrMinus = '+';
	//	}
	//	else
	//	{
	//		CVarPlusOrMinus = ' ';
	//	}

	//	cout << CVarPlusOrMinus << m_c;
	//}
}
