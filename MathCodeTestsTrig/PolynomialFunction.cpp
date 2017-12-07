#include "PolynomialFunction.h"


#include <iostream>



void PolynomialFunction::PrintEndBehaviours() const
{
	std::cout << "as x goes to " << GetEndBehaviourNegDir().first << " f(x) goes to " << GetEndBehaviourNegDir().second << std::endl;
	std::cout << "as x goes to " << GetEndBehaviourPosDir().first << " f(x) goes to " << GetEndBehaviourPosDir().second << std::endl;
}

void PolynomialFunction::PrintDomain() const
{
	std::cout << "Domain: ";
	PrintInterval(GetDomainInterval());
}

void PolynomialFunction::PrintRange() const
{
	std::cout << "Range: ";
	PrintInterval(GetRangeInterval());
}

void PolynomialFunction::SetPolynomialEndBehaviours()
{
	// As x goes to the negative direction is the first interval.
	//Interval AsXGoesToNegativeDirectionInterval = m_EndBehaviorIntervals.first;
	//float AsXGoesNegStartInterval = 0; 
	//float AsXGoesNegEndInterval = 0; // std::get<1>(AsXGoesToNegativeDirectionInterval);
	//IntervalType AsXGoesNegIntervalType = IntervalType::IT_UNASSIGNED; // std::get<2>(AsXGoesToNegativeDirectionInterval);

	//// As x goes to the positive direction is the Second interval.
	////Interval AsXGoesToPositiveDirectionInterval = m_EndBehaviorIntervals.second;
	//float AsXGoesPosStartInterval = 0;// std::get<0>(AsXGoesToPositiveDirectionInterval);
	//float AsXGoesPosEndInterval = 0; // std::get<1>(AsXGoesToPositiveDirectionInterval);
	//IntervalType AsXGoesPosIntervalType = IntervalType::IT_UNASSIGNED;// std::get<2>(AsXGoesToPositiveDirectionInterval);

	bool bIsLeadingCoefficentPositive = (m_LeadingCoefficent > 0);
	bool bIsFunctionDegreeEven = isEven(m_Degree);

	if (bIsLeadingCoefficentPositive)
	{
		// Leading coefficent is positive
		if (bIsFunctionDegreeEven)
		{
			//// leading coefficent positive and function degree even
			//float AsXGoesNegDirectionToStartInterval = 0; 
			//float FofXGoesToNegDir = 0; 

			m_XGoesNegDir = std::make_pair(NEGINFINITY, INFINITY);
			m_XGoesPosDir = std::make_pair(INFINITY, INFINITY);

			//IntervalType AsXGoesNegIntervalType = IntervalType::IT_OPEN;

			//float AsXGoesPosDirectionToEndInterval = 0;
			//float FofXGoesPosDir = 0;
			//IntervalType AsXGoesPosIntervalType = IntervalType::IT_UNASSIGNED;

		}
		else
		{
			// leading coefficent positive and function degree odd
			m_XGoesNegDir = std::make_pair(NEGINFINITY, NEGINFINITY);
			m_XGoesPosDir = std::make_pair(INFINITY, INFINITY);

		}
	}
	else
	{
		// Leading Coefficent is negative
		if (bIsFunctionDegreeEven)
		{
			// leading coefficent negative and function degree even

			m_XGoesNegDir = std::make_pair(NEGINFINITY, NEGINFINITY);
			m_XGoesPosDir = std::make_pair(INFINITY, NEGINFINITY);
		}
		else
		{
			// leading coefficent negative and function degree odd

			m_XGoesNegDir = std::make_pair(NEGINFINITY, INFINITY);
			m_XGoesPosDir = std::make_pair(INFINITY, NEGINFINITY);
		}
	}
}