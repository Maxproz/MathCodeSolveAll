#include "PolynomialFunction.h"


#include <iostream>


using std::cout;
using std::endl;


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

void PolynomialFunction::PrintIncreasingDecreasingIntervals() const
{
	std::cout << "Printing Increasing intervals of function:\n";
	for (const auto& Intervals : m_IncreasingIntervals)
	{
		std::cout << "from (" << Intervals.first << "," << Intervals.second << ")" << " the function is increasing\n";
	} std::cout << std::endl;


	std::cout << "Printing Decreasing intervals of function:\n";
	for (const auto& Intervals : m_DecreasingIntervals)
	{
		std::cout << "from (" << Intervals.first << "," << Intervals.second << ")" << " the function is decreasing\n";
	} std::cout << std::endl;


}

void PolynomialFunction::PrintRelatedCriticalPointData()
{
	cout << "Printing all Critical Points\n";
	for (int i = 0; i < m_AllCriticalPoints.size(); ++i)
	{
		cout << "x = " << m_AllCriticalPoints[i] << endl;
	}
	cout << endl;

	cout << "Printing Critical Points that are not a localmax or localmin \n";
	for (int i = 0; i < m_CriticalPointsNotLocalMaxOrMins.size(); ++i)
	{
		cout << "x = " << m_CriticalPointsNotLocalMaxOrMins[i].first << endl;
	}
	if (m_CriticalPointsNotLocalMaxOrMins.size() < 1)
	{
		cout << "None" << std::endl;
	}
	cout << endl;



	cout << "Printing all LocalMax Points\n";
	for (int i = 0; i < m_LocalMaximumPoints.size(); ++i)
	{
		cout << "(" << m_LocalMaximumPoints[i].first << "," << m_LocalMaximumPoints[i].second << ")" << endl;
	}
	cout << endl;

	cout << "Printing all LocalMin Points\n";
	for (int i = 0; i < m_LocalMinimumPoints.size(); ++i)
	{
		cout << "(" << m_LocalMinimumPoints[i].first << "," << m_LocalMinimumPoints[i].second << ")" << endl;
	}
	cout << endl;
}

void PolynomialFunction::SetPolynomialEndBehaviours()
{

	m_bIsLeadingCoefficentPositive = (m_LeadingCoefficent > 0);
	m_bIsEvenFunction = isEven(m_Degree);

	if (m_bIsLeadingCoefficentPositive)
	{
		// Leading coefficent is positive
		if (m_bIsEvenFunction)
		{
			m_XGoesNegDir = std::make_pair(NEGINFINITY, INFINITY);
			m_XGoesPosDir = std::make_pair(INFINITY, INFINITY);

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
		if (m_bIsEvenFunction)
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