#include "CubicFunction.h"

#include "Derivative.h"
#include "MiscMathEquations.h"

#include <algorithm>
#include <functional>
#include <iostream>

using std::cout;
using std::endl;

void AutoSetCubicDerivativeFunction(CubicFunction& InFunc)
{
	//auto Vars = InFunc.GetABCD();
	//auto a = std::get<0>(Vars);
	//auto b = std::get<1>(Vars);
	//auto c = std::get<2>(Vars);
	//auto d = std::get<3>(Vars);

	//InFunc.PrintFunction();

	//CubicFunction CubicCopy(a, b, c, d);
	Derivative<CubicFunction, QuadraticFunction> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());

}

void CubicFunction::FindCriticalPoints()
{
	QuadraticFunction DerivativeFunc = GetDerivativeFunction();
	//DerivativeFunc.PrintFunction();
	std::vector<double> CriticalPoints = DerivativeFunc.GetAllZerosVec();

	std::sort(CriticalPoints.begin(), CriticalPoints.end(), std::less<double>());

	// if f'(c) is undefined its a critical point of f
	// AKA (This Derivative is defined for all real numbers) so this doesnt matter
	//if (std::isnan(CriticalValue))
	//{
	//
	//}


	int AmountOfCriticalPoints = CriticalPoints.size();
	cout << "Amount of critical Points: " << AmountOfCriticalPoints << endl;


	for (int i = 0; i < AmountOfCriticalPoints; ++i)
	{
		double CriticalPoint = CriticalPoints[i];
		//double CriticalValue = DerivativeFunc(CriticalPoint);
		double OutputOfCubicFromQuadraticZero = operator()(CriticalPoint);

		Point PointOfIntrest(CriticalPoint, OutputOfCubicFromQuadraticZero);



		if (CriticalPointIsALocalMinimum(DerivativeFunc, CriticalPoint))
		{
			m_LocalMinimumPoints.push_back(PointOfIntrest);
			m_AllPotentialExtrenums.push_back(PointOfIntrest);
		}

		if (CriticalPointIsALocalMaximum(DerivativeFunc, CriticalPoint))
		{
			// goes from positive to negative (local max)
			m_LocalMaximumPoints.push_back(PointOfIntrest);
			m_AllPotentialExtrenums.push_back(PointOfIntrest);
		}

		if (CriticalPointIsNeitherALocalMaxOrLocalMin(DerivativeFunc, CriticalPoint))
		{
			m_CriticalPointsNotLocalMaxOrMins.push_back(PointOfIntrest);
		}


		m_AllCriticalPoints.push_back(CriticalPoint);
		//double CriticalValueForNonDerivativeFunction = operator()(CriticalPoint);
		//m_AllPotentialExtrenums.push_back(PointOfIntrest);

	}


	SetIncreasingDecreasingIntervals();
	PrintRelatedCriticalPointData();


}




void CubicFunction::SetDefaultDomainInterval()
{
	// Cubic functions have all real number domain/ranges
	m_DomainInterval = std::make_tuple(NEGINFINITY, INFINITY, IntervalType::IT_OPEN);
}

void CubicFunction::SetDefaultRangeInterval()
{
	// Cubic functions have all real number domain/ranges
	m_RangeInterval = std::make_tuple(NEGINFINITY, INFINITY, IntervalType::IT_OPEN);
}

void CubicFunction::SetIncreasingDecreasingIntervals()
{
	if (m_bHasDomainBeenRestricted == false)
	{
		QuadraticFunction DerivativeFunc = GetDerivativeFunction();
		std::vector<double> CriticalPoints = GetAllCriticalPoints();
		std::sort(CriticalPoints.begin(), CriticalPoints.end(), std::less<double>());

		// If the domain has not been restricted its (-inf,criticalpoint[0]) U (criticalpoint[0], criticalpoint[1]) U (criticalpoint[1], inf)

		// The deriviative of a cubic is a quadratic which should always have 2 critical points/zeros

		//If the domain has not been restricted we just pick a point to the left of the first critical point and test it, since it will be negative infinity to that point,
		if (IsPositiveToTheLeftOfCriticalPoint(DerivativeFunc, CriticalPoints[0]))
		{
			m_IncreasingIntervals.push_back(Point(NEGINFINITY, CriticalPoints[0]));
		}
		if (IsNegativeToTheLeftOfCriticalPoint(DerivativeFunc, CriticalPoints[0]))
		{
			m_DecreasingIntervals.push_back(Point(NEGINFINITY, CriticalPoints[0]));
		}

		double TestValueBetweenCriticalsInterval = (CriticalPoints[0] + CriticalPoints[1]) / 2.0;
		if (IsPositiveToTheLeftOfCriticalPoint(DerivativeFunc, TestValueBetweenCriticalsInterval))
		{
			m_IncreasingIntervals.push_back(Point(CriticalPoints[0], CriticalPoints[1]));
		}
		if (IsNegativeToTheLeftOfCriticalPoint(DerivativeFunc, TestValueBetweenCriticalsInterval))
		{
			m_DecreasingIntervals.push_back(Point(CriticalPoints[0], CriticalPoints[1]));
		}



		if (IsPositiveToTheRightOfCriticalPoint(DerivativeFunc, CriticalPoints[1]))
		{
			m_IncreasingIntervals.push_back(Point(CriticalPoints[1], INFINITY));
		}
		if (IsNegativeToTheRightOfCriticalPoint(DerivativeFunc, CriticalPoints[1]))
		{
			m_DecreasingIntervals.push_back(Point(CriticalPoints[1], INFINITY));
		}

	


	}
	else
	{
		throw std::logic_error("Not prepared for this logic in SetIncreasingDecreasingIntervals()");


	}




}

void CubicFunction::FindGlobalExtremums()
{
	IntervalType DomainInterval = std::get<2>(GetDomainInterval());

	// NOTE: If we have closed the interval from a normal domain or neg inf - pos inf 
	// we are guarnteed to have a max and min and they could be at those end points,
	// so we need to do a separate test for them and add them.

	// If we still have the open interval, we dont need to do that,
	// and we just test from the max and min points that we already have

	if (DomainInterval == IntervalType::IT_CLOSED)
	{
		cout << "Closed Interval" << endl;

		double StartOfCriticalPointsInterval = std::get<0>(GetDomainInterval());
		if (std::isinf(StartOfCriticalPointsInterval) * (-1))
		{
			throw std::logic_error("You cant have a closed intervals with infinitys in it");
		}
		else
		{
			double StartOfIntervalCriticalPointValue = operator()(StartOfCriticalPointsInterval);
			m_AllPotentialExtrenums.push_back(Point(StartOfCriticalPointsInterval, StartOfIntervalCriticalPointValue));
		}


		double EndOfCriticalPointsInterval = std::get<1>(GetDomainInterval());
		if (std::isinf(EndOfCriticalPointsInterval))
		{
			throw std::logic_error("You cant have a closed intervals with infinitys in it");
		}
		else
		{
			double EndOfIntervalCriticalPointValue = operator()(EndOfCriticalPointsInterval);
			m_AllPotentialExtrenums.push_back(Point(EndOfCriticalPointsInterval, EndOfIntervalCriticalPointValue));
		}

		cout << "Closed Part2" << endl;


		if (m_AllPotentialExtrenums.size() > 1)
		{
			cout << "Test1" << endl;

			cout << m_AllPotentialExtrenums.size();

			std::sort(m_AllPotentialExtrenums.begin(), m_AllPotentialExtrenums.end(), [](auto &left, auto &right)
			{
				return left.second < right.second;
			});

			cout << "Printing all potentials for Global Externums Points\n";
			for (int i = 0; i < m_AllPotentialExtrenums.size(); ++i)
			{
				cout << "(" << m_AllPotentialExtrenums[i].first << "," << m_AllPotentialExtrenums[i].second << ")" << endl;
			}

			SetAbsoluteMinimum(m_AllPotentialExtrenums.front().first);
			SetAbsoluteMaximum(m_AllPotentialExtrenums.back().first);
		}
		else if (m_AllPotentialExtrenums.size() == 1)
		{
			cout << "Test2" << endl;

			throw std::exception("Shouldnt happen1");

		}
		else
		{
			//// There is no global minimum
			//SetAbsoluteMinimum(NAN);
			throw std::exception("Shouldnt happen2");


		}


		cout << "Printing Absolute Maximum Point\n";
		cout << m_AbsoluteMaximum << std::endl;

		cout << endl;

		cout << "Printing Absolute Minimum Point\n";
		cout << m_AbsoluteMinimum << std::endl;


		cout << endl;

		return;

	}
	else if (DomainInterval == IntervalType::IT_OPEN)
	{
		cout << "Open Interval" << endl;


		if (m_LocalMinimumPoints.size() > 1)
		{
			cout << "Test1" << endl;

			cout << m_LocalMinimumPoints.size();

			std::sort(m_LocalMinimumPoints.begin(), m_LocalMinimumPoints.end(), [](auto &left, auto &right)
			{
				return left.second < right.second;
			});


			SetAbsoluteMinimum(m_LocalMinimumPoints[0].first);
		}
		else if (m_LocalMinimumPoints.size() == 1)
		{
			cout << "Test2" << endl;
			SetAbsoluteMinimum(m_LocalMinimumPoints[0].first);
		}
		else
		{
			// There is no global minimum
			SetAbsoluteMinimum(NAN);
		}


		if (m_LocalMaximumPoints.size() > 1)
		{
			//cout << m_LocalMaximumPoints.size();
			std::sort(m_LocalMaximumPoints.begin(), m_LocalMaximumPoints.end(), [](auto &left, auto &right)
			{
				return left.second < right.second;
			});

			cout << "Test3" << endl;
			SetAbsoluteMaximum(m_LocalMaximumPoints.back().first);
		}
		else if (m_LocalMaximumPoints.size() == 1)
		{
			cout << "Test4" << endl;
			SetAbsoluteMaximum(m_LocalMaximumPoints[0].first);
		}
		else
		{
			// (m_LocalMaximumPoints.size()  == 0)
			SetAbsoluteMaximum(NAN);
		}





		cout << "Printing Absolute Maximum Point\n";
		cout << m_AbsoluteMaximum << std::endl;

		cout << endl;

		cout << "Printing Absolute Minimum Point\n";
		cout << m_AbsoluteMinimum << std::endl;


		cout << endl;
	}
}

std::string CubicFunction::GetFunctionString() const
{
	std::string OutString; 
	OutString.append(std::to_string(m_a));
	OutString.append("x^3");
	OutString.append(" + ");
	OutString.append(std::to_string(m_b));
	OutString.append("x^2 + ");
	OutString.append(std::to_string(m_c));
	OutString.append("x + ");
	OutString.append(std::to_string(m_d));

	return OutString;
}


void CubicFunction::PrintFunction() const
{
	cout << m_a << "x^3";
	cout << " + " << m_b << "x^2";
	cout << " + " << m_c << "x";
	cout << " + " << m_d;

}

void CubicFunction::PrintHorizontalTangetLineXValues() const
{
	std::cout << "Printing all HorizontalTangentLine x values of the function\n";

	for (const auto& zero : m_HorizontalTangentLines)
	{
		std::cout << "x = " << zero << std::endl;
	}
	std::cout << "Done printing the HorizontalTangentLine x values";

}
