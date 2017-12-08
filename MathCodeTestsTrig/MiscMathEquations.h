#pragma once

#ifndef MISCMATHEQUATIONS_H
#define MISCMATHEQUATIONS_H


#include <utility>
#include <string>

#include "MathConstants.h"
#include "FunctionEnums.h"

using std::pair;
using std::string;


// http://www.rapidtables.com/calc/math/Percentage_Calculator.htm

// Percentage of a value calculation
double PercentageOfAValue(const double& Percent, const double& InVal);

// Percentage calculation
double PercentageCalculation(const double& TotalVal, const double& PartVal);

// Whole value calculation
double WholeValueWithPercent(const double& Value, const double& Percentage);

// Percentage change calculation
double PercentChangeCalculation(const double& OldValue, const double& NewValue);

double HeronsFormulaCalculateS(const double& SideA, const double& SideB,
	const double& SideC);

double HeronsFormulaFindArea(const double& SideA, const double& SideB,
	const double& SideC);

// This forward declaration is needed
long mygcd(long a, long b);


void PrintDecimalAsFraction(double input);

std::pair<double, double> OutputDecimalAsFract(const double& input);

// Maybe I will have a use for this later
template<typename T>
double evaluate_at(double x, const T& Function)
{
	return Function(x);
}

// Maybe I will have a use for this later
//typedef double (RationalFunction::* memfunptr)(const double&);

template <typename T>
inline bool IsPositive(const T& Num)
{
	if (Num > 0)
	{
		return true;
	}
	else if (Num < 0)
	{
		return false;
	}
	else
	{
		// Num == 0
		throw std::logic_error("error in IsPositive Function");
	}
}

template <typename T>
inline void FlipSign(T& Num)
{
	Num = Num * (-1);

}

inline bool isEven(int n) // add to helper function?
{
	if (n % 2 == 0)
		return true;
	else
		return false;
}


void PrintAllFilesInDirectory(const std::string& path);

void PrintInterval(const Interval& InInterval);


template <typename Func, typename DerivativeFunc>
inline bool IsIncreasingOverInterval(const Func& InFunc,
	const double& ClosedIntervalStart,
	const double& ClosedIntervalEnd)
{
	if (ClosedIntervalStart == NEGINFINITY)
	{
		ClosedIntervalStart = 10 * (-1);
	}

	if (ClosedIntervalEnd == INFINITY)
	{
		ClosedIntervalEnd = 10;
	}

	double OpenIntervalStart = ClosedIntervalStart + 0.1;
	double OpenIntervalEnd = ClosedIntervalEnd - 0.1;


	for ( ; OpenIntervalStart <= OpenIntervalEnd; ;)
	{
		if (!(DerivativeFunc(OpenIntervalStart) > 0))
		{
			//std::cout << "Not Increasing over interval" << std::endl;
			return false;
		}

		OpenIntervalStart = OpenIntervalStart + 0.1;
	}

	// Is Concave up

	return true;
}

template <typename Func, typename DerivativeFunc>
inline bool IsDecreasingOverInterval(const Func& InFunc,
	const double& ClosedIntervalStart,
	const double& ClosedIntervalEnd)
{
	//double ClosedIntervalStart = std::get<0>(Func.GetDomainInterval());
	//double ClosedIntervalEnd = std::get<1>(Func.GetDomainInterval());
	DerivativeFunc FirstDerivative = InFunc.GetDerivativeFunction();

	if (ClosedIntervalStart == NEGINFINITY)
	{
		ClosedIntervalStart = 10 * (-1);
	}

	if (ClosedIntervalEnd == INFINITY)
	{
		ClosedIntervalEnd = 10;
	}

	double OpenIntervalStart = ClosedIntervalStart + 0.1;
	double OpenIntervalEnd = ClosedIntervalEnd - 0.1;


	for (; OpenIntervalStart <= OpenIntervalEnd; ;)
	{
		if (!(FirstDerivative(OpenIntervalStart) < 0))
		{
			//std::cout << "Not decreasing over interval" << std::endl;
			return false;
		}

		OpenIntervalStart = OpenIntervalStart + 0.1;
	}


	// Is Concave down

	return true;
}

template <typename Func, typename FirstDerivativeFunc, typename SecondDerivativeFunc>
inline bool IsConcaveUpOverInterval(const Func& InFunc,
	const double& ClosedIntervalStart,
	const double& ClosedIntervalEnd)
{
	FirstDerivativeFunc FirstDerivative = InFunc.GetDerivativeFunction();
	SecondDerivativeFunc SecondDerivative = FirstDerivative.GetDerivativeFunction();

	if (ClosedIntervalStart == NEGINFINITY)
	{
		ClosedIntervalStart = 10 * (-1);
	}

	if (ClosedIntervalEnd == INFINITY)
	{
		ClosedIntervalEnd = 10;
	}

	double OpenIntervalStart = ClosedIntervalStart + 0.1;
	double OpenIntervalEnd = ClosedIntervalEnd - 0.1;


	for (; OpenIntervalStart <= OpenIntervalEnd; ;)
	{
		if (!(SecondDerivative(OpenIntervalStart) > 0))
		{
			//std::cout << "Not Increasing over interval" << std::endl;
			return false;
		}

		OpenIntervalStart = OpenIntervalStart + 0.1;
	}

	// Is Concave up

	return true;
}

template <typename Func, typename FirstDerivativeFunc, typename SecondDerivativeFunc>
inline bool IsConcaveDownOverInterval(const Func& InFunc,
	const double& ClosedIntervalStart,
	const double& ClosedIntervalEnd)
{
	//double ClosedIntervalStart = std::get<0>(Func.GetDomainInterval());
	//double ClosedIntervalEnd = std::get<1>(Func.GetDomainInterval());

	FirstDerivativeFunc FirstDerivative = InFunc.GetDerivativeFunction();
	SecondDerivativeFunc SecondDerivative = FirstDerivative.GetDerivativeFunction();

	if (ClosedIntervalStart == NEGINFINITY)
	{
		ClosedIntervalStart = 10 * (-1);
	}

	if (ClosedIntervalEnd == INFINITY)
	{
		ClosedIntervalEnd = 10;
	}

	double OpenIntervalStart = ClosedIntervalStart + 0.1;
	double OpenIntervalEnd = ClosedIntervalEnd - 0.1;


	for (; OpenIntervalStart <= OpenIntervalEnd; ;)
	{
		if (!(SecondDerivative(OpenIntervalStart) < 0))
		{
			//std::cout << "Not decreasing over interval" << std::endl;
			return false;
		}

		OpenIntervalStart = OpenIntervalStart + 0.1;
	}


	// Is Concave down

	return true;
}

template <typename Func>
inline bool CriticalPointIsALocalMinimum(const Func& Function, const double& CriticalPoint)
{
	// sign changes from negative to positive at the critical point
	if (IsNegativeToTheLeftOfCriticalPoint(Function, CriticalPoint) && (IsPositiveToTheRightOfCriticalPoint(Function, CriticalPoint)))
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename Func>
inline bool CriticalPointIsALocalMaximum(const Func& Function, const double& CriticalPoint)
{
	// sign changes from negative to positive at the critical point
	if (IsPositiveToTheLeftOfCriticalPoint(Function, CriticalPoint) && (IsNegativeToTheRightOfCriticalPoint(Function, CriticalPoint)))
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename Func>
inline bool CriticalPointIsNeitherALocalMaxOrLocalMin(const Func& Function, const double& CriticalPoint)
{
	// sign changes from negative to positive at the critical point
	if (CriticalPointIsALocalMaximum(Function, CriticalPoint) == false && (CriticalPointIsALocalMinimum(Function, CriticalPoint) == false))
	{
		return true;
	}
	else
	{
		return false;
	}
}


template <typename Func>
inline bool IsPositiveToTheRightOfCriticalPoint(const Func& Function, const double& CriticalPoint)
{
	// if sign of the function evaluated to the left of c is negative
	if ((Function(CriticalPoint + 0.1) > 0.0))
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename Func>
inline bool IsPositiveToTheLeftOfCriticalPoint(const Func& Function, const double& CriticalPoint)
{
	// if sign of the function evaluated to the left of c is negative
	if ((Function(CriticalPoint - 0.1) > 0.0))
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename Func>
inline bool IsNegativeToTheLeftOfCriticalPoint(const Func& Function, const double& CriticalPoint)
{
	// if sign of the function evaluated to the left of c is negative
	if ((Function(CriticalPoint - 0.1) < 0.0))
	{
		return true;
	}
	else
	{
		return false;
	}
}


template <typename Func>
inline bool IsNegativeToTheRightOfCriticalPoint(const Func& Function, const double& CriticalPoint)
{
	// if sign of the function evaluated to the left of c is negative
	if ((Function(CriticalPoint + 0.1) < 0.0))
	{
		return true;
	}
	else
	{
		return false;
	}
}


#endif