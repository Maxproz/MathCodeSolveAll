#pragma once

#include <utility>

using std::pair;

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

std::pair<double, double> OutputDecimalAsFract(const double& input);


// Maybe I will have a use for these later
template<typename T>
double evaluate_at(double x, const T& Function)
{
	return Function(x);
}

//typedef double (RationalFunction::* memfunptr)(const double&);