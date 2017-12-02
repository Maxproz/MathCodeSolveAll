#pragma once

#ifndef MISCMATHEQUATIONS_H
#define MISCMATHEQUATIONS_H


#include <utility>
#include <string>


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


void PrintAllFilesInDirectory(const std::string& path);


#endif