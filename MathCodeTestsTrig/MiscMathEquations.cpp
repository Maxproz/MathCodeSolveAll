#include "stdafx.h"
#include "MiscMathEquations.h"


#include <cmath>
#include <iostream>
#include <filesystem>
#include <tuple>

using std::endl;
using std::cout;
//namespace fs = std::experimental::filesystem;

double PercentageOfAValue(const double& Percent, const double& InVal)
{
	// what is 20 percent of 60?
	// 20% is mutiplied by $60:
	double OutVal = 0.0;
	OutVal = ((Percent / 100.0) * InVal);
	return OutVal;
}

// Percentage calculation
double PercentageCalculation(const double& TotalVal, const double& PartVal)
{
	// 12 is what percent of 60?
	// 12 is divided by 60 and multiplied by 100%
	double OutVal{ 0.0 };
	OutVal = (PartVal / TotalVal) * 100.0;

	// Returns OutVal as a percentage
	return OutVal;
}

// Whole value calculation
double WholeValueWithPercent(const double& Value, const double& Percentage)
{
	// 12 is 20% of what?
	// 12 is divided by 20%
	double OutVal{ 0.0 };
	OutVal = (Value / Percentage) * 100.0;

	// returns OutVal as a whole number dollar amount.
	return OutVal;
}

// Percentage change calculation
double PercentChangeCalculation(const double& OldValue, const double& NewValue)
{
	// What is the percentage change from $40 to $50 ?
	// The difference between $50 and $40 is divided by $40 and multiplied by 100 %

	double OutVal{ 0.0 };
	OutVal = ((NewValue - OldValue) / OldValue) * 100.0;

	// return OutVal as a percent
	return OutVal;
}


double HeronsFormulaCalculateS(const double& SideA, const double& SideB,
	const double& SideC)
{
	return (SideA + SideB + SideC) / 2.0;
}

double HeronsFormulaFindArea(const double & SideA, const double & SideB, const double & SideC)
{
	const double s = HeronsFormulaCalculateS(SideA, SideB, SideC);

	double OutArea = (((s - SideA)*(s - SideB)*(s - SideC)) * s);
	OutArea = std::sqrt(OutArea);
	
	// return the area in square units.
	return OutArea;
}

// TODO: Move this to misc math functions
long mygcd(long a, long b)
{
	if (a == 0)
		return b;
	else if (b == 0)
		return a;

	if (a < b)
		return mygcd(a, b % a);
	else
		return mygcd(b, a % b);
}

std::pair<double, double> OutputDecimalAsFract(const double& input)
{
	double integral = std::floor(input);
	double frac = input - integral;

	const long precision = 10; // This is the accuracy. // was 1 million or so

	long gcd_ = mygcd(round(frac * precision), precision);

	long denominator = precision / gcd_;
	long numerator = round(frac * precision) / gcd_;

	/*
	if (integral != 0)
	{
	std::cout << integral << " + ";
	}*/
	if (input < 0)
		numerator = numerator*-1;

	//std::cout << numerator << "/" << denominator << std::endl;
	return std::pair<double, double>(numerator, denominator);
}

void PrintDecimalAsFraction(double input)
{
	double integral = std::floor(input);
	double frac = input - integral;

	const long precision = 100; // This is the accuracy. // was 1 million or so

	long gcd_ = mygcd(round(frac * precision), precision);

	long denominator = precision / gcd_;
	long numerator = round(frac * precision) / gcd_;

	/*
	if (integral != 0)
	{
	std::cout << integral << " + ";
	}*/

	std::cout << numerator << "/" << denominator << std::endl;
}

// example format path "C:/Users/Maxpro/Downloads"
void PrintAllFilesInDirectory(const std::string& path)
{
	for (auto & p : std::experimental::filesystem::directory_iterator(path))
		std::cout << p << std::endl;
}

void PrintInterval(const Interval& InInterval)
{
	float a = std::get<0>(InInterval);
	float b = std::get<1>(InInterval);
	IntervalType Type = std::get<2>(InInterval);

	float IntervalStart = a;
	float IntervalEnd = b;

	switch (Type)
	{
		case IntervalType::IT_LEFT_CLOSED:
		{
			cout << "[" << IntervalStart << "," << IntervalEnd <<  ")";

			return;
		}
		//case IntervalType::IT_LEFT_OPEN:
		//{
		//	cout << "(" << IntervalStart << "," << IntervalEnd << ")";

		//	return;
		//}
		case IntervalType::IT_RIGHT_CLOSED:
		{
			cout << "(" << IntervalStart << "," << IntervalEnd << "]";

			return;
		}
		//case IntervalType::IT_RIGHT_OPEN:
		//{
		//	cout << "(" << IntervalStart << "," << IntervalEnd << ")";

		//	return;
		//}
		case IntervalType::IT_CLOSED:
		{
			cout << "[" << IntervalStart << "," << IntervalEnd << "]";

			return;
		}
		case IntervalType::IT_OPEN:
		{
			cout << "(" << IntervalStart << "," << IntervalEnd << ")";

			return;
		}
	}

}
