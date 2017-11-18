#include "stdafx.h"
#include "MiscMathEquations.h"

#include <cmath>

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
