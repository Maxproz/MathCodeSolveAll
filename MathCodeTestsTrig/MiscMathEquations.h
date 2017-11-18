#pragma once



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