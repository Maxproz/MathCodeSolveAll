#pragma once

#include <utility>

std::pair<double, double> OutputDecimalAsFraction(double input);
void PrintDecimalAsFraction(double input);

typedef std::pair<int, int> Combination;

double GetProbabilityFromCombinations(const Combination& Numerator,
										const Combination& Denominator);

// You negate when you are reversing the probability to get the chance something WONT happen
double GetProbabilityFromCombinations(const Combination& Numerator,
	const Combination& Denominator,
	const bool& bShouldNegate = false);


//double GetProbabilityFromCombinations(const Combination& RealNumerator, const Combination& Numerator, const Combination& SecondNum,
//	const Combination& Denominator);
