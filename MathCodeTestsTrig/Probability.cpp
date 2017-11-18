#include "Probability.h"

#include <iostream>

#include "CountingPrinciples.h"

double GetProbabilityFromCombinations(const Combination & Numerator, const Combination & Denominator)
{
	double LocalNum = GetPossibleCombinations(Numerator.first, Numerator.second);
	double LocalDenom = GetPossibleCombinations(Denominator.first, Denominator.second);


	
	return(LocalNum / LocalDenom);

}

double GetProbabilityFromCombinations(const Combination & Numerator, const Combination & Denominator,
	const bool& bShouldNegate)
{
	double LocalNum = GetPossibleCombinations(Numerator.first, Numerator.second);
	double LocalDenom = GetPossibleCombinations(Denominator.first, Denominator.second);

	if (bShouldNegate)
		return 1 - (LocalNum / LocalDenom);

	return (LocalNum / LocalDenom);

}

// TODO: no idea what I was trying to do...
//double GetProbabilityFromCombinations(const Combination& RealNumerator,
//	const Combination & Numerator, const Combination & SecondNum, const Combination & Denominator)
//{
//	double LocalFirstNum =
//		GetPossibleCombinations(RealNumerator.first, RealNumerator.second);
//
//
//	double LocalNum =
//		GetPossibleCombinations(Numerator.first, Numerator.second)
//		* 
//		GetPossibleCombinations(SecondNum.first, SecondNum.second);
//	double LocalDenom = GetPossibleCombinations(Denominator.first, Denominator.second);
//
//	
//
//	// If Mutally Exclusive just add
//	double FirstFract = LocalFirstNum / LocalDenom;
//	double SecondFract = LocalNum / LocalDenom;
//
//	std::cout << FirstFract << " " << SecondFract << std::endl;
//	PrintDecimalAsFraction(FirstFract);
//	PrintDecimalAsFraction(SecondFract);
//
//	double LocalResFirst = FirstFract + SecondFract;
//
//	double LocalResSecond = 1.0 - LocalResFirst;
//
//	return LocalResSecond;
//}
