#include "CountingPrinciples.h"



double CPFactorial(double n)
{
	return (n == 1 || n == 0) ? 1 : CPFactorial(n - 1) * n;
}

int GetNumberOfSubsetsInSet(const double& nDistinctObjects)
{
	return std::pow(2, nDistinctObjects); 
}


// FORMULA FOR PERMUTATIONS OF N DISTINCT OBJECTS
int GetPossiblePermutations(const double& nTotalAmount, const double& rSelectAmount)
{
	double Numerator = CPFactorial(nTotalAmount);
	double Denominator = CPFactorial(nTotalAmount - rSelectAmount);

	return Numerator / Denominator;
}

int GetPossibleCombinations(const double & nTotalAmount, const double & rSelectAmount)
{

	double Numerator = CPFactorial(nTotalAmount);
	double Denominator = CPFactorial(rSelectAmount) * (CPFactorial(nTotalAmount - rSelectAmount));

	return Numerator / Denominator;
}

int GetNumberOfPermutations(const double & nTotalAmount, const double & DuplicateAmountFirst, const double & DuplicateAmountSecond)
{
	double Numerator = CPFactorial(nTotalAmount);
	double Denominator = CPFactorial(DuplicateAmountFirst) * CPFactorial(DuplicateAmountSecond);

	return Numerator / Denominator;
}

int GetNumberOfPermutations(const double & nTotalAmount, const double & DuplicateAmountFirst, const double & DuplicateAmountSecond, const double & DuplicateAmountThird)
{
	double Numerator = CPFactorial(nTotalAmount);
	double Denominator = CPFactorial(DuplicateAmountFirst) * CPFactorial(DuplicateAmountSecond) * CPFactorial(DuplicateAmountThird);

	return Numerator / Denominator;

}
