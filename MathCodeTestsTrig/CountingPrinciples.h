#pragma once

#include <cmath>

double CPFactorial(double n);

// Finding the Number of Subsets of a Set
int GetNumberOfSubsetsInSet(const double& nDistinctObjects); 

// FORMULA FOR PERMUTATIONS OF N DISTINCT OBJECTS
int GetPossiblePermutations(const double& nTotalAmount, const double& rSelectAmount);

// FORMULA FOR COMBINATIONS OF N DISTINCT OBJECTS
int GetPossibleCombinations(const double& nTotalAmount, const double& rSelectAmount);

// FORMULA FOR FINDING THE NUMBER OF PERMUTATIONS OF N NON-DISTINCT OBJECTS
int GetNumberOfPermutations(const double& nTotalAmount, 
	const double& DuplicateAmountFirst,
	const double& DuplicateAmountSecond = 0);

int GetNumberOfPermutations(const double& nTotalAmount,
	const double& DuplicateAmountFirst,
	const double& DuplicateAmountSecond,
	const double& DuplicateAmountThird);


