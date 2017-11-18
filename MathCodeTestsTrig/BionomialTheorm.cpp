#include "BionomialTheorm.h"

#include <exception>
#include <sstream>  
#include <vector>
#include <iostream>

uint64_t BTFactorial(uint64_t n)
{
	return (n == 1 || n == 0) ? 1 : BTFactorial(n - 1) * n;
}

uint64_t GetPossibleCombinationsBT(const int & nTotalAmount, const int & rSelectAmount)
{
	if (nTotalAmount < 0 || rSelectAmount < 0)
		throw std::exception("n and r have to be greater than or equal to 0");

	if (nTotalAmount < rSelectAmount)
		throw std::exception("n has to be greater than or equal to r");


	uint64_t Numerator = BTFactorial(nTotalAmount);
	uint64_t Denominator = BTFactorial(rSelectAmount) * (BTFactorial(nTotalAmount - rSelectAmount));
	return Numerator / Denominator;
}

// uses (x + y)^n
std::string DisplaySimpleBionomialTherom(const int & n, const bool& bIsPlus)
{
	// TODO: maybe use stringstream here somehow?
	//std::stringstream StringStream;
	//StringStream << LocalDirectrix;

	//std::string XorY;
	//std::string EqualSign;
	//double DirectrixValue;

	//StringStream >> XorY >> EqualSign >> DirectrixValue;

	std::vector<std::string> Coefficents;

	for (int i = 0; i <= n; ++i)
	{
		Coefficents.push_back(std::to_string(GetPossibleCombinationsBT(n, i)));
	}

	int CurrentXExponent = n;
	int CurrentYExponent = 0;

	std::string FullForm;
	
	// Variable used when the input variable bIsPlus is set to false
	// The default is to have it set to true
	bool bNextIsMinus = true;


	for (int i = 0; i <= n; ++i)
	{
		

		if (Coefficents[i] != "1")
			FullForm.append(Coefficents[i]);
		

		if (CurrentXExponent > 1)
		{
			FullForm.append("x^");
			FullForm.append(std::to_string(CurrentXExponent));
		}

		if (CurrentXExponent == 1)
		{
			FullForm.append("x");
		}

		if (CurrentYExponent > 1)
		{
			FullForm.append("y^");
			FullForm.append(std::to_string(CurrentYExponent));
		}

		if (CurrentYExponent == 1)
		{
			FullForm.append("y");
		}

		CurrentXExponent = CurrentXExponent - 1;
		CurrentYExponent = CurrentYExponent + 1;

		if (bIsPlus)
		{
			if (i < n)
			{
				FullForm.append(" + ");
			}
		}
		else
		{
			if (i < n)
			{
				// Here we rotate between plus and minus on a subtraction bionomial
				if (bNextIsMinus)
				{
					FullForm.append(" - ");
					bNextIsMinus = false;
				}
				else
				{
					FullForm.append(" + ");
					bNextIsMinus = true;
				}
				
			}
		}

	}

	return FullForm;
}


// n == amount of terms
// r == one less than desired term
std::string FindASingleTerm(const int& n, const int& r)
{
	// example find 10th term, r = 9

	// use r+1 = 10
	// r = 9
	std::string OutString;

	uint64_t Comb = GetPossibleCombinationsBT(n, r);
	std::string CoffecientStr = std::to_string(Comb);
	std::cout << Comb << std::endl;
	
	OutString.append(CoffecientStr);
	OutString.append("x^");

	int XPow = n - r;
	std::string XPowString = std::to_string(XPow);
	OutString.append(XPowString);

	OutString.append("y^");
	OutString.append(std::to_string(r));

	// TODO: maybe add later.
	// If there is Coefficents on X or Y i would multiply them into the coefficent string here
	// to get the final result
	return OutString;
}

