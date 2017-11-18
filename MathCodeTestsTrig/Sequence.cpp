#include "Sequence.h"

#include <cmath>

std::vector<int> GetSigmaNums(const Sigma& InSigma)
{
	std::vector<int> OutVec;

	switch (InSigma.GetFunctionFormat())
	{
		case FunctionFormat::Multiple:
		{
			// Function is in the form Num*k
			for (int i = InSigma.GetLowerLimit(); i <= InSigma.GetUpperLimit(); ++i)
			{
				OutVec.push_back(InSigma.GetFuncVariable() * i);
			}

			return OutVec;
		}
		case FunctionFormat::Power:
		{
			// Function is in the form k^Num
			for (int i = InSigma.GetLowerLimit(); i <= InSigma.GetUpperLimit(); ++i)
			{
				OutVec.push_back(std::pow(i, InSigma.GetFuncVariable()));
			}

			return OutVec;
		}
		case FunctionFormat::SquareRoot:
		{
			// Not Implemented


			return OutVec;
		}
	}
}

double SumArithmeticSeriesUnknownTermAmount(const double& FirstTerm,
	const double& SecondTerm,
	const double& LastTerm)
{

	double LHS = LastTerm;

	if (FirstTerm > 0)
	{
		LHS = LHS - FirstTerm;
	}
	else
	{
		LHS = LHS + FirstTerm;
	}

	double D = SecondTerm - FirstTerm;

	LHS = LHS / D;

	double NumberOfTerms = LHS + 1;



	return (((FirstTerm + LastTerm) * NumberOfTerms) / 2.0);

}

double SumArithmeticSeriesGivenFirstTermAndD(const double & FirstTerm, const double & D, const double & NumOfTerms)
{
	double LastTerm{ 0.0 };
	LastTerm = (FirstTerm + (D*(NumOfTerms - 1)));

	return SumArithmeticSeries(FirstTerm, LastTerm, NumOfTerms);
}

double SumGeometricSeries(const double & FirstTerm, const double & SecondTerm, const int & NumberOfTerms)
{
	double CommonRatioR = SecondTerm / FirstTerm;

	if (CommonRatioR == 1)
		throw std::exception("r != 1");

	double Numerator = FirstTerm * (1.0 - std::pow(CommonRatioR, NumberOfTerms));
	double Denominator =  1.0 - CommonRatioR;
	
	return Numerator / Denominator;
}

double SumGeometricSeriesWithR(const double & FirstTerm, const double & CommonRatioR, const int & NumberOfTerms)
{
	//double CommonRatioR = SecondTerm / FirstTerm;

	if (CommonRatioR == 1)
		throw std::exception("r != 1");

	double Numerator = FirstTerm * (1.0 - std::pow(CommonRatioR, NumberOfTerms));
	double Denominator = 1.0 - CommonRatioR;

	return Numerator / Denominator;
}

// The sum of an infinite series is defined if the series is geometric and −1<r<1.
bool IsInfiniteSeriesDefined(const double & FirstTerm, const double & SecondTerm, const double & ThirdTerm, double& OutCommonRatio)
{
	double a, b;

	a = SecondTerm / FirstTerm;
	b = ThirdTerm / SecondTerm;
	
	bool BetweenNegOneAndOne = ((a > -1) && (a < 1));

	// geometric and -1<r<1
	if (a == b && BetweenNegOneAndOne)
	{
		OutCommonRatio = a;
		return true;
	}

	return false;
}

double SumInfiniteGeometricSeries(const double & FirstTerm, const double & SecondTerm, const double & ThirdTerm)
{
	double OutSum{ 0.0 };

	double CommonRatio{ 0.0 };
	bool bIsDefined = IsInfiniteSeriesDefined(FirstTerm, SecondTerm, ThirdTerm, CommonRatio);

	// geometric and -1<r<1
	if (!bIsDefined)
		throw std::exception("not geometric or r is not -1<r<1");

	OutSum = (FirstTerm / (1 - CommonRatio));

	return OutSum;
}


