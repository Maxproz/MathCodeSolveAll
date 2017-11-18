#pragma once


#include <map>
#include <iostream>
#include <vector>
#include <numeric>


enum class FunctionFormat
{
	Multiple,
	//MultipleAndConstant,
	Power,
	SquareRoot
};

class Sigma
{
private:
	int64_t m_UpperLimit;
	int64_t m_LowerLimit;

	FunctionFormat m_FuncFormat;
	int64_t m_FunctionVariable;

public:

	Sigma(int64_t UpLim, int64_t LowLim, FunctionFormat FuncForm, int64_t Num)
		: m_UpperLimit(UpLim), m_LowerLimit(LowLim), m_FuncFormat(FuncForm),
		m_FunctionVariable(Num)
	{

	}

	int64_t GetUpperLimit() const { return m_UpperLimit; }
	int64_t GetLowerLimit() const { return m_LowerLimit; }
	FunctionFormat GetFunctionFormat() const { return m_FuncFormat; }
	int64_t GetFuncVariable() const { return m_FunctionVariable; }


};

std::vector<int> GetSigmaNums(const Sigma& InSigma);

template <typename T>
inline void printSigma(const T& CONT)
{
	for (auto& num : CONT)
	{
		std::cout << num << " ";
	}
	int sum = std::accumulate(CONT.begin(), CONT.end(), 0);
	std::cout << "= " << sum << std::endl;

}

//  sum of the first nn terms of an arithmetic series.
inline double SumArithmeticSeries(const double& FirstTerm, 
						   const double& LastTerm,
						   const uint32_t& NumberOfTerms)
{
	return (((FirstTerm + LastTerm) * NumberOfTerms) / 2.0);
}

//  formula for the general term of an arithmetic sequence to find n.
double SumArithmeticSeriesUnknownTermAmount(const double& FirstTerm,
	const double& SecondTerm,
	const double& LastTerm);

//  formula for the general term of an arithmetic sequence to find n.
double SumArithmeticSeriesGivenFirstTermAndD(const double& FirstTerm,
	const double& D,
	const double& NumOfTerms);


// Given a geometric series, find the sum of the first n terms.
double SumGeometricSeries(const double& FirstTerm,
	const double& SecondTerm,
	const int& NumberOfTerms);

// Given a geometric series, find the sum of the first n terms.
double SumGeometricSeriesWithR(const double& FirstTerm,
	const double& CommonRatioR,
	const int& NumberOfTerms);

// Determine whether the sum of the infinite series is defined.
bool IsInfiniteSeriesDefined(const double& FirstTerm,
	const double& SecondTerm,
	const double& ThirdTerm,
	double& OutCommonRatio);


// As n gets very large, r^n gets very small.
// as n increases without bound, r^n approaches 0 
//  As r^n approaches 0, 1−r^n approaches 1.

// Given an infinite geometric series, find its sum.
// Given a geometric series, find the sum of the first n terms.
double SumInfiniteGeometricSeries(const double& FirstTerm,
	const double& SecondTerm,
	const double& ThirdTerm);
