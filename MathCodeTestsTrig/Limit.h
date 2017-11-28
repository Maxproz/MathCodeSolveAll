#pragma once

#ifndef LIMIT_H
#define LIMIT_H

// TODO: Next I am going to push the finished refactoring changes so I will have a saved working version before I try to template this class
// TODO: I should consider trying to shorten some of the code before trying to template maybe? Worth considering

#include <vector>
#include <utility>
#include <functional>
#include <iomanip>

#include "LinearFunction.h"
#include "RootFunction.h"
#include "RationalFunction.h"
#include "PiecewiseFunction.h"
#include "QuadraticFunction.h"
#include "MiscMathEquations.h"
#include "MathConstants.h" // for NEGINFINITY
#include "ConstantFunction.h"

using std::vector;
using std::pair;
using std::function;

// The limit of a constant is a constant 
// The limit of x as x approaches a is a


// TODO: need to add functionality for: ~ Evaluating a Limit by Simplifying a Complex Fraction
// TODO: need to add functionality for: ~ Evaluating a Limit When the Limit Laws Do Not Apply 
// Example: In this case, we find the limit by performing addition and then applying one of our previous strategies. 
// - https://cnx.org/contents/i4nRcikn@2.66:-xC--8XH@6/The-Limit-Laws#fs-id1170571669713
// TODO: maybe make a limit class file later?

template <typename Function>
inline void RunFunctionFromPosAndNegDirections(
	std::vector<std::pair<double, double>>& PosDirVec,
	std::vector<std::pair<double, double>>& NegDirVec,
	Function& m_Function,
	const double& m_a)
{

	std::pair<double, double> PosPairFirst;
	PosPairFirst.first = m_a + 0.1;
	PosPairFirst.second = m_Function(m_a + 0.1);

	std::pair<double, double> PosPairSecond;
	PosPairSecond.first = m_a + 0.01;
	PosPairSecond.second = m_Function(m_a + 0.01);

	std::pair<double, double> PosPairThird;
	PosPairThird.first = m_a + 0.001;
	PosPairThird.second = m_Function(m_a + 0.001);

	std::pair<double, double> PosPairFourth;
	PosPairFourth.first = m_a + 0.0001;
	PosPairFourth.second = m_Function(m_a + 0.0001);

	PosDirVec.push_back(PosPairFirst);
	PosDirVec.push_back(PosPairSecond);
	PosDirVec.push_back(PosPairThird);
	PosDirVec.push_back(PosPairFourth);


	// neg direction
	std::pair<double, double> NegPairFirst;
	NegPairFirst.first = m_a - 0.1;
	NegPairFirst.second = m_Function(m_a - 0.1);

	std::pair<double, double> NegPairSecond;
	NegPairSecond.first = m_a - 0.01;
	NegPairSecond.second = m_Function(m_a - 0.01);

	std::pair<double, double> NegPairThird;
	NegPairThird.first = m_a - 0.001;
	NegPairThird.second = m_Function(m_a - 0.001);

	std::pair<double, double> NegPairFourth;
	NegPairFourth.first = m_a - 0.0001;
	NegPairFourth.second = m_Function(m_a - 0.0001);

	NegDirVec.push_back(NegPairFirst);
	NegDirVec.push_back(NegPairSecond);
	NegDirVec.push_back(NegPairThird);
	NegDirVec.push_back(NegPairFourth);

}

template <typename Function>
class Limit
{
private:

	//// TempVariable
	//bool m_bIsLinearFuncLimit = false;
	//LinearFunction m_LinearFunction;
	//template <typename Function>

	Function m_Function;


	double m_a;
	//std::function<double(const double&)> m_Function;

	double m_L;

	double m_LimitFromNegDir = 0;
	double m_LimitFromPosDir = 0;

	// triggers to true if a root function is detected in the numerator or denominator
	// This disables the solution finding by factor/canceling
	bool m_bTryConjSolution = false;

	//double SimpifyComplexFraction(const ComplexFraction& InComplexFract);

	double EvaluateRootFuncLimit(const RootFunction& InRootFunc);

	//inline RationalFunction SolveByConjugateMultiplication(const RootFunction& Numerator, const LinearFunction& Denominator)
	//{
	//	auto NumeratorVars = Numerator.GetNABC();

	//	auto N = std::get<0>(NumeratorVars);

	//	auto A = std::get<1>(NumeratorVars);

	//	auto B = std::get<2>(NumeratorVars);

	//	auto C = std::get<3>(NumeratorVars);

	//	LinearFunction NumeratorRes;
	//	RootFunction DenominatorRes;

	//	bool bBInputIsPositive = (B > 0);


	//	LinearFunction TopRes;
	//	RootFunction BottomRes;

	//	// find conjugate
	//	if (C < 0)
	//	{
	//		// add C in conjugate
	//		//double Conjugate = (A*(std::pow(m_a - B, (1.0 / N)))) + C;
	//		if (bBInputIsPositive)
	//		{
	//			TopRes = LinearFunction(1, B*(-1) + C*(-C));// == x + TopRes
	//			BottomRes = RootFunction(N, A, B*(1), (C*(-1))); // Bottom Res = 
	//		}
	//		else
	//		{
	//			TopRes = LinearFunction(1, B*(-1) + C*(-C));// == x + TopRes
	//			BottomRes = RootFunction(N, A, B*(1), (C*(-1))); // Bottom Res = 
	//		}

	//		std::cout << TopRes.GetA() << std::endl;
	//		std::cout << Denominator.GetA() << std::endl;
	//		std::cout << TopRes.GetB() << std::endl;
	//		std::cout << Denominator.GetB() << std::endl;

	//		if (TopRes.GetA() == Denominator.GetA() && TopRes.GetB() == Denominator.GetB())
	//		{


	//			std::cout << "TEST INSIDE CANCEL FUNCS" << std::endl;
	//			// Pretend these canceled Out
	//			NumeratorRes = LinearFunction(0, 1);


	//			DenominatorRes = BottomRes;
	//		}
	//	}
	//	else
	//	{
	//		// subtract C

	//	}

	//	return RationalFunction(NumeratorRes, DenominatorRes);

	//}

	// helper function to help evaluate a limit if its a linear function
	inline double EvaluateLinearFuncLimit(LinearFunction& InLinearFunc)
	{
		double Tempa = InLinearFunc.GetA();
		double Tempb = InLinearFunc.GetB();

		double FirstLimit = ApplyLimitConstantMultipleLaw(Tempa);
		double SecondLimit = ApplyBasicLimitRuleForConstants(Tempb);

		return FirstLimit + SecondLimit;
	}


	inline double EvaluateQuadraticFuncLimit(const QuadraticFunction& InQuadraticFunc)
	{
		// TODO: maybe change this to use limit rules later
		// Instead of doing limit rules I will just plug in the value
		return InQuadraticFunc(m_a);
	}




		//if ((NumeratorRes == 0 && DenominatorRes == 0) && m_bTryConjSolution == false)
		//{
		//	LinearFunction NewFactoredNumerator;
		//	LinearFunction NewFactoredDenominator;

		//	std::vector<double> NumeratorZeros = Numerator.GetRealNumberZerosVec();
		//	std::vector<double> DenominatorZeros;
		//	DenominatorZeros.push_back(((Denominator.GetB() * (-1)) / Denominator.GetA()));

		//	std::vector<double> NewNumerator;

		//	for (int i = 0; i < NumeratorZeros.size(); ++i)
		//	{
		//		auto result1 = std::find(std::begin(DenominatorZeros), std::end(DenominatorZeros), NumeratorZeros[i]);

		//		if (result1 != std::end(DenominatorZeros))
		//		{
		//			std::cout << "Numerator contains: " << NumeratorZeros[i] << '\n';
		//		}
		//		else
		//		{
		//			std::cout << "Numerator does not contain: " << NumeratorZeros[i] << '\n';
		//			NewNumerator.push_back(NumeratorZeros[i]);
		//		}
		//	}

		//	for (int i = 0; i < NewNumerator.size(); ++i)
		//	{
		//		std::cout << NewNumerator[i] << " ";
		//		std::cout << std::endl;
		//	}

		//	//double One = std::abs(std::floor(NewNumerator[0]));
		//	//double Two = NewNumerator[0];

		//	//std::cout << "abs " << One << std::endl;
		//	//std::cout << Two << std::endl;

		//	if (/*std::abs*/(std::floor(NewNumerator[0])) == NewNumerator[0])
		//	{
		//		std::cout << "NewNumerator[0] is whole\n";


		//		if (NewNumerator[0] == 0.0)
		//		{
		//			NewFactoredNumerator = LinearFunction(1, 0);
		//		}
		//		else
		//		{
		//			if (NewNumerator[0] < 0.0)
		//			{
		//				NewFactoredNumerator = LinearFunction(1, NewNumerator[0] * (-1));
		//			}
		//			else
		//			{
		//				NewFactoredNumerator = LinearFunction(1, NewNumerator[0]);
		//			}
		//		}
		//	}
		//	else
		//	{
		//		std::cout << "NewNumerator[0] is not whole\n";

		//		std::pair<double, double> NumeratorFract = OutputDecimalAsFract(NewNumerator[0]);

		//		std::cout << NumeratorFract.first << "/" << NumeratorFract.second << std::endl;

		//		if (NumeratorFract.first < 0)
		//		{
		//			NewFactoredNumerator = LinearFunction(NumeratorFract.second, NumeratorFract.first);
		//			//std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
		//		}
		//	}



		//	std::vector<double> NewDenominator;

		//	for (int i = 0; i < DenominatorZeros.size(); ++i)
		//	{
		//		auto result1 = std::find(std::begin(NumeratorZeros), std::end(NumeratorZeros), DenominatorZeros[i]);

		//		if (result1 != std::end(NumeratorZeros))
		//		{
		//			std::cout << "Denominator contains: " << DenominatorZeros[i] << '\n';
		//		}
		//		else
		//		{
		//			std::cout << "Denominator does not contain: " << DenominatorZeros[i] << '\n';
		//			NewDenominator.push_back(DenominatorZeros[i]);
		//		}
		//	}

		//	if (NewDenominator.empty())
		//	{
		//		NewFactoredNumerator.PrintLinearFunctionInfo();

		//		double FinalFactoredNumerator = EvaluateLinearFuncLimit(NewFactoredNumerator);

		//		std::cout << "FF Num: " << FinalFactoredNumerator << std::endl;

		//		return FinalFactoredNumerator;

		//	}
		//	else
		//	{

		//		for (int i = 0; i < NewDenominator.size(); ++i)
		//		{
		//			std::cout << NewDenominator[i] << " ";
		//			std::cout << std::endl;
		//		}

		//		if (std::abs(std::floor(NewDenominator[0])) == NewDenominator[0])
		//		{
		//			std::cout << "NewDenominator[0] is whole\n";

		//			if (NewDenominator[0] == 0.0)
		//			{
		//				NewFactoredDenominator = LinearFunction(1, 0);
		//			}
		//			else
		//			{
		//				if (NewDenominator[0] < 0.0)
		//				{
		//					NewFactoredDenominator = LinearFunction(1, NewDenominator[0]);
		//				}
		//				else
		//				{
		//					NewFactoredDenominator = LinearFunction(1, NewDenominator[0] * (-1));
		//				}
		//			}

		//		}
		//		else
		//		{
		//			std::cout << "NewDenominator[0] is not whole\n";

		//			std::pair<double, double> DenominatorFract = OutputDecimalAsFract(NewDenominator[0]);
		//			if (DenominatorFract.first < 0)
		//			{
		//				std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
		//				NewFactoredDenominator = LinearFunction(DenominatorFract.second, DenominatorFract.first*(-1));

		//			}
		//		}

		//		NewFactoredNumerator.PrintLinearFunctionInfo();
		//		NewFactoredDenominator.PrintLinearFunctionInfo();

		//		double FinalFactoredNumerator = EvaluateLinearFuncLimit(NewFactoredNumerator);
		//		double FinalFactoredDenominator = EvaluateLinearFuncLimit(NewFactoredDenominator);

		//		std::cout << "FF Num: " << FinalFactoredNumerator << std::endl;
		//		std::cout << "FF Den: " << FinalFactoredDenominator << std::endl;

		//		return FinalFactoredNumerator / FinalFactoredDenominator;

		//	}



			//// TODO: remove debug code
			//std::cout << "Numerator:\t " << NumeratorRes << std::endl;
			//std::cout << "Denominator:\t " << DenominatorRes << std::endl;


			//return NumeratorRes / DenominatorRes;
		


	inline double EvaluateFuncLimit(RationalFunction<LinearFunction, LinearFunction>& InFunction)
	{

		LinearFunction Numerator = InFunction.GetNumeratorFunction();
		LinearFunction Denominator = InFunction.GetDenominatorFunction();

		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, InFunction(Numerator, Denominator), m_a);



		std::cout << "Evaluating Limit: Please Wait...\n";

		for (auto & num : PosDirVec)
		{

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		double TopRes = PosDirVec[3].second;
		double BottomRes = NegDirVec[3].second;
		/*std::cout << TopRes << std::endl;
		std::cout << BottomRes << std::endl;*/


		std::string TopResStr = std::to_string(TopRes);
		// TODO: remove debug code later
		std::cout << TopResStr << std::endl;

		std::string BottomResStr = std::to_string(BottomRes);
		// TODO: remove debug code later
		std::cout << BottomResStr << std::endl;

		double LocalPosRes = std::stod(TopResStr) * 100;

		double LocalNegRes = std::stod(BottomResStr) * 100;

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;



		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;


		// TODO: I need a better way to check for infinity here its not detecting a small increase to a limit at 1

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > (PosDirVec[1].second * 9))
			{
				bIsPosDirPosInfinity = true;
			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < (NegDirVec[1].second * 9))
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < (PosDirVec[1].second * 9))
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second >(NegDirVec[1].second * 9))
			{
				bIsNegDirPosInfinity = true;
			}
		}


		double LimitFromPosDir = 0;
		double LimitFromNegDir = 0;




		if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
		{
			std::cout << "We have a working limit result! 1\n";


			if (bIsPosDirPosInfinity)
			{
				if (bIsNegDirPosInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << INFINITY << std::endl;
					return INFINITY;
				}
			}

			if (bIsPosDirNegInfinity)
			{
				if (bIsNegDirNegInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << NEGINFINITY << std::endl;
					return NEGINFINITY;
				}
			}

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			//	// return either
			//return std::ceil(LocalPosRes) / 100;
		}

		// if you have two results that are returning negatives
		// you need to do some sort of swap with floor / ceil
		if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
		{
			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(LocalPosRes) / 100;

		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;

			LimitFromPosDir = INFINITY;
		}
		else if (bIsPosDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;

			LimitFromPosDir = NEGINFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << std::nearbyint(PosDirVec[3].second) << std::endl;


			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);
		}

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;

			LimitFromNegDir = NEGINFINITY;
		}
		else if (bIsNegDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;

			LimitFromNegDir = INFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " ";
			std::cout << (std::nearbyint(NegDirVec[3].second)) << std::endl;

			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);

		}


		m_LimitFromPosDir = LimitFromPosDir;
		m_LimitFromNegDir = LimitFromNegDir;

		std::cout << std::endl;

		// TODO: What do I return here, does it matter?
		//return std::floor(LocalPosRes) / 100;
		return InFunction.GetLastCalculatedRes();

	}
	
	inline double EvaluateFuncLimit(RationalFunction<QuadraticFunction, LinearFunction>& InFunction)
	{
		//double NumeratorRes{ 0 };
		//double DenominatorRes{ 0 };

		//// Make sure the denominator is not equal to 0

		//QuadraticFunction Numerator = InFunction.GetNumeratorFunction();
		//LinearFunction Denominator = InFunction.GetDenominatorFunction();

		//NumeratorRes = Numerator(m_a);
		//DenominatorRes = Denominator(m_a);

		//if (DenominatorRes == 0)
		//	throw std::domain_error("Denominator function cannot == 0");

		
		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, std::move(InFunction)/*(Numerator, Denominator)*/, m_a);

		std::cout << "Evaluating Limit: Please Wait...\n\n";

		std::cout << "From Positive Direction\n";
		for (auto & num : PosDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		std::cout << "From Negative Direction\n";
		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		} 
		
		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;


		// TODO: I need a better way to check for infinity here its not detecting a small increase to a limit at 1

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > (PosDirVec[1].second * 9))
			{
				bIsPosDirPosInfinity = true;
			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < (NegDirVec[1].second * 9))
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < (PosDirVec[1].second * 9))
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second >(NegDirVec[1].second * 9))
			{
				bIsNegDirPosInfinity = true;
			}
		}


		double LimitFromPosDir = 0;
		double LimitFromNegDir = 0;

		double PositiveDirectionApproaches = PosDirVec[3].second;
		double NegativeDirectionApporoaches = NegDirVec[3].second;


		// Check to see if both directions go to infinity.
		if (bIsPosDirPosInfinity)
		{
			if (bIsNegDirPosInfinity)
			{
				std::cout << "We have a working limit result! INFINITY: \n";
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << INFINITY << std::endl;
				return INFINITY;
			}
		}

		// Check to see if both directions go to neg infinity
		if (bIsPosDirNegInfinity)
		{
			if (bIsNegDirNegInfinity)
			{
				std::cout << "We have a working limit result! NEGINFINITY: \n";
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << NEGINFINITY << std::endl;
				return NEGINFINITY;
			}
		}
		

		// Do I not need this part?
		//// Check to see if the limit from both sides is the same
		//if (std::ceil(PositiveDirectionApproaches) == std::floor(NegativeDirectionApporoaches))
		//{
		//	// x approaches -2
		//	// We are dealing with 2 positive numbers? 
		//	// -1.999, -2.111
		//	// ceil -1.9

		//	std::cout << "We have a working limit result! Ceil pos dir, floor \n";



		//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

		//	//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

		//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


		//	//	// return either
		//	//return std::ceil(LocalPosRes) / 100;
		//}


		if (std::floor(PositiveDirectionApproaches) == std::ceil(NegativeDirectionApporoaches))
		{
			// x approaches 2
			// We are dealing with 2 positive numbers? 
			// 2.111, 1.999
			// floor: 2.111 == 2;
			// ceil:  1.999 == 2;

			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(PosDirVec[3].second) /*/ 100*/ << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(PosDirVec[3].second);// / 100;

		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		// DNE~ From the postive direction

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
			LimitFromPosDir = INFINITY;
		}
		else if (bIsPosDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
			LimitFromPosDir = NEGINFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << std::nearbyint(PosDirVec[3].second) << std::endl;


			// http://en.cppreference.com/w/cpp/numeric/math/nearbyint
			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);
		}

		// DNE~ From the negative direction

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
			LimitFromNegDir = NEGINFINITY;
		}
		else if (bIsNegDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
			LimitFromNegDir = INFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " ";
			std::cout << (std::nearbyint(NegDirVec[3].second)) << std::endl;


			// http://en.cppreference.com/w/cpp/numeric/math/nearbyint
			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);

		}

		m_LimitFromPosDir = LimitFromPosDir;
		m_LimitFromNegDir = LimitFromNegDir;

		return NAN;

		//return NumeratorRes / DenominatorRes;

				
				//}

		//		break;
		//	}
		//	case PolynomialFunctionType::LINEAR:
		//	{


		//		break;
		//	}
		//	case PolynomialFunctionType::QUADRATIC:
		//	{


		//		break;
		//	}
		//	case PolynomialFunctionType::CUBIC:
		//	{


		//		break;
		//	}
		//}

		//// TODO: remove debug code
		//std::cout << "Res:\t " << Result << std::endl;


		//// TODO: For fraction functions that contain square roots. Conjugate multiplication. SET IT UP
		//if (InRationalFunc.GetNumeratorFunctionType() == PolynomialFunctionType::ROOT &&
		//	InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::LINEAR)
		//{
		//	RootFunction Numerator = InRationalFunc.GetNumeratorRoot();

		//	NumeratorRes = Numerator(m_a);

		//	//std::cout << NumeratorRes << std::endl;

		//	LinearFunction Denominator = InRationalFunc.GetDenominatorLinear();

		//	DenominatorRes = Denominator(m_a);

		//	// root function detected, disable factor/cancel solution option
		//	m_bTryConjSolution = true;


		//	if ((NumeratorRes == 0 && DenominatorRes == 0) && m_bTryConjSolution == true)
		//	{

		//		std::cout << "Testttttttttttttt" << std::endl;
		//		RationalFunction NewFactoredFunction = SolveByConjugateMultiplication(Numerator, Denominator);


		//		//NumeratorRes = NewFactoredFunction.GetNumeratorLinear()(m_a);
		//		DenominatorRes = NewFactoredFunction.GetDenominatorRoot()(m_a);




		//		if (NewFactoredFunction.GetNumeratorLinear().IsBOnlyForm())
		//		{
		//			NumeratorRes = NewFactoredFunction.GetNumeratorLinear().GetB();
		//		}
		//		//std::cout << "Test" <<  NumeratorRes << std::endl;

		//		//LinearFunction NewFactoredNumerator;
		//		//RootFunction NewFactoredDenominator;
		//	}
		//}

		//if (InRationalFunc.GetNumeratorFunctionType() == PolynomialFunctionType::QUADRATIC &&
		//	InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
		//{
		//	QuadraticFunction Numerator = InRationalFunc.GetNumeratorQuadratic();

		//	if (Numerator.IsABForm())
		//	{
		//		std::tuple<double, double> AB = Numerator.GetAB();

		//		double A = std::get<0>(AB);
		//		double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		//		double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);
		//		double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		//		double B = std::get<1>(AB);
		//		double FinalB = ApplyBasicLimitRuleForX() * B;
		//		NumeratorRes = FinalA + FinalB;
		//	}
		//	else
		//	{

		//		// a b c form
		//		// you now know the form break up the limits
		//		std::tuple<double, double, double> ABC = Numerator.GetABC();

		//		double A = std::get<0>(ABC);

		//		double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		//		double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

		//		double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		//		double B = std::get<1>(ABC);
		//		double FinalB = ApplyBasicLimitRuleForX() * B;

		//		double C = std::get<2>(ABC);
		//		double FinalC = ApplyBasicLimitRuleForConstants(C);

		//		NumeratorRes = FinalA + FinalB + FinalC;
		//	}


		//	//if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
		//	//{
		//	QuadraticFunction Denominator = InRationalFunc.GetDenominatorQuadratic();

		//	if (Denominator.IsACForm()) // TODO: Set this up next
		//	{
		//		std::tuple<double, double> AC = Denominator.GetAC();

		//		double A = std::get<0>(AC);
		//		double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		//		double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);
		//		double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		//		//std::cout << FinalA << std::endl;

		//		double C = std::get<1>(AC);
		//		double FinalC = ApplyBasicLimitRuleForConstants(C);
		//		//std::cout << FinalC << std::endl;

		//		DenominatorRes = FinalA + FinalC;
		//		//std::cout << DenominatorRes << std::endl;
		//	}
		//	else
		//	{
		//		// a b c form
		//		// you now know the form break up the limits
		//		std::tuple<double, double, double> ABC = Denominator.GetABC();

		//		double A = std::get<0>(ABC);

		//		double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		//		double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

		//		double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		//		double B = std::get<1>(ABC);
		//		double FinalB = ApplyBasicLimitRuleForX() * B;

		//		double C = std::get<2>(ABC);
		//		double FinalC = ApplyBasicLimitRuleForConstants(C);

		//		DenominatorRes = FinalA + FinalB + FinalC;
		//	}


		//	//}


		//	if ((NumeratorRes == 0 && DenominatorRes == 0) && m_bTryConjSolution == false)
		//	{
		//		LinearFunction NewFactoredNumerator;
		//		LinearFunction NewFactoredDenominator;

		//		std::vector<double> NumeratorZeros = Numerator.GetRealNumberZerosVec();
		//		std::vector<double> DenominatorZeros = Denominator.GetRealNumberZerosVec();

		//		std::vector<double> NewNumerator;

		//		for (int i = 0; i < NumeratorZeros.size(); ++i)
		//		{
		//			auto result1 = std::find(std::begin(DenominatorZeros), std::end(DenominatorZeros), NumeratorZeros[i]);

		//			if (result1 != std::end(DenominatorZeros))
		//			{
		//				std::cout << "Numerator contains: " << NumeratorZeros[i] << '\n';
		//			}
		//			else
		//			{
		//				std::cout << "Numerator does not contain: " << NumeratorZeros[i] << '\n';
		//				NewNumerator.push_back(NumeratorZeros[i]);
		//			}
		//		}

		//		for (int i = 0; i < NewNumerator.size(); ++i)
		//		{
		//			std::cout << NewNumerator[i] << " ";
		//			std::cout << std::endl;
		//		}

		//		//double One = std::abs(std::floor(NewNumerator[0]));
		//		//double Two = NewNumerator[0];

		//		//std::cout << "abs " << One << std::endl;
		//		//std::cout << Two << std::endl;

		//		if (/*std::abs*/(std::floor(NewNumerator[0])) == NewNumerator[0])
		//		{
		//			std::cout << "NewNumerator[0] is whole\n";


		//			if (NewNumerator[0] == 0.0)
		//			{
		//				NewFactoredNumerator = LinearFunction(1, 0);
		//			}
		//			else
		//			{
		//				if (NewNumerator[0] < 0.0)
		//				{
		//					NewFactoredNumerator = LinearFunction(1, NewNumerator[0] * (-1));
		//				}
		//				else
		//				{
		//					NewFactoredNumerator = LinearFunction(1, NewNumerator[0]);
		//				}
		//			}
		//		}
		//		else
		//		{
		//			std::cout << "NewNumerator[0] is not whole\n";

		//			std::pair<double, double> NumeratorFract = OutputDecimalAsFract(NewNumerator[0]);

		//			std::cout << NumeratorFract.first << "/" << NumeratorFract.second << std::endl;

		//			if (NumeratorFract.first < 0)
		//			{
		//				NewFactoredNumerator = LinearFunction(NumeratorFract.second, NumeratorFract.first);
		//				//std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
		//			}
		//		}



		//		std::vector<double> NewDenominator;

		//		for (int i = 0; i < NumeratorZeros.size(); ++i)
		//		{
		//			auto result1 = std::find(std::begin(NumeratorZeros), std::end(NumeratorZeros), DenominatorZeros[i]);

		//			if (result1 != std::end(NumeratorZeros))
		//			{
		//				std::cout << "Numerator contains: " << DenominatorZeros[i] << '\n';
		//			}
		//			else
		//			{
		//				std::cout << "Numerator does not contain: " << DenominatorZeros[i] << '\n';
		//				NewDenominator.push_back(DenominatorZeros[i]);
		//			}
		//		}

		//		for (int i = 0; i < NewDenominator.size(); ++i)
		//		{
		//			std::cout << NewDenominator[i] << " ";
		//			std::cout << std::endl;
		//		}

		//		if (std::abs(std::floor(NewDenominator[0])) == NewDenominator[0])
		//		{
		//			std::cout << "NewDenominator[0] is whole\n";

		//			if (NewDenominator[0] == 0.0)
		//			{
		//				NewFactoredDenominator = LinearFunction(1, 0);
		//			}
		//			else
		//			{
		//				if (NewDenominator[0] < 0.0)
		//				{
		//					NewFactoredDenominator = LinearFunction(1, NewDenominator[0]);
		//				}
		//				else
		//				{
		//					NewFactoredDenominator = LinearFunction(1, NewDenominator[0] * (-1));
		//				}
		//			}

		//		}
		//		else
		//		{
		//			std::cout << "NewDenominator[0] is not whole\n";

		//			std::pair<double, double> DenominatorFract = OutputDecimalAsFract(NewDenominator[0]);
		//			if (DenominatorFract.first < 0)
		//			{
		//				std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
		//				NewFactoredDenominator = LinearFunction(DenominatorFract.second, DenominatorFract.first*(-1));

		//			}
		//		}

		//		NewFactoredNumerator.PrintLinearFunctionInfo();
		//		NewFactoredDenominator.PrintLinearFunctionInfo();

		//		double FinalFactoredNumerator = EvaluateLinearFuncLimit(NewFactoredNumerator);
		//		double FinalFactoredDenominator = EvaluateLinearFuncLimit(NewFactoredDenominator);

		//		std::cout << "FF Num: " << FinalFactoredNumerator << std::endl;
		//		std::cout << "FF Den: " << FinalFactoredDenominator << std::endl;

		//		return FinalFactoredNumerator / FinalFactoredDenominator;

		//	}
		//}


		//if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::CUBIC)
		//{
		//	CubicFunction Denominator = InRationalFunc.GetDenominatorCubic();

		//	if (Denominator.GetIsFuncInAAndDForm())
		//	{

		//		std::pair<double, double> AAndD = Denominator.GetAAndDCubicFuncForm();

		//		double A = AAndD.first;

		//		double TempAPowerLawCubic = ApplyLimitRuleForPowers(3);
		//		double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

		//		double FinalA = TempAConstMultipleLaw*TempAPowerLawCubic;

		//		double D = AAndD.second;
		//		double TempDConstRule = ApplyBasicLimitRuleForConstants(D);

		//		double FinalD = TempDConstRule;

		//		DenominatorRes = FinalA + FinalD;
		//	}
		//}

		////if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::LINEAR)
		////{
		////	LinearFunction Denominator = InRationalFunc.GetDenominatorLinear();

		////	DenominatorRes = EvaluateLinearFuncLimit(Denominator);
		////}

		////if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
		////{
		////	QuadraticFunction Denominator = InRationalFunc.GetDenominatorQuadratic();

		////	// a b c form
		////	// you now know the form break up the limits
		////	std::tuple<double, double, double> ABC = Denominator.GetABC();

		////	double A = std::get<0>(ABC);

		////	double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		////	double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

		////	double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		////	double B = std::get<1>(ABC);
		////	double FinalB = ApplyBasicLimitRuleForX() * B;

		////	double C = std::get<2>(ABC);
		////	double FinalC = ApplyBasicLimitRuleForConstants(C);

		////	NumeratorRes = FinalA + FinalB + FinalC;
		////}

		//// TODO: remove debug code
		//std::cout << "Numerator:\t " << NumeratorRes << std::endl;
		//std::cout << "Denominator:\t " << DenominatorRes << std::endl;


		//if (NumeratorRes == 0 && DenominatorRes == 0)
		//{

		//
		//	std::pair<double, double> Fract = OutputDecimalAsFract(NumeratorRes / DenominatorRes);

		//	std::cout << Fract.first << " / " << Fract.second << std::endl;

		//}


		//return NumeratorRes / DenominatorRes;
	}

	inline double EvaluateFuncLimit(RationalFunction<ConstantFunction, QuadraticFunction>& InFunction)
	{

		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, std::move(InFunction)/*(Numerator, Denominator)*/, m_a);

		std::cout << "Evaluating Limit: Please Wait...\n\n";

		std::cout << "From Positive Direction\n";
		for (auto & num : PosDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		std::cout << "From Negative Direction\n";
		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;


		// TODO: I need a better way to check for infinity here its not detecting a small increase to a limit at 1

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > (PosDirVec[1].second * 9))
			{
				bIsPosDirPosInfinity = true;
			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < (NegDirVec[1].second * 9))
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < (PosDirVec[1].second * 9))
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second >(NegDirVec[1].second * 9))
			{
				bIsNegDirPosInfinity = true;
			}
		}


		double LimitFromPosDir = 0;
		double LimitFromNegDir = 0;

		double PositiveDirectionApproaches = PosDirVec[3].second;
		double NegativeDirectionApporoaches = NegDirVec[3].second;


		// Check to see if both directions go to infinity.
		if (bIsPosDirPosInfinity)
		{
			if (bIsNegDirPosInfinity)
			{
				std::cout << "We have a working limit result! INFINITY: \n";
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << INFINITY << std::endl;
				return INFINITY;
			}
		}

		// Check to see if both directions go to neg infinity
		if (bIsPosDirNegInfinity)
		{
			if (bIsNegDirNegInfinity)
			{
				std::cout << "We have a working limit result! NEGINFINITY: \n";
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << NEGINFINITY << std::endl;
				return NEGINFINITY;
			}
		}


		// Do I not need this part?
		//// Check to see if the limit from both sides is the same
		//if (std::ceil(PositiveDirectionApproaches) == std::floor(NegativeDirectionApporoaches))
		//{
		//	// x approaches -2
		//	// We are dealing with 2 positive numbers? 
		//	// -1.999, -2.111
		//	// ceil -1.9

		//	std::cout << "We have a working limit result! Ceil pos dir, floor \n";



		//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

		//	//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

		//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


		//	//	// return either
		//	//return std::ceil(LocalPosRes) / 100;
		//}


		if (std::floor(PositiveDirectionApproaches) == std::ceil(NegativeDirectionApporoaches))
		{
			// x approaches 2
			// We are dealing with 2 positive numbers? 
			// 2.111, 1.999
			// floor: 2.111 == 2;
			// ceil:  1.999 == 2;

			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(PosDirVec[3].second) /*/ 100*/ << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(PosDirVec[3].second);// / 100;

		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		// DNE~ From the postive direction

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
			LimitFromPosDir = INFINITY;
		}
		else if (bIsPosDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
			LimitFromPosDir = NEGINFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << std::nearbyint(PosDirVec[3].second) << std::endl;


			// http://en.cppreference.com/w/cpp/numeric/math/nearbyint
			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);
		}

		// DNE~ From the negative direction

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
			LimitFromNegDir = NEGINFINITY;
		}
		else if (bIsNegDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
			LimitFromNegDir = INFINITY;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " ";
			std::cout << (std::nearbyint(NegDirVec[3].second)) << std::endl;


			// http://en.cppreference.com/w/cpp/numeric/math/nearbyint
			LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
			LimitFromNegDir = std::nearbyint(NegDirVec[3].second);

		}

		m_LimitFromPosDir = LimitFromPosDir;
		m_LimitFromNegDir = LimitFromNegDir;

		return NAN;

	}

	inline double EvaluateCubicFunctionLimit(const CubicFunction& InCubicFunc)
	{
		if (InCubicFunc.GetIsFuncInACDForm())
		{
			std::tuple<double, double, double> ACDVariables = InCubicFunc.GetACDForm();

			double A = std::get<0>(ACDVariables);
			double C = std::get<1>(ACDVariables);
			double D = std::get<2>(ACDVariables);

			const double x = ApplyBasicLimitRuleForX();

			double PowerLawA = ApplyLimitRuleForPowers(3);
			double FinalA = PowerLawA * A;

			double FinalC = C*x;

			double FinalD = ApplyBasicLimitRuleForConstants(D);

			return FinalA + FinalC + FinalD;

		}

		return 0.0;

	}

	// TODO: add evaluations for other limit functions // logaritm/rationals etc.

	// This functions assumes you put into it the multiple of x of a linear function for the moment
	inline double ApplyLimitConstantMultipleLaw(const double& InA)
	{
		// m_a assigned during the constructor
		// This is moved outside of the limit here
		double c = InA;

		// Now inside of the limit all we have is an x
		// taking the limit of just an x is a!
		double LimitOfX = ApplyBasicLimitRuleForX();

		return c*LimitOfX;
	}

	inline double ApplyLimitRuleForPowers(const double& InN)
	{
		double XLimit = ApplyBasicLimitRuleForX();
		return std::pow(XLimit, InN);
	}

	inline double ApplyBasicLimitRuleForConstants(const double& InConstant)
	{
		double LimitOfConstant = InConstant;
		return LimitOfConstant;
	}

	inline double ApplyBasicLimitRuleForX()
	{
		return m_a;
	}

	template <typename FirstFunc, typename SecondFunc>
	inline double EvaluatePiecewiseFuncLimit(const PiecewiseFunction<typename FirstFunc, typename SecondFunc>& InFunc)
	{

		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		//PiecewiseFunction<FirstFunc, SecondFunc> Piecewise = InFunc;

		// TODO: no idea how to fix this deleted function problem
		//RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, InFunc, m_a);


		//std::vector<std::pair<double, double>> PosDirVec;
		//// pos direction
		//std::pair<double, double> PosPairFirst;
		//PosPairFirst.first = m_a + 0.1;
		//PosPairFirst.second = InFunc(m_a + 0.1);

		//std::pair<double, double> PosPairSecond;
		//PosPairSecond.first = m_a + 0.01;
		//PosPairSecond.second = InFunc(m_a + 0.01);

		//std::pair<double, double> PosPairThird;
		//PosPairThird.first = m_a + 0.001;
		//PosPairThird.second = InFunc(m_a + 0.001);

		//std::pair<double, double> PosPairFourth;
		//PosPairFourth.first = m_a + 0.0001;
		//PosPairFourth.second = InFunc(m_a + 0.0001);

		//PosDirVec.push_back(PosPairFirst);
		//PosDirVec.push_back(PosPairSecond);
		//PosDirVec.push_back(PosPairThird);
		//PosDirVec.push_back(PosPairFourth);

		//std::vector<std::pair<double, double>> NegDirVec;
		//// neg direction
		//std::pair<double, double> NegPairFirst;
		//NegPairFirst.first = m_a - 0.1;
		//NegPairFirst.second = InFunc(m_a - 0.1);

		//std::pair<double, double> NegPairSecond;
		//NegPairSecond.first = m_a - 0.01;
		//NegPairSecond.second = InFunc(m_a - 0.01);

		//std::pair<double, double> NegPairThird;
		//NegPairThird.first = m_a - 0.001;
		//NegPairThird.second = InFunc(m_a - 0.001);

		//std::pair<double, double> NegPairFourth;
		//NegPairFourth.first = m_a - 0.0001;
		//NegPairFourth.second = InFunc(m_a - 0.0001);

		//NegDirVec.push_back(NegPairFirst);
		//NegDirVec.push_back(NegPairSecond);
		//NegDirVec.push_back(NegPairThird);
		//NegDirVec.push_back(NegPairFourth);

		std::cout << "Evaluating Limit: Please Wait...\n";

		for (auto & num : PosDirVec)
		{

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		double TopRes = PosDirVec[3].second;
		double BottomRes = NegDirVec[3].second;
		/*std::cout << TopRes << std::endl;
		std::cout << BottomRes << std::endl;*/


		std::string TopResStr = std::to_string(TopRes);
		// TODO: remove debug code later
		std::cout << TopResStr << std::endl;

		std::string BottomResStr = std::to_string(BottomRes);
		// TODO: remove debug code later
		std::cout << BottomResStr << std::endl;

		double LocalPosRes = std::stod(TopResStr) * 100;

		double LocalNegRes = std::stod(BottomResStr) * 100;

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;



		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;


		// TODO: I need a better way to check for infinity here its not detecting a small increase to a limit at 1

		//if (PosDirVec[1].second > PosDirVec[0].second)
		//{
		//	if (PosDirVec[2].second > PosDirVec[1].second)
		//	{
		//		bIsPosDirPosInfinity = true;

		//	}
		//}

		//if (NegDirVec[1].second < NegDirVec[0].second)
		//{
		//	if (NegDirVec[2].second < NegDirVec[1].second)
		//	{
		//		bIsNegDirNegInfinity = true;
		//	}
		//}

		//if (PosDirVec[1].second < PosDirVec[0].second)
		//{
		//	if (PosDirVec[2].second < PosDirVec[1].second)
		//	{
		//		bIsPosDirNegInfinity = true;
		//	}
		//}
		//if (NegDirVec[1].second > NegDirVec[0].second)
		//{
		//	if (NegDirVec[2].second > NegDirVec[1].second)
		//	{
		//		bIsNegDirPosInfinity = true;
		//	}
		//}







		if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
		{
			std::cout << "We have a working limit result! 1\n";


			if (bIsPosDirPosInfinity)
			{
				if (bIsNegDirPosInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << INFINITY << std::endl;
					return INFINITY;
				}
			}

			if (bIsPosDirNegInfinity)
			{
				if (bIsNegDirNegInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << NEGINFINITY << std::endl;
					return NEGINFINITY;
				}
			}

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			//	// return either
			//return std::ceil(LocalPosRes) / 100;
		}

		// if you have two results that are returning negatives
		// you need to do some sort of swap with floor / ceil
		if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
		{
			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(LocalPosRes) / 100;
		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << std::nearbyint(PosDirVec[3].second) << std::endl;
		}

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
		}
		else
		{
			//std::fesetround(FE_TONEAREST);
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " ";
			std::cout << (std::nearbyint(NegDirVec[3].second)) << std::endl;
		}


		auto LimitFromPosDir = std::nearbyint(PosDirVec[3].second);
		auto LimitFromNegDir = std::nearbyint(NegDirVec[3].second);

		m_LimitFromPosDir = LimitFromPosDir;
		m_LimitFromNegDir = LimitFromNegDir;

		std::cout << std::endl;

		return std::floor(LocalPosRes) / 100;
	}

	template <typename FirstFunc, typename SecondFunc, int ThirdFunctionConstant>
	inline double EvaluatePiecewiseFuncLimitThreeFunctions(
		const PiecewiseFunctionThreeFunctions<typename FirstFunc, typename SecondFunc, ThirdFunctionConstant>& InFunc)
	{

		std::vector<std::pair<double, double>> PosDirVec;
		// pos direction
		std::pair<double, double> PosPairFirst;
		PosPairFirst.first = m_a + 0.1;
		PosPairFirst.second = InFunc(m_a + 0.1);

		std::pair<double, double> PosPairSecond;
		PosPairSecond.first = m_a + 0.01;
		PosPairSecond.second = InFunc(m_a + 0.01);

		std::pair<double, double> PosPairThird;
		PosPairThird.first = m_a + 0.001;
		PosPairThird.second = InFunc(m_a + 0.001);

		std::pair<double, double> PosPairFourth;
		PosPairFourth.first = m_a + 0.0001;
		PosPairFourth.second = InFunc(m_a + 0.0001);

		PosDirVec.push_back(PosPairFirst);
		PosDirVec.push_back(PosPairSecond);
		PosDirVec.push_back(PosPairThird);
		PosDirVec.push_back(PosPairFourth);

		std::vector<std::pair<double, double>> NegDirVec;
		// neg direction
		std::pair<double, double> NegPairFirst;
		NegPairFirst.first = m_a - 0.1;
		NegPairFirst.second = InFunc(m_a - 0.1);

		std::pair<double, double> NegPairSecond;
		NegPairSecond.first = m_a - 0.01;
		NegPairSecond.second = InFunc(m_a - 0.01);

		std::pair<double, double> NegPairThird;
		NegPairThird.first = m_a - 0.001;
		NegPairThird.second = InFunc(m_a - 0.001);

		std::pair<double, double> NegPairFourth;
		NegPairFourth.first = m_a - 0.0001;
		NegPairFourth.second = InFunc(m_a - 0.0001);

		NegDirVec.push_back(NegPairFirst);
		NegDirVec.push_back(NegPairSecond);
		NegDirVec.push_back(NegPairThird);
		NegDirVec.push_back(NegPairFourth);

		std::cout << "Evaluating Limit: Please Wait...\n";

		for (auto & num : PosDirVec)
		{

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		double TopRes = PosDirVec[3].second;
		double BottomRes = NegDirVec[3].second;
		/*std::cout << TopRes << std::endl;
		std::cout << BottomRes << std::endl;*/


		std::string TopResStr = std::to_string(TopRes);
		// TODO: remove debug code later
		std::cout << TopResStr << std::endl;

		std::string BottomResStr = std::to_string(BottomRes);
		// TODO: remove debug code later
		std::cout << BottomResStr << std::endl;

		double LocalPosRes = std::stod(TopResStr) * 100;

		double LocalNegRes = std::stod(BottomResStr) * 100;

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;



		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;


		// TODO: I need a better way to check for infinity here its not detecting a small increase to a limit at 1

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > PosDirVec[1].second)
			{
				bIsPosDirPosInfinity = true;

			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < NegDirVec[1].second)
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < PosDirVec[1].second)
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second > NegDirVec[1].second)
			{
				bIsNegDirPosInfinity = true;
			}
		}

		if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
		{
			std::cout << "We have a working limit result! 1\n";


			if (bIsPosDirPosInfinity)
			{
				if (bIsNegDirPosInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << INFINITY << std::endl;
					return INFINITY;
				}
			}

			if (bIsPosDirNegInfinity)
			{
				if (bIsNegDirNegInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << NEGINFINITY << std::endl;
					return NEGINFINITY;
				}
			}

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			//	// return either
			//return std::ceil(LocalPosRes) / 100;
		}

		// if you have two results that are returning negatives
		// you need to do some sort of swap with floor / ceil
		if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
		{
			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(LocalPosRes) / 100;
		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
		}

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
		}

		std::cout << std::endl;


		std::cout << std::ceil(LocalPosRes) / 100 << std::endl;
		// TODO: Temporary mesaure for this one function problem :: generalize later
		return std::ceil(LocalPosRes) / 100;
	}

	bool GetIsLinearFunctionLimit() const { return m_bIsLinearFuncLimit; }

public:

	LinearFunction GetLinearFunctionIfExists() const
	{
		if (GetIsLinearFunctionLimit() == false)
			throw std::logic_error("This limit does not hold a single linear function");

		return m_LinearFunction;
	}


	double GetLimitFromPosDir() const { return m_LimitFromPosDir; }
	double GetLimitFromNegDir() const { return m_LimitFromNegDir; }

	template <typename NumeratorFunc, typename DenominatorFunc>
	void CheckAndSetRationalFuncDiscontinunities(RationalFunction<NumeratorFunc,DenominatorFunc>& InRationalFunc, const double& x)
	{
		if (!(std::isnan(m_L)))
		{
			// is a real number? maybe not sure exactly what it defines "as a number"
			//	Since f is discontinuous at 2 and limx→2 exists, f has a removable discontinuity at x = 2.

			InRationalFunc.SetDiscontinunityPtr(DiscontinunityType::REMOVEABLE, m_L);
		}

		//std::cout << " TESTING " << m_L << std::endl;

		bool LimitFromPosDirPosOrNegInfinity = (std::isinf(m_LimitFromPosDir));
		bool LimitFromNegDirPosOrNegInfinity = (std::isinf(m_LimitFromNegDir));

		if (LimitFromPosDirPosOrNegInfinity && LimitFromNegDirPosOrNegInfinity)
			InRationalFunc.SetDiscontinunityPtr(DiscontinunityType::INFINITE, x);

	}

	void CheckAndSetPiecewiseFuncQuadLinearDiscontinunities(PiecewiseFunction<QuadraticFunction, LinearFunction>& InPiecewiseFunc, const double& x)
	{
		//std::cout << " TESTING " << m_L << std::endl;

		bool LimitFromPosDirExists = !(std::isnan(m_LimitFromPosDir));
		bool LimitFromNegDirExists = !(std::isnan(m_LimitFromNegDir));

		if (LimitFromPosDirExists && LimitFromNegDirExists)
		{
			InPiecewiseFunc.SetDiscontinunityPtr(DiscontinunityType::JUMP, x);
		}
	}




	 explicit Limit(Function& InFunction, const double& a)
		: /* m_Function(InFunction), */m_a(a)
	{
		m_Function = std::move(InFunction);

		double TestIfNan = m_Function(m_a);

		if (std::isnan(TestIfNan))
		{
			// If its not a number we got to evaluate from both sides
			m_L = EvaluateFuncLimit(InFunction);
		}
		else
		{
			// If its a real number this is our limit
			m_L = TestIfNan;
		}

		// automatically run the limit on construction
		//m_L = operator()(a);

		//std::cout << "Is NumeratorFuncType Quadratic?: ";
		//auto NumeratorFuncType = m_Function.GetNumeratorFunctionType();
		//
		//if (NumeratorFuncType == PolynomialFunctionType::QUADRATIC)
		//{
		//	std::cout << "Move Success" << endl;
		//}
		//else
		//{
		//	cout << "Fail" << endl;
		//}

		// TODO: remove debug code
		DisplayLimitResult();
	}

	
	 // TODO: Set your evaluate func limit function up so it works for all types of template instantiations

	//explicit Limit(std::function<double(const double&)> InFunc, const double& a)
	//	: m_Function(InFunc), m_a(a)
	//{

	//	// automatically run the limit on construction
	//	m_L = operator()(a);
	//}

	//explicit Limit(LinearFunction& InLinearFunc, const double& a)
	//	: m_a(a)
	//{
	//	m_L = EvaluateLinearFuncLimit(InLinearFunc);
	//	m_bIsLinearFuncLimit = true;
	//	m_LinearFunction = InLinearFunc;

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}

	//explicit Limit(const CubicFunction& InCubicFunc, const double& a)
	//	: m_a(a)
	//{
	//	m_L = EvaluateCubicFunctionLimit(InCubicFunc);

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}

	//template <typename NumeratorFunc, typename DenominatorFunc>
	//explicit Limit(const RationalFunction<NumeratorFunc, DenominatorFunc>& InRationalFunc, const double& a)
	//	: m_a(a)
	//{
	//	m_L = EvaluateRationalFuncLimit(InRationalFunc);

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}



	//explicit Limit(const QuadraticFunction& InQuadraticFunc, const double& a)
	//	: m_a(a)
	//{

	//	m_L = EvaluateQuadraticFuncLimit(InQuadraticFunc);

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}


	//explicit Limit(const ComplexFraction& InComplexFract, const double& a)
	//	: m_a(a)
	//{

	//	m_L = SimpifyComplexFraction(InComplexFract);

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}

	//explicit Limit(const RootFunction& InRootFunc, const double& a)
	//	: m_a(a)
	//{

	//	m_L = EvaluateRootFuncLimit(InRootFunc);

	//	// TODO: remove debug code
	//	DisplayLimitResult();
	//}

	//template <typename FirstFunc, typename SecondFunc>
	//explicit Limit(const PiecewiseFunction<typename FirstFunc, typename SecondFunc>& InPiecewiseFunc, const double& a)
	//	: m_a(a)
	//{

	//	m_L = EvaluatePiecewiseFuncLimit(InPiecewiseFunc);



	//	// TODO: remove debug code
	//	//DisplayLimitResult();
	//}

	//template <typename FirstFunc, typename SecondFunc, int ThirdFunctionConstant>
	//explicit Limit(const PiecewiseFunctionThreeFunctions<typename FirstFunc, typename SecondFunc, ThirdFunctionConstant>& InPiecewiseFunc, const double& a)
	//	: m_a(a)
	//{

	//	m_L = EvaluatePiecewiseFuncLimitThreeFunctions(InPiecewiseFunc);


	//	// TODO: remove debug code
	//	//DisplayLimitResult();
	//}




	// TODO: Clean this function up by adding helper functions
	double operator()(double x) const
	{

		std::vector<std::pair<double, double>> PosDirVec;
		// pos direction
		std::pair<double, double> PosPairFirst;
		PosPairFirst.first = m_a + 0.1;
		PosPairFirst.second = m_Function(m_a + 0.1);

		std::pair<double, double> PosPairSecond;
		PosPairSecond.first = m_a + 0.01;
		PosPairSecond.second = m_Function(m_a + 0.01);

		std::pair<double, double> PosPairThird;
		PosPairThird.first = m_a + 0.001;
		PosPairThird.second = m_Function(m_a + 0.001);

		std::pair<double, double> PosPairFourth;
		PosPairFourth.first = m_a + 0.0001;
		PosPairFourth.second = m_Function(m_a + 0.0001);

		PosDirVec.push_back(PosPairFirst);
		PosDirVec.push_back(PosPairSecond);
		PosDirVec.push_back(PosPairThird);
		PosDirVec.push_back(PosPairFourth);

		std::vector<std::pair<double, double>> NegDirVec;
		// neg direction
		std::pair<double, double> NegPairFirst;
		NegPairFirst.first = m_a - 0.1;
		NegPairFirst.second = m_Function(m_a - 0.1);

		std::pair<double, double> NegPairSecond;
		NegPairSecond.first = m_a - 0.01;
		NegPairSecond.second = m_Function(m_a - 0.01);

		std::pair<double, double> NegPairThird;
		NegPairThird.first = m_a - 0.001;
		NegPairThird.second = m_Function(m_a - 0.001);

		std::pair<double, double> NegPairFourth;
		NegPairFourth.first = m_a - 0.0001;
		NegPairFourth.second = m_Function(m_a - 0.0001);

		NegDirVec.push_back(NegPairFirst);
		NegDirVec.push_back(NegPairSecond);
		NegDirVec.push_back(NegPairThird);
		NegDirVec.push_back(NegPairFourth);

		std::cout << "Evaluating Limit: Please Wait...\n";

		for (auto & num : PosDirVec)
		{

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		double TopRes = PosDirVec[3].second;
		double BottomRes = NegDirVec[3].second;
		/*std::cout << TopRes << std::endl;
		std::cout << BottomRes << std::endl;*/


		std::string TopResStr = std::to_string(TopRes);
		// TODO: remove debug code later
		std::cout << TopResStr << std::endl;

		std::string BottomResStr = std::to_string(BottomRes);
		// TODO: remove debug code later
		std::cout << BottomResStr << std::endl;

		double LocalPosRes = std::stod(TopResStr) * 100;

		double LocalNegRes = std::stod(BottomResStr) * 100;

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;



		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > PosDirVec[1].second)
			{
				bIsPosDirPosInfinity = true;
			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < NegDirVec[1].second)
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < PosDirVec[1].second)
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second > NegDirVec[1].second)
			{
				bIsNegDirPosInfinity = true;
			}
		}

		if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
		{
			std::cout << "We have a working limit result! 1\n";


			if (bIsPosDirPosInfinity)
			{
				if (bIsNegDirPosInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << INFINITY << std::endl;
					return INFINITY;
				}
			}

			if (bIsPosDirNegInfinity)
			{
				if (bIsNegDirNegInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << NEGINFINITY << std::endl;
					return NEGINFINITY;
				}
			}

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			//	// return either
			//return std::ceil(LocalPosRes) / 100;
		}

		// if you have two results that are returning negatives
		// you need to do some sort of swap with floor / ceil
		if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
		{
			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(LocalPosRes) / 100;
		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
		}

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
		}

		std::cout << std::endl;

		return 0;
	}

	//void TakeOneSidedLimits()


	//inline void TakeLimit(double x)
	//{
	//	m_L = operator()(x);
	//}

	inline void DisplayLimitResult()
	{
		std::cout << "Limit of f(x) as x --> " << m_a << " is " << m_L << std::endl;
	}


	double GetLimitResult() const { return m_L; }



};


#endif

