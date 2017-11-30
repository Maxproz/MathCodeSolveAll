#pragma once


#ifndef DERIVATIVE_H
#define DERIVATIVE_H


#include "Limit.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "ConstantFunction.h"

template <typename InFunction, typename OutFunction>
class Derivative
{
private:
	InFunction m_InFunction;
	OutFunction m_OutFunction;


	// Evaluate a quadratic function derivative get a linear function
	inline LinearFunction EvaluateFunctionDerivative(QuadraticFunction& InFunction)
	{
		// When getting derivative here the c variable doesnt matter
		auto AB = InFunction.GetAB();

		double a = std::get<0>(AB);
		double b = std::get<1>(AB);

		// need to change for higher powered functions
		a = 2 * a;

		// loses x variable
		b = b;

		LinearFunction OutFunc(a, b);

		return OutFunc;
	}

	inline ConstantFunction EvaluateFunctionDerivative(LinearFunction& InFunction)
	{
		// When getting derivative here the c variable doesnt matter
		auto a = InFunction.GetA();
		auto b = InFunction.GetB();

		// a loses x variable
		a = a;

		// b is dropped from derivative function
		b = 0;

		ConstantFunction OutFunc(a);

		return OutFunc;
	}

	inline TrigometricFunction<MPCOS> EvaluateFunctionDerivative(TrigometricFunction<MPSIN>& InFunction)
	{
		// TODO: which variables do I grab for these transfers? most online examples only show generic function
		auto AllVars =  InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);
		
		a = a * b;
		b = b;
		c = c;
		d = 0;

		TrigometricFunction<MPCOS> OutFunc(a,b,c,d);

		return OutFunc;
	}


public:

	explicit Derivative(InFunction& InFunc)
	{
		m_InFunction = std::move(InFunc);

		m_OutFunction = EvaluateFunctionDerivative(InFunc);
	}

	inline OutFunction GetDerivativeFunction() const { return m_OutFunction; }

	//double GetDerivative(const double& x)
	//{
	//	const unsigned int h = 0;

	//	//Limit<Function> LocalLimit(InFunc, h);
	//	double Numerator = m_Function(x + h) - m_Function(x);
	//	double Denominator = h;

	//	return Numerator / Denominator;

	//}


	// public function that can be called, assumes that everything went ok in the constructors with assigning the m_OutFunction template
	inline double EstimateDerivative(const int& x)
	{
		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, std::move(m_OutFunction)/*(Numerator, Denominator)*/, x);

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

		Limit<OutFunction> LocalLimit(m_OutFunction, x);
		return LocalLimit.GetLimitResult();

	}

	// Overloaded because I needed a double to check exact prices with decimals during rate of change production tests 
	// see CalculusFunction header -(ShouldProductionBeIncreasedUsingRateOfChange(const QuadraticFunction& Profit, const double& PriceOfItem)
	// public function that can be called, assumes that everything went ok in the constructors with assigning the m_OutFunction template
	inline double EstimateDerivative(const double& x)
	{
		std::vector<std::pair<double, double>> PosDirVec;
		std::vector<std::pair<double, double>> NegDirVec;

		RunFunctionFromPosAndNegDirections(PosDirVec, NegDirVec, std::move(m_OutFunction)/*(Numerator, Denominator)*/, x);

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

		Limit<OutFunction> LocalLimit(m_OutFunction, x);
		return LocalLimit.GetLimitResult();

	}

};



#endif