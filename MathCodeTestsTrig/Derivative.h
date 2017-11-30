#pragma once


#ifndef DERIVATIVE_H
#define DERIVATIVE_H


#include "Limit.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "ConstantFunction.h"
#include "PowerFunction.h"


// NOTES:
// Let f(x) be a function and "a" be in its domain. If f(x) is differentiable at a, then f is continuous at a.
// We conclude that in order to be differentiable at a point, a function must be “smooth” at that point.
// We saw in f(x) = cuberoot(x). A function fails to be differentiable at a point where there is a vertical tangent line.
// a function may fail to be differentiable at a point in more complicated ways as well.

// TODO: Go into more in-depth reserach and understanding so you can make better use of continuity, and differentiability.
// TODO: Add functionality for taking the derivative of a cubic function x^3
// TODO: Add functionality to handle the power rule/power functions above cubic
// TODO: Add functionality for the product rule
// TODO: Add functionality for The Quotient Rule
// TODO: Figure out how I will handle derivatives for like secx = secxtanx
// TODO: Setup a function for the constant multiple rule 
// TODO: my solution for the sum and muiltiple rule seems garbage at the moment doesn't it? not sure to tired.


template <typename InFunction, typename OutFunction>
class Derivative
{
private:
	InFunction m_InFunction;
	OutFunction m_OutFunction;


	// Evaluate a quadratic function derivative get a linear function
	inline double EvaluateFunctionDerivative(ConstantFunction& InFunction)
	{
		return 0;
	}

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

	template <int Exponent>
	inline PowerFunction<Exponent - 1> EvaluateFunctionDerivative(PowerFunction<Exponent>& InFunction)
	{
		auto GetAllInputVariables = InFunction.GetNAKDC();
		double n, a, k, d, c;
		n = std::get<0>(GetAllInputVariables);
		a = std::get<1>(GetAllInputVariables);
		k = std::get<2>(GetAllInputVariables);
		d = std::get<3>(GetAllInputVariables);
		c = std::get<4>(GetAllInputVariables);

		// take out the constant a
		double OutSideAConstant = a;

		// Now we have to take the derivative of InsideFuncDerivative: u =  k(x-d) 

		double FirstFuncA = k;
		double FirstFuncB = k * d;
		LinearFunction InsideFunc(FirstFuncA, FirstFuncB);
		Derivative<LinearFunction, ConstantFunction> InsideFuncDerivative(InsideFunc);
		ConstantFunction InsideDerivativeFunc = InsideFuncDerivative.GetDerivativeFunction();

		// And the derivative of OutsideFuncDerivative  f = 1u^n 
		
		double SecondFuncA = 1 * Exponent;
		double SecondFuncNewExponent = Exponent - 1;
		// SecondFuncA * u * InsideDerivativeFunc.m_b();

		PowerFunction<Exponent - 1> OutsideDerivativeFunc(SecondFuncA, 1, 0, 0);

		// after that I should have

		// OutSideAConstant * OutsideDerivativeFunc * InsideDerivativeFunc
		// Substitude u for k(x-d)
		OutsideDerivativeFunc = PowerFunction<Exponent - 1>(SecondFuncA, k, d, 0);
		// evaluate and simplify
		double NewK = std::pow(std::get<2>(OutsideDerivativeFunc.GetNAKDC()), Exponent - 1);
		double NewAK = NewK * SecondFuncA;
		double NewAKInsideConst = NewAK * InsideDerivativeFunc.GetB();
		double FullNewFuncConst = NewAKInsideConst * OutSideAConstant;

		PowerFunction<Exponent - 1> OutFunc(FullNewFuncConst, 1, d, 0);

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

	inline RootFunction<-2> EvaluateFunctionDerivative(RootFunction<2>& InFunction)
	{
		// TODO: which variables do I grab for these transfers? most online examples only show generic function
		auto AllVars = InFunction.GetNABC();
		double N = std::get<0>(AllVars);
		double a = std::get<1>(AllVars);
		double b = std::get<2>(AllVars);
		double c = std::get<3>(AllVars);
		
		const double RootExponent = 1.0 / N;
		
		a = a * RootExponent;
		b = b;
		c = 0;

		//double ExponentAfterSubtraction = RootExponent - 1.0;
		// dont need i guess but I could convert to fraction and add the -2 that way.

		RootFunction<-2> OutFunc(a,b,c);

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

// Let f(x)and g(x) be differentiable functions and k be a constant. Then each of the following equations holds.
template <typename F, typename FPrime, 
			typename G, typename GPrime/*,
				typename OutFunc*/>
inline void /*OutFunc*/ ApplyDerivativeDifferenceRule(F& FirstFunction, G& SecondFunction)
{
	// for j(x)=f(x)+g(x),j′(x)=f′(x)+g′(x).
	Derivative<F, FPrime> FirstDerivative(FirstFunction);
	FPrime FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();

	Derivative<G, GPrime> SecondDerivative(SecondFunction);
	GPrime SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();

	FirstDerivativeFunction.PrintFunction();
	cout << " - ";
	SecondDerivativeFunction.PrintFunction();

	return;

}

template <typename F, typename FPrime,
	typename G, typename GPrime/*,
							   typename OutFunc*/>
	inline void /*OutFunc*/ ApplyDerivativeSumRule(F& FirstFunction, G& SecondFunction)
{
	// for j(x)=f(x)+g(x),j′(x)=f′(x)+g′(x).
	Derivative<F, FPrime> FirstDerivative(FirstFunction);
	FPrime FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();

	Derivative<G, GPrime> SecondDerivative(SecondFunction);
	GPrime SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();

	FirstDerivativeFunction.PrintFunction();
	cout << " + ";
	SecondDerivativeFunction.PrintFunction();

	return;

}



template <typename FirstFuncForm, typename SecondFuncForm, typename ThirdFuncForm>
inline ThirdFuncForm GetSecondDerivativeFunction(FirstFuncForm& InFunction)
{
	Derivative<FirstFuncForm, SecondFuncForm> FirstDerivative(InFunction);
	SecondFuncForm FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();

	Derivative<SecondFuncForm, ThirdFuncForm> SecondDerivative(FirstDerivativeFunction);
	ThirdFuncForm SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();

	return SecondDerivativeFunction;
}

//template <typename InFunc, typename SecondDerivativeFunc>
//inline SecondDerivativeFunc GetSecondDerivativeFunction(InFunc& InFunction)
//{
//	// TODO: Fill out later
//	throw std::exception("Shouldn't call this template");
//
//	SecondDerivativeFunc OutFunc;
//	return OutFunc;
//}
//
//template <>
//inline ConstantFunction GetSecondDerivativeFunction<QuadraticFunction, ConstantFunction>(QuadraticFunction& InFunction)
//{
//	Derivative<QuadraticFunction, LinearFunction> FirstDerivative(InFunction);
//	LinearFunction FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();
//	Derivative<LinearFunction, ConstantFunction> SecondDerivative(FirstDerivativeFunction);
//	ConstantFunction SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();
//
//	return SecondDerivativeFunction;
//}


#endif