#pragma once


#ifndef DERIVATIVE_H
#define DERIVATIVE_H



#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "RationalFunction.h"
#include "QuarticFunction.h"


#include "ConstantFunction.h"

#include "PowerFunction.h"

#include "CubicFunction.h"


#include "TrigonometricFunction.h"

//class TrigometricFunction;

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
// TODO: Figure out a way to create classes for functions that are the result of calculations between polynomials and trigometric functions
// Maybe make a function that first prints the variables out and then work from there.



template <typename InFunction, typename OutFunction>
class Derivative
{
private:
	InFunction m_InFunction;
	OutFunction m_OutFunction;


	// Evaluate a quadratic function derivative get a linear function
	double EvaluateFunctionDerivative(const ConstantFunction& InFunction);
	// Evaluate a quadratic function derivative get a linear function
	LinearFunction EvaluateFunctionDerivative(const QuadraticFunction& InFunction);
	ConstantFunction EvaluateFunctionDerivative(const LinearFunction& InFunction);
	QuadraticFunction EvaluateFunctionDerivative(const CubicFunction& InFunction);
	RootFunction<-2> EvaluateFunctionDerivative(const RootFunction<2>& InFunction);
	
	//
	//TrigometricFunction<MPCOS, 1> EvaluateFunctionDerivative(const TrigometricFunction<MPSIN, 1>& InFunction);
	//TrigometricFunction<MPNEGSIN, 1> EvaluateFunctionDerivative(const TrigometricFunction<MPCOS, 1>& InFunction);
	//TrigometricFunction<MPSEC, 2> EvaluateFunctionDerivative(const TrigometricFunction<MPTAN, 1>& InFunction);
	//TrigometricFunction<MPNEGCSC, 2> EvaluateFunctionDerivative(const TrigometricFunction<MPCOT, 1>& InFunction);
	//TrigometricFunction<MPSECTAN, 1> EvaluateFunctionDerivative(const TrigometricFunction<MPSEC, 1>& InFunction);
	//TrigometricFunction<MPNEGCSCCOT, 1> EvaluateFunctionDerivative(const TrigometricFunction<MPCSC, 1>& InFunction);

	//template<TrigometricFunction<MPSIN, 1>, TrigometricFunction<MPCOS, 1>>
	 inline MPCOS<1> EvaluateFunctionDerivative(const MPSIN<1>& InFunction)
	{
		// TODO: which variables do I grab for these transfers? most online examples only show generic function
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);

		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPCOS<1> OutFunc(a, b, c, d);

		return OutFunc;
	}

	//template<TrigometricFunction<MPCOS, 1>, TrigometricFunction<MPNEGSIN, 1>>
	inline MPNEGSIN<1> EvaluateFunctionDerivative(const MPCOS<1>& InFunction)
	{
		// TODO: which variables do I grab for these transfers? most online examples only show generic function
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);


		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPNEGSIN<1> OutFunc(a, b, c, d);

		return OutFunc;
	}

	//template<TrigometricFunction<MPTAN, 1>, TrigometricFunction<MPSEC, 2>>
	inline MPSEC<2> EvaluateFunctionDerivative(const MPTAN<1>& InFunction)
	{
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);


		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPSEC<2> OutFunc(a, b, c, d);

		return OutFunc;
	}

	//template<TrigometricFunction<MPCOT, 1>, TrigometricFunction<MPNEGCSC, 2>>
	inline MPNEGCSC<2> EvaluateFunctionDerivative(const MPCOT<1>& InFunction)
	{
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);


		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPNEGCSC<2> OutFunc(a, b, c, d);

		return OutFunc;
	}

	////template<TrigometricFunction<MPSEC, 1>, TrigometricFunction<MPSECTAN, 1>>
	inline MPSECTAN<1> EvaluateFunctionDerivative(const MPSEC<1>& InFunction)
	{
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);


		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPSECTAN<1> OutFunc(a, b, c, d);

		return OutFunc;
	}

	////template<TrigometricFunction<MPCSC, 1>, TrigometricFunction<MPNEGCSCCOT, 1>>
	inline MPNEGCSCCOT<1> EvaluateFunctionDerivative(const MPCSC<1>& InFunction)
	{
		auto AllVars = InFunction.GetABCD();
		double a = std::get<0>(AllVars);
		double b = std::get<1>(AllVars);
		double c = std::get<2>(AllVars);
		double d = std::get<3>(AllVars);


		a = a * b;
		b = b;
		c = c;
		d = 0;

		MPNEGCSCCOT<1> OutFunc(a, b, c, d);

		return OutFunc;
	}

	//template <int HighestExponent, int NumberOfTerms>
	//inline PowerFunction EvaluateFunctionDerivative(PowerFunction<HighestExponent>& InFunction)
	//{
	//	// When getting derivative here the c variable doesnt matter
	//	auto a = InFunction.GetA();
	//	auto b = InFunction.GetB();

	//	// a loses x variable
	//	a = a;

	//	// b is dropped from derivative function
	//	b = 0;

	//	ConstantFunction OutFunc(a);

	//	return OutFunc;
	//}

	// To evaluate a power function with form // y = a[k(x – d)]^n + c
	template <int Exponent>
	inline PowerFunction<Exponent - 1> EvaluateFunctionDerivative(const PowerFunction<Exponent>& InFunction)
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




public:

	explicit Derivative() = default;

	//explicit Derivative(InFunction& InFunc)
	//{
	//	m_InFunction = std::move(InFunc);

	//	m_OutFunction = EvaluateFunctionDerivative(InFunc);
	//}

	explicit Derivative(const InFunction& InFunc)
		: m_InFunction(InFunc)
	{
		// m_InFunction = std::move(InFunc);

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
	double EstimateDerivative(const int& x);

	// Overloaded because I needed a double to check exact prices with decimals during rate of change production tests 
	// see CalculusFunction header -(ShouldProductionBeIncreasedUsingRateOfChange(const QuadraticFunction& Profit, const double& PriceOfItem)
	// public function that can be called, assumes that everything went ok in the constructors with assigning the m_OutFunction template
	double EstimateDerivative(const double& x);

};

// These Two sum and difference rule functions below are garbage and cant really return anything....

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

	// TODO: How do I return a function here?
	// The problem is that I have no idea how to analyze what templates I put in without adding a ton of specific code.

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

	// TODO: How do I return a function here?
	// The problem is that I have no idea how to analyze what templates I put in without adding a ton of specific code.

	FirstDerivativeFunction.PrintFunction();
	cout << " + ";
	SecondDerivativeFunction.PrintFunction();

	return;

}

// Maxpro: Good function but really specific functionality.
// Maxpro: The operator overloading really worked exactly as intended, it was very nice.
QuarticFunction ApplyDerivativeProductRule(QuadraticFunction& FirstFunction, CubicFunction& SecondFunction);

RationalFunction<QuadraticFunction, QuadraticFunction>
ApplyDerivativeQuotientRule(QuadraticFunction& Numerator, LinearFunction& Denominator);

RationalFunction<LinearFunction, QuadraticFunction>
ApplyDerivativeQuotientRule(LinearFunction& Numerator, LinearFunction& Denominator);


void ApplyDerivativePowerRules(const double& a, const double& n, double& OutA, double& OutN);

double ApplyDerivativeConstantRule(const double& a);


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