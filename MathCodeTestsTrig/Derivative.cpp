#include "Derivative.h"


#include "Limit.h"

#include "RootFunction.h"



// Evaluate a quadratic function derivative get a linear function
template<typename InFunction, typename OutFunction>
double Derivative<InFunction, OutFunction>::EvaluateFunctionDerivative(const ConstantFunction& InFunction)
{
	return 0;
}

// Evaluate a quadratic function derivative get a linear function
template<typename InFunction, typename OutFunction>
LinearFunction Derivative<InFunction, OutFunction>::EvaluateFunctionDerivative(const QuadraticFunction& InFunction)
{
	// When getting derivative here the c variable doesnt matter
	auto ABC = InFunction.GetABC();

	double a = std::get<0>(ABC);
	double b = std::get<1>(ABC);

	// need to change for higher powered functions
	a = 2 * a;

	// loses x variable
	b = b;

	LinearFunction OutFunc(a, b);

	return OutFunc;
}

template<typename InFunction, typename OutFunction>
ConstantFunction Derivative<InFunction, OutFunction>::EvaluateFunctionDerivative(const LinearFunction& InFunction)
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

template<typename InFunction, typename OutFunction>
QuadraticFunction Derivative<InFunction, OutFunction>::EvaluateFunctionDerivative(const CubicFunction& InFunction)
{
	// When getting derivative here the c variable doesnt matter
	auto Vars = InFunction.GetABCD();
	double a = std::get<0>(Vars);
	double b = std::get<1>(Vars);
	double c = std::get<2>(Vars);
	double d = std::get<3>(Vars);

	double OutA(0);
	double OutAN(0);
	ApplyDerivativePowerRules(a, 3, OutA, OutAN);

	double OutB(0);
	double OutBN(0);
	ApplyDerivativePowerRules(b, 2, OutB, OutBN);

	double OutC(0);
	double OutCN(0);
	ApplyDerivativePowerRules(c, 1, OutC, OutCN);

	// taking the derivative of a cubic function makes the d constant equal zero
	d = 0;

	return QuadraticFunction(OutA, OutB, OutC);
}


template<typename InFunction, typename OutFunction>
RootFunction<-2> Derivative<InFunction, OutFunction>::EvaluateFunctionDerivative(const RootFunction<2>& InFunction)
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

	RootFunction<-2> OutFunc(a, b, c);

	return OutFunc;
}



// public function that can be called, assumes that everything went ok in the constructors with assigning the m_OutFunction template
template<typename InFunction, typename OutFunction>
double Derivative<InFunction, OutFunction>::EstimateDerivative(const int& x)
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
template<typename InFunction, typename OutFunction>
double Derivative<InFunction, OutFunction>::EstimateDerivative(const double& x)
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


// TODO: this was a failure to create a foundational way to do the product rule with a trigomeric and a cubic
//QuarticFunction ApplyDerivativeProductRule(const CubicFunction& FirstFunction, const TrigometricFunction<MPSIN>& SecondFunction)
//{
//	Derivative<CubicFunction, QuadraticFunction> FirstDerivative(FirstFunction);
//	QuadraticFunction FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();
//	QuadraticTrigometric<MPSIN> FirstPart = FirstDerivativeFunction * SecondFunction;
//
//
//	Derivative<TrigometricFunction<MPSIN>, TrigometricFunction<MPCOS>> SecondDerivative(SecondFunction);
//	TrigometricFunction<MPCOS> SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();
//	QuarticFunction SecondPart = SecondDerivativeFunction * FirstFunction;
//
//
//	QuarticFunction OutFunction = (FirstPart + SecondPart);
//	return OutFunction;
//}

QuarticFunction ApplyDerivativeProductRule(QuadraticFunction& FirstFunction, CubicFunction& SecondFunction)
{
	Derivative<QuadraticFunction, LinearFunction> FirstDerivative(FirstFunction);
	LinearFunction FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();
	QuarticFunction FirstPart = FirstDerivativeFunction * SecondFunction;


	Derivative<CubicFunction, QuadraticFunction> SecondDerivative(SecondFunction);
	QuadraticFunction SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();
	QuarticFunction SecondPart = SecondDerivativeFunction * FirstFunction;


	QuarticFunction OutFunction = (FirstPart + SecondPart);
	return OutFunction;
}

RationalFunction<QuadraticFunction, QuadraticFunction>
ApplyDerivativeQuotientRule(QuadraticFunction& Numerator, LinearFunction& Denominator)
{
	Derivative<QuadraticFunction, LinearFunction> FirstDerivative(Numerator);
	LinearFunction FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();
	QuadraticFunction NewNumeratorFirstPart = FirstDerivativeFunction * Denominator;

	Derivative<LinearFunction, ConstantFunction> SecondDerivative(Denominator);
	ConstantFunction SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();
	QuadraticFunction NewNumeratorSecondPart = SecondDerivativeFunction * Numerator;


	QuadraticFunction FullNewNumerator = (NewNumeratorFirstPart - NewNumeratorSecondPart);
	QuadraticFunction FullNewDenominator = (Denominator.GetSquaredFunction());

	return RationalFunction<QuadraticFunction, QuadraticFunction>(FullNewNumerator, FullNewDenominator);
}

RationalFunction<LinearFunction, QuadraticFunction>
ApplyDerivativeQuotientRule(LinearFunction& Numerator, LinearFunction& Denominator)
{
	Derivative<LinearFunction, ConstantFunction> FirstDerivative(Numerator);
	ConstantFunction FirstDerivativeFunction = FirstDerivative.GetDerivativeFunction();
	LinearFunction NewNumeratorFirstPart = FirstDerivativeFunction * Denominator;

	Derivative<LinearFunction, ConstantFunction> SecondDerivative(Denominator);
	ConstantFunction SecondDerivativeFunction = SecondDerivative.GetDerivativeFunction();
	LinearFunction NewNumeratorSecondPart = SecondDerivativeFunction * Numerator;


	LinearFunction FullNewNumerator = (NewNumeratorFirstPart - NewNumeratorSecondPart);
	QuadraticFunction FullNewDenominator = (Denominator.GetSquaredFunction());

	return RationalFunction<LinearFunction, QuadraticFunction>(FullNewNumerator, FullNewDenominator);
}

void ApplyDerivativePowerRules(const double& a, const double& n, double& OutA, double& OutN)
{
	OutA = a * n;
	OutN = n - 1;

	return;
}

double ApplyDerivativeConstantRule(const double& a)
{
	return 0;
}


