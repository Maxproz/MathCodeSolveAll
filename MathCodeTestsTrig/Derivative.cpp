#include "Derivative.h"



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


