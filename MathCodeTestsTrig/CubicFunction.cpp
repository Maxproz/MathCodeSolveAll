#include "CubicFunction.h"
#include "Derivative.h"

#include <iostream>

using std::cout;
using std::endl;

void AutoSetDerivativeFunction(CubicFunction& InFunc)
{
	auto Vars = InFunc.GetABCD();
	auto a = std::get<0>(Vars);
	auto b = std::get<1>(Vars);
	auto c = std::get<2>(Vars);
	auto d = std::get<3>(Vars);

	//InFunc.PrintFunction();

	//CubicFunction CubicCopy(a, b, c, d);
	Derivative<CubicFunction, QuadraticFunction> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());

}

void CubicFunction::PrintFunction() const
{
	cout << m_a << "x^3";
	cout << " + " << m_b << "x^2";
	cout << " + " << m_c << "x";
	cout << " + " << m_d;

}
