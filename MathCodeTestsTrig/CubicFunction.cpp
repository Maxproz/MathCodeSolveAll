#include "CubicFunction.h"
#include "Derivative.h"



#include <iostream>

using std::cout;
using std::endl;

void AutoSetCubicDerivativeFunction(CubicFunction& InFunc)
{
	//auto Vars = InFunc.GetABCD();
	//auto a = std::get<0>(Vars);
	//auto b = std::get<1>(Vars);
	//auto c = std::get<2>(Vars);
	//auto d = std::get<3>(Vars);

	//InFunc.PrintFunction();

	//CubicFunction CubicCopy(a, b, c, d);
	Derivative<CubicFunction, QuadraticFunction> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());

}

std::string CubicFunction::GetFunctionString() const
{
	std::string OutString; 
	OutString.append(std::to_string(m_a));
	OutString.append("x^3");
	OutString.append(" + ");
	OutString.append(std::to_string(m_b));
	OutString.append("x^2 + ");
	OutString.append(std::to_string(m_c));
	OutString.append("x + ");
	OutString.append(std::to_string(m_d));

	return OutString;
}


void CubicFunction::PrintFunction() const
{
	cout << m_a << "x^3";
	cout << " + " << m_b << "x^2";
	cout << " + " << m_c << "x";
	cout << " + " << m_d;

}

void CubicFunction::PrintHorizontalTangetLineXValues() const
{
	std::cout << "Printing all HorizontalTangentLine x values of the function\n";

	for (const auto& zero : m_HorizontalTangentLines)
	{
		std::cout << "x = " << zero << std::endl;
	}
	std::cout << "Done printing the HorizontalTangentLine x values";

}
