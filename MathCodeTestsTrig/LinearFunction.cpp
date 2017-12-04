#include "LinearFunction.h"

#include "QuarticFunction.h"
#include "QuadraticFunction.h"
#include "CubicFunction.h"


QuarticFunction LinearFunction::operator*(CubicFunction const & rhs) const
{
	double ThisA = m_a;
	double ThisB = m_b;

	auto Vars = rhs.GetABCD();
	auto RhsA = std::get<0>(Vars);
	auto RhsB = std::get<1>(Vars);
	auto RhsC = std::get<2>(Vars);
	auto RhsD = std::get<3>(Vars);

	double OutQuarticA = 0;
	double OutQuarticB = 0;
	double OutQuarticC = 0;
	double OutQuarticD = 0;
	double OutQuarticE = 0;

	if (ThisA != 0)
	{
		// Multiply each variable in the Cubic by it 
		if (RhsA != 0)
		{
			OutQuarticA = ThisA * RhsA;
		}

		if (RhsB != 0)
		{
			OutQuarticB = ThisA * RhsB;
		}

		if (RhsC != 0)
		{
			OutQuarticC = ThisA * RhsC;
		}

		if (RhsD != 0)
		{
			OutQuarticD = ThisA * RhsD;
		}

	}

	if (ThisB != 0)
	{
		// add logic

		if (RhsA != 0)
		{
			OutQuarticB = OutQuarticB + (ThisB * RhsA);
		}

		if (RhsB != 0)
		{
			OutQuarticC = OutQuarticB + (ThisB * RhsB);
		}

		if (RhsC != 0)
		{
			OutQuarticD = OutQuarticD + (ThisB * RhsC);
		}
		if (RhsD != 0)
		{
			OutQuarticE = OutQuarticE + (ThisB * RhsD);
		}
	}

	return QuarticFunction(
		OutQuarticA,
		OutQuarticB,
		OutQuarticC,
		OutQuarticD,
		OutQuarticE);

}

QuadraticFunction LinearFunction::operator*(LinearFunction const & rhs) const
{
	double OutQuadA(0);
	double OutQuadB(0);
	double OutQuadC(0);

	OutQuadA = m_a * rhs.m_a;
	
	double Outside = m_a * rhs.m_b;
	double Inside = m_b * rhs.m_a;

	OutQuadB = Outside + Inside;

	OutQuadC = m_b * rhs.m_b;

	return QuadraticFunction(OutQuadA, OutQuadB, OutQuadC);
}

QuadraticFunction LinearFunction::GetSquaredFunction() const
{
	// return (ax + b)(ax + b)
	double First = m_a * m_a;

	double Outside = m_a * m_b;
	double Inside = m_b * m_a;

	double QuadraticB = Outside + Inside;

	double Last = m_b * m_b;

	return QuadraticFunction(First, QuadraticB, Last);
}

void LinearFunction::PrintAllZeros() const
{
	std::cout << "Printing all zeros of the function\n";

	for (const auto& zero : m_AllZeros)
	{
		std::cout << zero << std::endl;
	}
	std::cout << "Done printing all zeros\n";

}

std::string LinearFunction::GetFunctionString() const
{
	std::string OutString;
	OutString.append(std::to_string(m_a));
	OutString.append("x + ");
	OutString.append(std::to_string(m_b));

	return OutString;
}
