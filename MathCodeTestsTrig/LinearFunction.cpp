#include "LinearFunction.h"

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
