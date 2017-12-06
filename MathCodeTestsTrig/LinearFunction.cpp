#include "LinearFunction.h"

#include "QuarticFunction.h"
#include "QuadraticFunction.h"
#include "CubicFunction.h"


void LinearFunction::FindAndStoreAllRealZeros()
{
	// form y = ax + b

	// set linear func == to 0 solve
	double LocalA = m_a;
	double LocalB = m_b;

	//bool bIsANegative = m_a < 0;
	bool bIsBNegative = m_b < 0.0;

	double LeftHandSide = LocalA; // * x

	double RightHandSide{ 0.0 };


	if (bIsBNegative)
	{
		RightHandSide = RightHandSide + (LocalB * (-1));
	}
	else
	{
		RightHandSide = RightHandSide + (LocalB * (-1));
	}


	double RealZero = RightHandSide / LocalA;

	// TODO: will this give me issues later?
	const bool bIsANumber = (!(std::isnan(RealZero)));
	if (bIsANumber)
	{
		//m_AllRealZeros.push_back(RealZero);
		m_AllZeros.push_back(RealZero);
	}

}

void LinearFunction::AutoSetXAndYIntercepts(const double & b)
{
	m_XIntercept.first = (-1 * b) / m_Slope;
	m_XIntercept.second = 0;

	m_YIntercept.first = 0;
	m_YIntercept.second = b;

}

void LinearFunction::AutoSetSlopeAndFunctionForm(const double & a)
{
	m_Slope = a;

	if (m_Slope == 0)
	{
		//m_Degree = 0;
		m_bIsConstantFunction = true;
	}
	else if (m_Slope != 0)
	{
		//m_Degree = 1;
	}
	else
	{
		throw std::exception("Something went wrong with assigning the Slope and the ConstantFunction = true; variable");
	}
}

void LinearFunction::AutoSetFunctionsLineBehaviour(const double & a)
{
	if (a > 0)
	{
		// Line rises as x increases
		m_LineBehavior = LineBehavior::Increasing;

	}
	else if (a < 0)
	{
		// line falls as x increases
		m_LineBehavior = LineBehavior::Decreasing;
	}
	else
	{
		// a == 0
		// horizontal line
		m_LineBehavior = LineBehavior::Horizontal;
	}
}

void LinearFunction::AutoCheckSetDegree(const double & a)
{
	if (a == 0)
	{
		SetDegree(0);// m_Degree = 0;

	}
	else if (a != 0)
	{
		SetDegree(1); // m_Degree = 1;
	}
	else
	{
		throw std::exception("Something went wrong with assigning the degree");
	}
}

void LinearFunction::AutoSetDomainInterval()
{
	m_DomainInterval = std::make_tuple(NEGINFINITY, INFINITY, IntervalType::IT_OPEN);
	
	// Linear functions have all real number domain/ranges
	// really dont need this variable in this class anymore but I can change it around later
	m_Domain = Domain::NegInfinityToPosInfinity;
}

void LinearFunction::AutoSetRangeInterval()
{
	m_RangeInterval = std::make_tuple(NEGINFINITY, INFINITY, IntervalType::IT_OPEN);

	// really dont need this variable in this class anymore but I can change it around later
	m_Range = Range::NegInfinityToPosInfinity;
}

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

LinearFunction LinearFunction::operator-(const LinearFunction& rhs) const
{
	double OutA(0);
	double OutB(0);

	OutA = m_a - rhs.m_a;
	OutB = m_b - rhs.m_b;

	return LinearFunction(OutA, OutB);
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

void LinearFunction::PrintFunction() const
{
	//cout << "f(x) = " << m_a << "x";
	if (m_a != 0)
		cout << m_a << "x";

	if (m_b == 0)
	{
		//cout << endl;
		return;
	}
	else
	{
		char PlusOrMinus;
		if (m_b < 0)
		{
			PlusOrMinus = ' ';
		}
		else
		{
			PlusOrMinus = '+';
		}

		cout << PlusOrMinus << m_b;
	}
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


