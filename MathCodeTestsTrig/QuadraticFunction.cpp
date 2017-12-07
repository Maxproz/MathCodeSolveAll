#include "QuadraticFunction.h"

#include "Derivative.h"
#include "QuarticFunction.h"
#include "LinearFunction.h"
#include "MiscMathEquations.h"
#include "TrigonometricFunction.h"
#include "LinearFunction.h"

#include <algorithm>
#include <functional>

QuadraticFunction::QuadraticFunction(const std::string& FuncForm)
{
	std::istringstream iss(FuncForm);

	double a = 0;
	char FirstParathensis;
	char x;
	char PlusOrMinusOpOne;
	double h;
	char PowerOp;
	char TwoChar;
	char SecondParathensis;
	char PlusOrMinusOpTwo;
	double k;


	

	// TODO: Add check for an implicit + k missing

	int MaybeA = iss.peek();
	if (MaybeA == EOF) return;
	if (!isdigit(MaybeA))
	{
		// Assume that we are dealing with an implicit 1 for the a variable
		a = 1;
		iss >> FirstParathensis >> x >> PlusOrMinusOpOne >> h >> PowerOp >> TwoChar >> SecondParathensis >>
			PlusOrMinusOpTwo >> k;
	}
	else
	{
		// otherwise read the function like normal
		iss >> a >> FirstParathensis >> x >> PlusOrMinusOpOne >> h >> PowerOp >> TwoChar >> SecondParathensis >>
			PlusOrMinusOpTwo >> k;
	}

	// TODO: remove debug code
	std::cout << a << std::endl;
	std::cout << FirstParathensis << std::endl;
	std::cout << x << std::endl;
	std::cout << PlusOrMinusOpOne << std::endl;
	std::cout << h << std::endl;
	std::cout << PowerOp << std::endl;
	std::cout << TwoChar << std::endl;
	std::cout << SecondParathensis << std::endl;
	std::cout << PlusOrMinusOpTwo << std::endl;
	std::cout << k << std::endl;

	double TermInOne, TermInTwo, TermInThree;

	if (PlusOrMinusOpOne == '-')
	{
		if (PlusOrMinusOpTwo == '+')
		{
			// y = a*(x-h)^2 + k
			TermInOne = a;
			TermInTwo = (((h) * (-1)) * 2);
			TermInThree = (std::pow(h, 2) * (a)) + k;

			m_a = TermInOne;
			m_b = TermInTwo;
			m_c = TermInThree;

			if (a == 0)
				throw std::exception("a cannot == 0 for quadratic func initalization (vertex form constructor)");

		}
	}


	// Qudratic degree is 2 and also even 
	SetDegree(2);
	SetLeadingCoefficent(m_a);
	SetPolynomialEndBehaviours();
	SetIsEvenFunction(true);
	SetParabolaOpeningDirection(a);
	AutoSetVertexFromVertexForm(h, k);

	m_bIsFunctionVertexForm = true;

	// Normal polynomials are continuous 
	SetIsContinuousFunction(true);

	AutoSetHowManyRealZeroVariables();

	SetZerosQuadraticFormula(*this);

	AutoSetDerivativeFunction(*this);

}


void QuadraticFunction::AutoSetHowManyRealZeroVariables()
{
	const double FourAC = 4 * m_a * m_c;
	const double BSquaredMinusFourAC = std::sqrt(std::pow(m_b, 2) - FourAC);

	if (BSquaredMinusFourAC > 0)
	{
		// there are two real numbers that satisfy the quadratic equation
		m_AmountOfRealNumberZeros = 2;

	}
	else if (BSquaredMinusFourAC < 0)
	{
		// no real numbers satisfy the quadratic equation.
		m_AmountOfRealNumberZeros = 0;

	}
	else
	{
		// == 0
		// This formula tells us there is only one solution, and it is a real number.
		m_AmountOfRealNumberZeros = 1;
	}
}

QuarticFunction QuadraticFunction::operator*(QuadraticFunction const & rhs) const
{
	double QuarticA(0);
	double QuarticB(0);
	double QuarticC(0);
	double QuarticD(0);
	double QuarticE(0);

	if ((m_bIsFunctionGeneralForm && m_b == 0) &&  rhs.m_bIsFunctionGeneralForm && rhs.m_b == 0)
	{
		QuarticA = m_a * rhs.m_a;
		QuarticC = m_a * rhs.m_c;

		QuarticC = QuarticC + (m_c * rhs.m_a);
		QuarticE = m_c * rhs.m_c;

	}
	else
	{
		throw std::logic_error("You need to finish your quadratic function multiplication operator for returning quartics");
	}

	return QuarticFunction(QuarticA, QuarticB, QuarticC, QuarticD, QuarticE);
}

QuadraticFunction QuadraticFunction::operator-(QuadraticFunction const & rhs) const
{
	double OutA(0);
	double OutB(0);
	double OutC(0);

	OutA = m_a - rhs.m_a;
	OutB = m_b - rhs.m_b;
	OutC = m_c - rhs.m_c;

	return QuadraticFunction(OutA, OutB, OutC);
}


void QuadraticFunction::SetTheMaxMinValue(double InNum)
{
	if (m_a < 0)
	{
		// Opens downward (mole hill) max
		m_MaxValueAtXIsEqualTo = InNum;
	}
	else if (m_a < 0)
	{
		m_MinValueAtXIsEqualTo = InNum;
	}
	else
	{
		// this shoudlnt be possible for a quadratic we should
		// have already thrown an exception in the constructor.
	}

}

void QuadraticFunction::PrintFunctionEndBehavior() const
{
	if (m_EndBehavior == EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToPosInfinity)
	{
		std::cout << "As X goes to positive or negative infinity f(x) goes to positive infinity\n";
	}
	else if (m_EndBehavior == EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToNegInfinity)
	{
		std::cout << "As X goes to positive or negative infinity f(x) goes to negative infinity\n";
	}
	else
	{
		// not implemented...

	}

}

void QuadraticFunction::PrintParabolaOpensDirection() const
{
	if (m_ParabolaOpens == ParabolaOpen::UP)
	{
		std::cout << "The parabola opens up\n";
	}
	else if (m_ParabolaOpens == ParabolaOpen::DOWN)
	{
		std::cout << "The parabola opens down\n";
	}
	else
	{
		// not implemented...

	}
}

inline void QuadraticFunction::SetParabolaOpeningDirection(const double & LeadingCoefficent)
{
	if (LeadingCoefficent > 0)
	{
		// m_EndBehavior = EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToPosInfinity;
		m_ParabolaOpens = ParabolaOpen::UP;
	}
	else
	{
		// a < 0 at this point 
		// m_EndBehavior = EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToNegInfinity;
		m_ParabolaOpens = ParabolaOpen::DOWN;
	}
}



void QuadraticFunction::AutoSetDomainInterval()
{
	m_DomainInterval = std::make_tuple(NEGINFINITY, INFINITY, IntervalType::IT_OPEN);

	// quadratic functions have all real number domain/ranges
	// really dont need this variable in this class anymore but I can change it around later
	m_Domain = Domain::NegInfinityToPosInfinity;

}

void QuadraticFunction::AutoSetRangeInterval()
{
	switch (m_ParabolaOpens)
	{
		case ParabolaOpen::UP:
		{
			// Range is all points >= the y cord of vertex
			m_RangeInterval = std::make_tuple(m_VertexY, INFINITY, IntervalType::IT_LEFT_CLOSED);

			return;

			// old unused - 
			//m_Range = Range::NegInfinityToPosInfinity;
		}
		case ParabolaOpen::DOWN:
		{
			// Range is all points <= the y cord of vertex
			m_RangeInterval = std::make_tuple(NEGINFINITY, m_VertexY, IntervalType::IT_RIGHT_CLOSED);

			return;

			// old unused - 
			//m_Range = Range::NegInfinityToPosInfinity;
		}
	}
}

void QuadraticFunction::AutoSetVertexFromVertexForm(const double& h, const double& k)
{
	m_VertexX = h;
	m_VertexY = k;

	return;
}

void QuadraticFunction::AutoSetVertexFromGeneralForm()
{
	double NumXVertex = -1 * m_b;
	double DenominatorXVertex = 2 * m_a;

	double XVertex = NumXVertex / DenominatorXVertex;

	// Plug the XVertex into the general form equation operator
	double YVertex = operator()(XVertex);

	m_VertexX = XVertex;
	m_VertexY = YVertex;

	return;
}

void QuadraticFunction::PrintAllZeros() const
{
	std::cout << "Printing all zeros of the function\n";

	for (const auto& zero : m_AllZeros)
	{
		std::cout << zero << std::endl;
	}
	std::cout << "Done printing all zeros\n";

}

//void QuadraticFunction::PrintAllRealNumberZeros() const
//{
//	std::cout << "Printing all real number zeros of the function\n";
//
//	for (const auto& zero : m_RealNumberZeros)
//	{
//		std::cout << zero << std::endl;
//	}
//	std::cout << "Done printing all real number zeros\n";
//
//
//}

void QuadraticFunction::PrintNumberOfRealNumberSoltions() const
{
std::cout << "The function has " << m_AmountOfRealNumberZeros << " real number solutions\n";
}

void QuadraticFunction::PrintBasicFunctionInfo() const
{
	std::cout << "Starting to print basic function info\n";

	PrintFunctionEndBehavior();
	PrintParabolaOpensDirection();
	std::cout << std::endl;
	PrintAllZeros();
	std::cout << std::endl;
	PrintNumberOfRealNumberSoltions();
	//PrintAllRealNumberZeros();
	std::cout << std::endl;

	std::cout << "End of data output\n";
}

void QuadraticFunction::PrintFunction() const
{
	//cout << "f(x) = " << m_a << "x^2";
	cout << m_a << "x^2";

	if (m_b == 0)
	{
		// do nothing
	}
	else
	{
		// Handle b
		bool bIsBPos = IsPositive<double>(m_b);
		char BVarPlusOrMinus(' ');

		if (bIsBPos)
		{
			BVarPlusOrMinus = ' + ';
		}
		else
		{
			BVarPlusOrMinus = ' ';
		}

		cout << BVarPlusOrMinus << m_b << "x";
	}

	if (m_c == 0)
	{
		// do nothing new line return
		//cout << endl;
		return;
	}
	else
	{
		// Handle c
		bool bIsCPos = IsPositive<double>(m_c);
		std::string CVarPlusOrMinus(" ");

		if (bIsCPos)
		{
			CVarPlusOrMinus = " + ";
		}
		else
		{
			CVarPlusOrMinus = " ";
		}

		cout << CVarPlusOrMinus << m_c;
	}
}

std::string QuadraticFunction::GetFunctionString() const
{
	std::string OutString;
	OutString.append(std::to_string(m_a));
	OutString.append("x^2");
	OutString.append(" + ");
	OutString.append(std::to_string(m_b));
	OutString.append("x + ");
	OutString.append(std::to_string(m_c));

	return OutString;
}

void QuadraticFunction::FindCriticalPoints()
{
	// if the interval is not bounded or closed, there is no guarntee the function will have global externums

	LinearFunction DerivativeFunc = GetDerivativeFunction();
	std::vector<double> CriticalPoints = DerivativeFunc.GetAllZerosVec();

	int AmountOfCriticalPoints = CriticalPoints.size();
	cout << "Amount of critical Points: " << AmountOfCriticalPoints << endl;


	for (int i = 0; i < AmountOfCriticalPoints; ++i)
	{
		double CriticalPoint = CriticalPoints[i];
		m_CriticalPoints.push_back(CriticalPoint);
		double CriticalValue = DerivativeFunc(CriticalPoint);

		Point PointOfIntrest(CriticalPoint, CriticalValue);

		bool bGoesFromNegativeToPositiveAroundCriticalPoint = ((DerivativeFunc(CriticalPoint - 0.1) < 0) && (DerivativeFunc(CriticalPoint + 0.1) > 0));
		
		if (bGoesFromNegativeToPositiveAroundCriticalPoint)
		{
			// Local min point
			m_LocalMinimumPoints.push_back(PointOfIntrest);

		}
		else
		{
			// goes from positive to negative (local max)
			m_LocalMaximumPoints.push_back(PointOfIntrest);
		}
	}


	// One can distinguish whether a critical point is a local maximum or local minimum by using the first derivative test, 
	// second derivative test, or higher-order derivative test, given sufficient differentiability.
	std::sort(m_CriticalPoints.begin(), m_CriticalPoints.end(), std::less<double>());

	cout << "Printing all Critical Points\n";
	for (int i = 0; i < m_CriticalPoints.size(); ++i)
	{
		cout << m_CriticalPoints[i] << endl;
	}
	cout << endl;

	cout << "Printing all LocalMax Points\n";
	for (int i = 0; i < m_LocalMaximumPoints.size(); ++i)
	{
		cout << "(" << m_LocalMaximumPoints[i].first << "," << m_LocalMaximumPoints[i].second << ")" << endl;
	}
	cout << endl;

	cout << "Printing all LocalMin Points\n";
	for (int i = 0; i < m_LocalMinimumPoints.size(); ++i)
	{
		cout << "(" << m_LocalMinimumPoints[i].first << "," << m_LocalMinimumPoints[i].second << ")" << endl;
	}
	cout << endl;


	FindGlobalExtremums();


}


// Function honestly shaping up decently, need to fix the top part of it though.
void QuadraticFunction::FindGlobalExtremums()
{

	if (m_bHasDomainBeenRestricted)
	{
		cout << "Restricted" << endl;

		double StartOfCriticalPointsInterval = std::get<0>(GetDomainInterval());

		if (std::isinf(StartOfCriticalPointsInterval))
		{

		}
		else
		{
			double StartOfIntervalCriticalPointValue = operator()(StartOfCriticalPointsInterval);
		}


		double EndOfCriticalPointsInterval = std::get<1>(GetDomainInterval());

		if (std::isinf(EndOfCriticalPointsInterval))
		{

		}
		else
		{
			double EndOfIntervalCriticalPointValue = operator()(EndOfCriticalPointsInterval);
		}
	}
	else
	{
		cout << "Not Restricted" << endl;


		if (m_LocalMinimumPoints.size() > 1)
		{
			cout << "Test1" << endl;

			cout << m_LocalMinimumPoints.size();

			std::sort(m_LocalMinimumPoints.begin(), m_LocalMinimumPoints.end(), [](auto &left, auto &right)
			{
				return left.second < right.second;
			});

	
			SetAbsoluteMinimum(m_LocalMinimumPoints[0].first);
		}
		else if (m_LocalMinimumPoints.size() == 1)
		{
			cout << "Test2" << endl;
			SetAbsoluteMinimum(m_LocalMinimumPoints[0].first);
		}
		else
		{
			// There is no global minimum
			SetAbsoluteMinimum(NAN);
		}


		if (m_LocalMaximumPoints.size() > 1)
		{
			//cout << m_LocalMaximumPoints.size();
			std::sort(m_LocalMaximumPoints.begin(), m_LocalMaximumPoints.end(), [](auto &left, auto &right)
			{
				return left.second < right.second;
			});
		
			cout << "Test3" << endl;
			SetAbsoluteMaximum(m_LocalMaximumPoints.back().first);
		}
		else if (m_LocalMaximumPoints.size() == 1)
		{
			cout << "Test4" << endl;
			SetAbsoluteMaximum(m_LocalMaximumPoints[0].first);
		}
		else
		{
			// (m_LocalMaximumPoints.size()  == 0)
			SetAbsoluteMaximum(NAN);
		}


		


		cout << "Printing Absolute Maximum Point\n";
		cout << m_AbsoluteMaximum << std::endl;

		cout << endl;

		cout << "Printing Absolute Minimum Point\n";
		cout << m_AbsoluteMinimum << std::endl;


		cout << endl;
	}



	


	//if (StartOfCriticalPointsInterval == INFINITY || StartOfCriticalPointsInterval == NEGINFINITY || EndOfCriticalPointsInterval == INFINITY || EndOfCriticalPointsInterval == NEGINFINITY)
	//	throw std::domain_error("You need to fix this behaviour max");


	// TODO IMPORTANT: Function needs work once my brain is rested

	//bool AreEndpointsDiscluded = (StartOfCriticalPointsInterval == INFINITY || StartOfCriticalPointsInterval == NEGINFINITY || EndOfCriticalPointsInterval == INFINITY || EndOfCriticalPointsInterval == NEGINFINITY);

	// TODO: setup evaluation for  the closed domain version I prepared for above.

	//if (AreEndpointsDiscluded)
	//{


	//}


}

LinearFunction QuadraticFunction::GetDerivativeFunction() const
{
	return m_DerivativeFunction;
}

void QuadraticFunction::SetDerivativeFunction(const LinearFunction& InFunc)
{
	m_DerivativeFunction = InFunc;
}



void SetZerosQuadraticFormula(QuadraticFunction& QuadraticFunc)
{
	std::tuple<double, double, double> ABC = QuadraticFunc.GetABC();

	double a = std::get<0>(ABC);
	double b = std::get<1>(ABC);
	double c = std::get<2>(ABC);


	std::vector<double> LocalVecOfZeros = GetZerosQuadraticFormula(a, b, c);
	
	// Shouldnt need to check for extranneous solutions 

	//std::vector<double> OutVecRealNumZeros;

	//std::cout << LocalVecOfZeros.size();

	//for (int i = 0; i < LocalVecOfZeros.size(); ++i)
	//{
	//	double TermOne, TermTwo, TermThree;

	//	TermOne = (std::pow(LocalVecOfZeros[i], 2)) * a;
	//	TermTwo = LocalVecOfZeros[i] * b;
	//	TermThree = c;

	//	double ZeroTest = TermOne + TermTwo + TermThree;

	//	//std::cout << TermOne << std::endl;
	//	//std::cout << TermTwo << std::endl;
	//	//std::cout << TermThree << std::endl;

	//	if (is_close_to_zero(ZeroTest))
	//	{
	//		OutVecRealNumZeros.push_back(LocalVecOfZeros[i]);
	//	}
	//}


	QuadraticFunc.SetAllZeroVec(LocalVecOfZeros);
	//QuadraticFunc.SetRealNumberZeroVec(OutVecRealNumZeros);

	// Line of symmetry is between the two zeros.
	QuadraticFunc.SetLineOfSymmetry(LocalVecOfZeros[1] / LocalVecOfZeros[0]);
	QuadraticFunc.SetTheMaxMinValue(LocalVecOfZeros[1] / LocalVecOfZeros[0]);


	// dont really need to return it but whatever...
	//return OutVecRealNumZeros;

}

std::vector<double> GetZerosQuadraticFormula(const double& a, const double& b, const double& c)
{

	std::vector<double> OutVec;

	const double NegativeB = -1 * b;
	const double FourAC = 4 * a * c;
	const double BSquaredMinusFourAC = std::sqrt(std::pow(b, 2) - FourAC);
	const double TwoA = 2 * a;

	double FirstZero = ((NegativeB + BSquaredMinusFourAC) / TwoA);
	double SecondZero = ((NegativeB - BSquaredMinusFourAC) / TwoA);

	// TODO: Separate these into real zeros from all zeros later
	OutVec.push_back(FirstZero);
	OutVec.push_back(SecondZero);


	return OutVec;
}

void AutoSetDerivativeFunction(QuadraticFunction& InFunc)
{
	//auto Vars = InFunc.GetABC();
	//auto a = std::get<0>(Vars);
	//auto b = std::get<1>(Vars);
	//auto c = std::get<2>(Vars);


	//InFunc.PrintFunction();

	//CubicFunction CubicCopy(a, b, c, d);

	Derivative<QuadraticFunction, LinearFunction> Derivative(InFunc);
	InFunc.SetDerivativeFunction(Derivative.GetDerivativeFunction());


}

