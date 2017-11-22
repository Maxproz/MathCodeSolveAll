#include "CalculusFunction.h"

#include <math.h>       /* isnan, sqrt */

int RationalFunction::m_AmountOfDiscontinunitiesFound = 0;

void PrintPointSlopeForm(const double& Slope, const Point& InPoint)
{
	double X1 = InPoint.first;
	double Y1 = InPoint.second;

	std::cout << "y - " << Y1 << " = " << Slope << "(x - " << X1 << ")" << std::endl;

}

double GetSlope(const Point& FirstPoint, const Point& SecondPoint)
{
	double X1, X2, Y1, Y2;

	X1 = FirstPoint.first;
	Y1 = FirstPoint.second;
	X2 = SecondPoint.first;
	Y2 = SecondPoint.second;

	double Numerator = Y2 - Y1;
	double Denominator = X2 - X1;

	return Numerator / Denominator;
}

void PrintSlopeInterceptForm(const Point& Point, const double& Slope)
{
	double RHS = Point.first;
	double LHS = Point.second;

	if (LHS > 0)
	{
		LHS = LHS*-1;
	}

	if (RHS > 0)
	{
		RHS = RHS*-1;
	}

	RHS = Slope*RHS;


	if (LHS < 0)
	{
		RHS = RHS + LHS;
	}
	else
	{
		// LHS number is needs subtracted ( > 0)
		RHS = RHS - LHS;
	}

	double B = RHS;

	std::cout << "y = " << Slope << "x + " << B << std::endl;
}



void QuadraticFunction::AutomaticSetRealZeroVariables()
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

void QuadraticFunction::PrintAllZeros() const
{
	std::cout << "Printing all zeros of the function\n";

	for (const auto& zero : m_AllZeros)
	{
		std::cout << zero << std::endl;
	}
	std::cout << "Done printing all zeros\n";

}

void QuadraticFunction::PrintAllRealNumberZeros() const
{
	std::cout << "Printing all real number zeros of the function\n";

	for (const auto& zero : m_RealNumberZeros)
	{
		std::cout << zero << std::endl;
	}
	std::cout << "Done printing all real number zeros\n";


}

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
	PrintAllRealNumberZeros();
	std::cout << std::endl;

	std::cout << "End of data output\n";
}


std::vector<double> GetZerosQuadraticFormula(QuadraticFunction& QuadraticFunc)
{
	std::tuple<double, double, double> ABC = QuadraticFunc.GetABC();

	double a = std::get<0>(ABC);
	double b = std::get<1>(ABC);
	double c = std::get<2>(ABC);


	std::vector<double> LocalVecOfZeros = GetZerosQuadraticFormula(a, b, c);
	std::vector<double> OutVecRealNumZeros;

	//std::cout << LocalVecOfZeros.size();

	for (int i = 0; i < LocalVecOfZeros.size(); ++i)
	{
		double TermOne, TermTwo, TermThree;

		TermOne = (std::pow(LocalVecOfZeros[i], 2)) * a;
		TermTwo = LocalVecOfZeros[i] * b;
		TermThree = c;

		double ZeroTest = TermOne + TermTwo + TermThree;

		//std::cout << ZeroTest << std::endl;

		if (is_close_to_zero(ZeroTest))
		{
			OutVecRealNumZeros.push_back(LocalVecOfZeros[i]);
		}
	}


	QuadraticFunc.SetAllZeroVec(LocalVecOfZeros);
	QuadraticFunc.SetRealNumberZeroVec(OutVecRealNumZeros);

	// Line of symmetry is between the two zeros.
	QuadraticFunc.SetLineOfSymmetry(LocalVecOfZeros[1] / LocalVecOfZeros[0]);
	QuadraticFunc.SetTheMaxMinValue(LocalVecOfZeros[1] / LocalVecOfZeros[0]);


	// dont really need to return it but whatever...
	return OutVecRealNumZeros;

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

void PrintFunctionTransformationInfo()
{
	// Transformation of f(c>0)	Effect on the graph of f
	/*	f(x)+c	Vertical shift up c units
	f(x)−c	Vertical shift down c units
	f(x+c)	Shift left by c units
	f(x−c)	Shift right by c units
	c*f(x)
	Vertical stretch if c>1;
	vertical compression if 0<c<1

	f(cx)	Horizontal stretch if 0<c<1; horizontal compression if c>1
	−f(x)Reflection about the x-axis
	f(−x)	Reflection about the y-axis*/
}

double TrigonometricFunction::FixAngleBetweenZeroAndTwoRad(const double & InAngle)
{
	double FixedAngle = InAngle;

	while (FixedAngle > (2 * M_PI))
	{
		FixedAngle = FixedAngle - (2 * M_PI);
	}
	while (FixedAngle < 0)
	{
		FixedAngle = FixedAngle + (2 * M_PI);
	}

	return FixedAngle;
}

double TrigonometricFunction::SinOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
{
	double FixedAngle = InAngle;

	while (FixedAngle > (2 * M_PI))
	{
		FixedAngle = FixedAngle - (2 * M_PI);
	}
	while (FixedAngle < 0)
	{
		FixedAngle = FixedAngle + (2 * M_PI);
	}

	auto Map = InUnitCircle.GetRadianSinCosMap();

	double SinOfAngle = Map[FixedAngle].second;

	return SinOfAngle;
}

double TrigonometricFunction::CosOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
{
	double FixedAngle = InAngle;

	while (FixedAngle > (2 * M_PI))
	{
		FixedAngle = FixedAngle - (2 * M_PI);
	}
	while (FixedAngle < 0)
	{
		FixedAngle = FixedAngle + (2 * M_PI);
	}

	auto Map = InUnitCircle.GetRadianSinCosMap();

	double CosOfAngle = Map[FixedAngle].first;

	return CosOfAngle;
}

double TrigonometricFunction::TanOfUnitCircleAngle(const UnitCircle & InUnitCircle, const double & InAngle)
{
	double FixedAngle = InAngle;

	while (FixedAngle > (2.0 * M_PI))
	{
		FixedAngle = FixedAngle - (2.0 * M_PI);
	}
	while (FixedAngle < 0)
	{
		FixedAngle = FixedAngle + (2.0 * M_PI);
	}

	auto Map = InUnitCircle.GetRadianSinCosMap();


	double CosOfAngle = Map[FixedAngle].first;
	double SinOfAngle = Map[FixedAngle].second;


	// TODO: used this as debug code (remove later once you memorize numeric limits stuff)
	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
	//	<< SevenM_PIOverFour << std::endl;

	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
	//	<< FixedAngle << std::endl;

	return SinOfAngle / CosOfAngle;
}


double FindCompoundIntrest(const double & P, const double & r, const double & t)
{
	const double Exponent = r*t;
	const double Base = Eulere;

	double TempResult = std::pow(Base, Exponent);
	TempResult = TempResult * P;

	return TempResult;
}

double FindCompoundIntrest(
	const double& P,
	const double& r,
	const double& n,
	const double& t)
{
	const double Exponent = n*t;
	const double InsideParanthesis = 1.0 + (r / n);

	double TempResult = std::pow(InsideParanthesis, Exponent);
	TempResult = TempResult * P;

	return TempResult;

}

double Limit::SimpifyComplexFraction(const ComplexFraction & InComplexFract)
{
	//ComplexFraction OutComplexFract;

	// Run the test so I can see the limitBothZeroCheck result used just below
	double LimitTest = InComplexFract(m_a);

	if ((InComplexFract.GetLimitBothZeroCheck()) == true && m_bTryConjSolution == false)
	{

		// First simplify through multipliction
		// Get common denominator
		LinearFunction FirstNeededTerm = InComplexFract.GetNumAData().second;
		int SecondNeededTerm = InComplexFract.GetNumConstFractData().second;



		LinearFunction NewFactoredNumerator;
		LinearFunction NewFactoredDenominator;

		

		std::vector<double> NumeratorZeros;
		NumeratorZeros.push_back(SecondNeededTerm);
		NumeratorZeros.push_back(FirstNeededTerm.GetB());

		
		std::vector<double> DenominatorZeros;
		DenominatorZeros.push_back(InComplexFract.GetNumAData().second.GetB());

		std::vector<double> NewNumerator;

		for (int i = 0; i < NumeratorZeros.size(); ++i)
		{
			auto result1 = std::find(std::begin(DenominatorZeros), std::end(DenominatorZeros), NumeratorZeros[i]);

			if (result1 != std::end(DenominatorZeros))
			{
				std::cout << "Numerator contains: " << NumeratorZeros[i] << '\n';
			}
			else
			{
				std::cout << "Numerator does not contain: " << NumeratorZeros[i] << '\n';
				NewNumerator.push_back(NumeratorZeros[i]);
			}
		}

		for (int i = 0; i < NewNumerator.size(); ++i)
		{
			std::cout << NewNumerator[i] << " ";
			std::cout << std::endl;
		}

		//double One = std::abs(std::floor(NewNumerator[0]));
		//double Two = NewNumerator[0];

		//std::cout << "abs " << One << std::endl;
		//std::cout << Two << std::endl;

		if (/*std::abs*/(std::floor(NewNumerator[0])) == NewNumerator[0])
		{
			std::cout << "NewNumerator[0] is whole\n";


			if (NewNumerator[0] == 0.0)
			{
				NewFactoredNumerator = LinearFunction(1, 0);
			}
			else
			{
				if (NewNumerator[0] < 0.0)
				{
					NewFactoredNumerator = LinearFunction(1, NewNumerator[0] * (-1));
				}
				else
				{
					NewFactoredNumerator = LinearFunction(1, NewNumerator[0]);
				}
			}
		}
		else
		{
			std::cout << "NewNumerator[0] is not whole\n";

			std::pair<double, double> NumeratorFract = OutputDecimalAsFract(NewNumerator[0]);

			std::cout << NumeratorFract.first << "/" << NumeratorFract.second << std::endl;

			if (NumeratorFract.first < 0)
			{
				NewFactoredNumerator = LinearFunction(NumeratorFract.second, NumeratorFract.first);
				//std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
			}
		}



		std::vector<double> NewDenominator;

		for (int i = 0; i < NumeratorZeros.size(); ++i)
		{
			auto result1 = std::find(std::begin(NumeratorZeros), std::end(NumeratorZeros), DenominatorZeros[i]);

			if (result1 != std::end(NumeratorZeros))
			{
				std::cout << "Numerator contains: " << DenominatorZeros[i] << '\n';
			}
			else
			{
				std::cout << "Numerator does not contain: " << DenominatorZeros[i] << '\n';
				NewDenominator.push_back(DenominatorZeros[i]);
			}
		}

		for (int i = 0; i < NewDenominator.size(); ++i)
		{
			std::cout << NewDenominator[i] << " ";
			std::cout << std::endl;
		}

		if (std::abs(std::floor(NewDenominator[0])) == NewDenominator[0])
		{
			std::cout << "NewDenominator[0] is whole\n";

			if (NewDenominator[0] == 0.0)
			{
				NewFactoredDenominator = LinearFunction(1, 0);
			}
			else
			{
				if (NewDenominator[0] < 0.0)
				{
					NewFactoredDenominator = LinearFunction(1, NewDenominator[0]);
				}
				else
				{
					NewFactoredDenominator = LinearFunction(1, NewDenominator[0] * (-1));
				}
			}

		}
		else
		{
			std::cout << "NewDenominator[0] is not whole\n";

			std::pair<double, double> DenominatorFract = OutputDecimalAsFract(NewDenominator[0]);
			if (DenominatorFract.first < 0)
			{
				std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
				NewFactoredDenominator = LinearFunction(DenominatorFract.second, DenominatorFract.first*(-1));

			}
		}

		NewFactoredNumerator.PrintLinearFunctionInfo();
		NewFactoredDenominator.PrintLinearFunctionInfo();

		double FinalFactoredNumerator = EvaluateLinearFuncLimit(NewFactoredNumerator);
		double FinalFactoredDenominator = EvaluateLinearFuncLimit(NewFactoredDenominator);

		std::cout << "FF Num: " << FinalFactoredNumerator << std::endl;
		std::cout << "FF Den: " << FinalFactoredDenominator << std::endl;

		return FinalFactoredNumerator / FinalFactoredDenominator;

	}



	//return OutComplexFract;
}

double Limit::EvaluateRootFuncLimit(const RootFunction & InRootFunc)
{
	bool bDoesNegDirLimitExist = true;


	double Limit = InRootFunc(m_a);

	double DomainStart = InRootFunc.GetStartingDomainNum();

	

	// TODO: Find a way to use the starting domain number to trigger a bool that will display 
	// Limit Does not Exist for the side that DNE
	// UPDATE: sort of did it needs more polishing though

	std::vector<std::pair<double, double>> PosDirVec;
	// pos direction
	std::pair<double, double> PosPairFirst;
	PosPairFirst.first = m_a + 0.1;
	PosPairFirst.second = InRootFunc(m_a + 0.1);

	std::pair<double, double> PosPairSecond;
	PosPairSecond.first = m_a + 0.01;
	PosPairSecond.second = InRootFunc(m_a + 0.01);

	std::pair<double, double> PosPairThird;
	PosPairThird.first = m_a + 0.001;
	PosPairThird.second = InRootFunc(m_a + 0.001);

	std::pair<double, double> PosPairFourth;
	PosPairFourth.first = m_a + 0.0001;
	PosPairFourth.second = InRootFunc(m_a + 0.0001);

	PosDirVec.push_back(PosPairFirst);
	PosDirVec.push_back(PosPairSecond);
	PosDirVec.push_back(PosPairThird);
	PosDirVec.push_back(PosPairFourth);

	std::vector<std::pair<double, double>> NegDirVec;
	// neg direction
	std::pair<double, double> NegPairFirst;
	NegPairFirst.first = m_a - 0.1;
	NegPairFirst.second = InRootFunc(m_a - 0.1);

	if (std::isnan(NegPairFirst.second))
	{
		bDoesNegDirLimitExist = false;
	}

	std::pair<double, double> NegPairSecond;
	NegPairSecond.first = m_a - 0.01;
	NegPairSecond.second = InRootFunc(m_a - 0.01);

	std::pair<double, double> NegPairThird;
	NegPairThird.first = m_a - 0.001;
	NegPairThird.second = InRootFunc(m_a - 0.001);

	std::pair<double, double> NegPairFourth;
	NegPairFourth.first = m_a - 0.0001;
	NegPairFourth.second = InRootFunc(m_a - 0.0001);

	NegDirVec.push_back(NegPairFirst);
	NegDirVec.push_back(NegPairSecond);
	NegDirVec.push_back(NegPairThird);
	NegDirVec.push_back(NegPairFourth);

	std::cout << "Evaluating Limit: Please Wait...\n";

	for (auto & num : PosDirVec)
	{

		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< std::setw(7) << num.first << " " << num.second << std::endl;
	}

	std::cout << std::endl;

	for (auto & num : NegDirVec)
	{
		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< std::setw(7) << num.first << " " << num.second << std::endl;
	}

	double TopRes = PosDirVec[3].second;
	double BottomRes = NegDirVec[3].second;
	/*std::cout << TopRes << std::endl;
	std::cout << BottomRes << std::endl;*/


	std::string TopResStr = std::to_string(TopRes);
	// TODO: remove debug code later
	std::cout << TopResStr << std::endl;

	std::string BottomResStr = std::to_string(BottomRes);
	// TODO: remove debug code later
	std::cout << BottomResStr << std::endl;

	double LocalPosRes = std::stod(TopResStr) * 100;

	double LocalNegRes = std::stod(BottomResStr) * 100;

	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
	//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;



	std::cout << std::endl;

	// Are these limits infinite?
	bool bIsPosDirPosInfinity = false;
	bool bIsPosDirNegInfinity = false;
	bool bIsNegDirNegInfinity = false;
	bool bIsNegDirPosInfinity = false;

	

	if (PosDirVec[1].second > PosDirVec[0].second)
	{
		if (PosDirVec[2].second > PosDirVec[1].second)
		{
			bIsPosDirPosInfinity = true;
		}
	}

	if (NegDirVec[1].second < NegDirVec[0].second)
	{
		if (NegDirVec[2].second < NegDirVec[1].second)
		{
			bIsNegDirNegInfinity = true;
		}
	}

	if (PosDirVec[1].second < PosDirVec[0].second)
	{
		if (PosDirVec[2].second < PosDirVec[1].second)
		{
			bIsPosDirNegInfinity = true;
		}
	}
	if (NegDirVec[1].second > NegDirVec[0].second)
	{
		if (NegDirVec[2].second > NegDirVec[1].second)
		{
			bIsNegDirPosInfinity = true;
		}
	}

	if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
	{
		std::cout << "We have a working limit result! 1\n";


		if (bIsPosDirPosInfinity)
		{
			if (bIsNegDirPosInfinity)
			{
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << INFINITY << std::endl;
				return INFINITY;
			}
		}

		if (bIsPosDirNegInfinity)
		{
			if (bIsNegDirNegInfinity)
			{
				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
				std::cout << NEGINFINITY << std::endl;
				return NEGINFINITY;
			}
		}

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

		//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


		//	// return either
		//return std::ceil(LocalPosRes) / 100;
	}

	// if you have two results that are returning negatives
	// you need to do some sort of swap with floor / ceil
	if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
	{
		std::cout << "We have a working limit result! 2\n";

		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

		std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


		// return either
		return std::floor(LocalPosRes) / 100;
	}


	std::cout << "Limit: DNE (does not exist): \n\n";

	std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";

	if (bIsPosDirPosInfinity)
	{
		std::cout << INFINITY << std::endl;
	}
	else
	{
		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
	}

	if (bDoesNegDirLimitExist == false)
	{
		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = DOES NOT EXIST";
	}
	else
	{
		// The limit Does Exist
		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";

		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
		}
	}
	
	

	std::cout << std::endl;

	return 0.0;
}



bool DetermineContinunityAtAPoint(const RationalFunction& InFunc, const int& InPoint)
{
	QuadraticFunction Numerator = InFunc.GetNumeratorQuadratic();
	LinearFunction Denominator = InFunc.GetDenominatorLinear();

	// Step 1: check to see if f(a) is defined
	double TestOne = Numerator(InPoint) / Denominator(InPoint);

	if (std::isnan(TestOne))
	{
		// failed 
		std::cout << "TestOneFailed: The function is not continuous at " << InPoint << "\n";
		return false;
	}
	else
	{
		//  If f(a) is defined, continue to step 2.

		// Step 2: Compute Limit from both sides
		// If Limit does not exist (that is, it is not a real number),
		// then the function is not continuous at a and the problem is solved.



		Limit TestTwoLimit(InFunc, InPoint);
		double TestTwo = TestTwoLimit.GetLimitResult();

		// if Limit Exists go to step 3
		// TODO: Need a rational function check to return bool to see if it exists
		// TODO: Need more example data in order to fix

		if (TestOne != TestTwo)
		{
			std::cout << "TestThreeFailed: The function is not continuous at " << InPoint << "\n";
			return false;
		}
		else
		{
			std::cout << "TestThreePassed: The function is continuous at " << InPoint << "\n";
			return true;
		}

	}
}
