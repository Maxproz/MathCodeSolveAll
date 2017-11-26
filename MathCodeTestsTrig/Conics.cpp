#include "Conics.h"

#include <exception>
#include <iostream>
#include <complex>
#include <cmath>
#include <cctype>
#include <locale>         // std::locale, std::tolower

const double MATHPI = 3.14;

ConicType IdentifyConicType(
	const double& A,
	const double& C,
	const double& D,
	const double& E,
	const double& F,
	const double& B)
{
	if (A == 0 && B == 0 && C == 0)
	{
		throw std::exception("A B C cant all == 0");
	}

	// general form
	// Ax^2+Bxy+Cy^2+Dx+Ey+F=0.

	bool AAndCBothNotZero = (A != 0 && C != 0);
	bool AAndCHaveSameSign = ((A < 0 && C < 0) || (A > 0 && C > 0));
	bool AAndCNotEqual = (A != C);
	if (AAndCBothNotZero && AAndCHaveSameSign && AAndCNotEqual)
	{
		// may be an ellipse
		return ConicType::Ellipse;
	}

	bool AAndCEqual = (A == C);

	if (AAndCEqual && AAndCBothNotZero && AAndCHaveSameSign)
	{
		// may be circle
		return ConicType::Circle;
	}

	if (AAndCBothNotZero && (!AAndCHaveSameSign))
	{
		// may be hyperbola 
		return ConicType::Hyperbola;
	}

	if (A == 0 || C == 0)
	{
		// may be parabola 
		return ConicType::Parabola;
	}

}

void PrintConicType(const ConicType& InType)
{
	switch (InType)
	{
		case ConicType::Circle:
		{
			std::cout << "Circle" << std::endl;
			break;
		}
		case ConicType::Ellipse:
		{
			std::cout << "Ellipse" << std::endl;
			break;
		}
		case ConicType::Hyperbola:
		{
			std::cout << "Hyperbola" << std::endl;
			break;
		}
		case ConicType::Parabola:
		{
			std::cout << "Parabola" << std::endl;
			break;
		}
		default:
		{
			throw std::exception("Something went wrong in PrintConicTypeFunction");
		}
	}
}

// Unfinished function
ConicType IdentifyConicType(const double & A, const double & C, const double & D, const double & E, const double & F, const double & RotationAngleRadians, const double & B)
{

	double X = std::cos(RotationAngleRadians) - std::sin(RotationAngleRadians);
	double Y = std::sin(RotationAngleRadians) + std::cos(RotationAngleRadians);


	double LocalA = A*X;
	double LocalD = D*X;

	double LocalC = C*Y;
	double LocalE = E*Y;



	// TODO: From here you need to FOIL
	// COMBINE LIKE TERMS
	// COMBINE LIKE TERMS
	// Multiply both sides by 2
	// Simplify 
	// Distribute
	// Set Equal to 1 (divide left side by the right side)
	// Then you should get x prime and y prime in standard form.
	
	// Not sure how to do this in c++ with the variables etc...
	// Would need to be done by hand or with an online calculator at
	// least until my skills improve.


	if (A == 0 && B == 0 && C == 0)
	{
		throw std::exception("A B C cant all == 0");
	}

	// general form
	// Ax^2+Bxy+Cy^2+Dx+Ey+F=0.

	bool AAndCBothNotZero = (A != 0 && C != 0);
	bool AAndCHaveSameSign = ((A < 0 && C < 0) || (A > 0 && C > 0));
	bool AAndCNotEqual = (A != C);
	if (AAndCBothNotZero && AAndCHaveSameSign && AAndCNotEqual)
	{
		// may be an ellipse
		return ConicType::Ellipse;
	}

	bool AAndCEqual = (A == C);

	if (AAndCEqual && AAndCBothNotZero && AAndCHaveSameSign)
	{
		// may be circle
		return ConicType::Circle;
	}

	if (AAndCBothNotZero && (!AAndCHaveSameSign))
	{
		// may be hyperbola 
		return ConicType::Hyperbola;
	}

	if (A == 0 || C == 0)
	{
		// may be parabola 
		return ConicType::Parabola;
	}
}

void WriteEquationOfRatatedConicStandardForm(
	const double& A,
	const double& B,
	const double& C,
	const double& Equals)
{
	// Basic Notes:
	// If  cot(2θ)>0, then 2θ  is in the first quadrant, and θ  θ  is between (0°,45°). (0°,45°).
	// If, cot(2θ)<0, then 2θ  is in the second quadrant, and θ  θ  is between (45°, 90°). (45°, 90°).
	// If A = C, A = C, then θ = 45°.


	// 1. Find cot(2θ).cot(2θ).
	double COT2Theta = (A - C) / B;

	// Numerator == adjacent
	// Denominator == opposite
	std::pair<double, double> COT2ThetaPair = OutputDecimalAsFract(COT2Theta);

	double Hypotenuse = std::sqrt(
		std::pow(COT2ThetaPair.first,2) + 
		std::pow(COT2ThetaPair.second, 2));

	// TODO: Remove debug code
	std::cout << Hypotenuse << std::endl;

	// 2. Find sin θ and cos θ
	// TODO: hmmm maybe change how i do this later
	//std::pair<double, double> COS2ThetaPair(COT2ThetaPair.first, Hypotenuse);
	double COS2Theta = COT2ThetaPair.first / Hypotenuse;

	double SinTheta = std::sqrt((1 - COS2Theta) / 2.0);
	double CosTheta = std::sqrt((1 + COS2Theta) / 2.0);

	// 3. Substitute  sin θ and cos θ into x = x′cos θ−y′sin θ and θ. y = x′sin θ + y′cos θ.
	// TODO: How the fuck am i going to show this
	//double X = 

	// 4. Substitute the expression for x and y into in the given equation, and then simplify.
	
	// 4. Notes
	// Can i use std::complex to trap the numerator and denominators here so I 
	// can return useful information somehow?

	// 5. Write the equations with  x′ and  y′ in the standard form with respect to the rotated axes.


}


ConicType UseDescrimentToFindConicType(const double & B, const double & A, const double & C)
{
	ConicType OutType;

	double Discriminant = (std::pow(B, 2) - (4 * A*C));

	if (Discriminant < 0)
	{
		OutType = ConicType::Ellipse;
	}
	else if (Discriminant == 0)
	{
		OutType = ConicType::Parabola;
	}
	else if (Discriminant > 0)
	{
		OutType = ConicType::Hyperbola;
	}
	else
	{
		throw std::exception("Something went wrong in Discriminent Evaluation func");
	}

	return OutType;

}

ConicType GetTypeFromEccentricity(const double & Eccentricity)
{
	ConicType OutType;

	if (Eccentricity == 0 || Eccentricity < 1)
	{
		OutType = ConicType::Ellipse;
	}
	else if (Eccentricity == 1)
	{
		OutType = ConicType::Parabola;
	}
	else if (Eccentricity > 1)
	{
		OutType = ConicType::Hyperbola;
	}
	else
	{
		throw std::exception("Something went wrong in GetTypeFromEccentricity");
	}

	return OutType;
}

// Take in a function provided the basic function info 
// Return the rewritten function
PolarEquation EvaluatePolarEquation(PolarEquation& InEquation)
{
	// 1. Multiply the numerator and denominator by the reciprocal
	// of the constant in the denominator.

	double Recipriocal = (1.0 / InEquation.m_DenominatorConstant);

	double LocalDenomConstant = InEquation.m_DenominatorConstant * Recipriocal;
	double LocalNumerator = InEquation.m_Numerator * Recipriocal;
	double LocalDenomMult = InEquation.m_DenominatorFuncMultiple * Recipriocal;

	double P{ 0.0 };
	double Eccentricity = LocalDenomMult;

	InEquation.m_Eccentricity = Eccentricity;

	P = (LocalNumerator * (1.0 / Eccentricity));
	
	std::string DirectrixString;
	ConicType FuncType = GetTypeFromEccentricity(Eccentricity);

	if (InEquation.m_Function == TrigFunctions::SIN)
	{
		// directrix: y = p
		if (InEquation.m_PoM == PlusOrMinusDenom::Plus)
		{
			DirectrixString = "y = ";
			
		}
		else
		{
			// Minus
			DirectrixString = "y = -";
		}

	}
	else if (InEquation.m_Function == TrigFunctions::COS)
	{
		// directrix: x = p
	
		if (InEquation.m_PoM == PlusOrMinusDenom::Plus)
		{
			DirectrixString = "x = ";
		}
		else
		{
			// Minus
			DirectrixString = "x = -";
		}

	}
	else
	{
		throw std::exception("Something went wrong in polar equation evaluation");
	}

	DirectrixString.append(std::to_string(P));
	while (DirectrixString.size() >= 9)
	{
		DirectrixString.pop_back();
	}

	InEquation.m_DirectrixString = DirectrixString;

	PrintPolarEquationInfo(DirectrixString, Eccentricity, FuncType);

	
	return PolarEquation(LocalNumerator, LocalDenomConstant, Eccentricity,
		InEquation.m_Function, InEquation.m_PoM);

}

void PrintPolarEquationInfo(
	const std::string& InDirectrixStr,
	const double& InEccentricity,
	const ConicType& InType)
{
	std::cout << "conic with focus at the origin" << std::endl;

	std::cout << "The conic type is a: ";
	PrintConicType(InType);

	std::cout << "eccentricity = " << InEccentricity << std::endl;
	//PrintDecimalAsFraction(InEccentricity);

	std::cout << "Directrix: " << InDirectrixStr << std::endl;
}

std::map<std::string, double> GetPolarGraphInfo(PolarEquation& InFunc)
{
	std::map<std::string, double> OutGraphMap;

	//ConicType FuncType = GetTypeFromEccentricity(EvaluatedFunc.m_Eccentricity);
	PolarEquation EvaluatedFunc = EvaluatePolarEquation(InFunc);

	std::string A = "0";
	std::string B = "pi/2";
	std::string C = "pi";
	std::string D = "3pi/2";

	double FirstInput = RunPolarFunction(EvaluatedFunc, 0.0);
	double SecondInput = RunPolarFunction(EvaluatedFunc, MATHPI / 2.0);
	double ThirdInput = RunPolarFunction(EvaluatedFunc, MATHPI);
	double FourthInput = RunPolarFunction(EvaluatedFunc, (3.0 * MATHPI) / 2.0);

	OutGraphMap[A] = FirstInput;
	OutGraphMap[B] = SecondInput;
	OutGraphMap[C] = ThirdInput;
	OutGraphMap[D] = FourthInput;

	return OutGraphMap;
}


// TODO: figure out a way to know if a function output is undefined instead
// of adding it to the data.
double RunPolarFunction(const PolarEquation& InFunction, const double& AngleInput)
{
	double Output{ 0.0 };

	if (InFunction.m_Function == TrigFunctions::SIN)
	{

		if (InFunction.m_PoM == PlusOrMinusDenom::Plus)
		{
			Output = InFunction.m_Numerator /
				((std::sin(AngleInput) * InFunction.m_DenominatorFuncMultiple) +
					InFunction.m_DenominatorConstant);

			return Output;
		}
		else
		{
			// Minus 
			Output = InFunction.m_Numerator /
				((std::sin(AngleInput) * InFunction.m_DenominatorFuncMultiple) -
					InFunction.m_DenominatorConstant);

			return Output;
		}

	}
	else if (InFunction.m_Function == TrigFunctions::COS)
	{
		if (InFunction.m_PoM == PlusOrMinusDenom::Plus)
		{
			Output = (InFunction.m_Numerator / 
				((std::cos(AngleInput) * InFunction.m_DenominatorFuncMultiple) +
					InFunction.m_DenominatorConstant));
		
			return Output;
		}
		else
		{
			// Minus 
			Output = InFunction.m_Numerator / 
				((std::cos(AngleInput) * InFunction.m_DenominatorFuncMultiple) -
					InFunction.m_DenominatorConstant);
			
			return Output;
		}
	}
	else
	{
		throw std::exception("Failed to run a polar function");
	}

	
}

// When eccentricity is a whole number
PolarEquation FindPolarFormGivenEccentricityAndDirectrix(
	const double& Eccentricity,
	const std::string& Directrix)
{
	std::string LocalDirectrix = Directrix;

	std::stringstream StringStream;
	StringStream << LocalDirectrix;

	std::string XorY;
	std::string EqualSign;
	double DirectrixValue;
	

	StringStream >> XorY >> EqualSign >> DirectrixValue;
	
	bool bIsPostiveDirectrix = DirectrixValue > 0.0;

	std::locale loc;
	for (auto elem : XorY)
	{
		std::tolower(elem, loc);
	}
	
	double Outp = std::abs(DirectrixValue);
	double Oute = Eccentricity;
	double OutNumerator = Outp*Oute;
	double OutDenomConst{ 1 };
	PlusOrMinusDenom OutPoM;
	TrigFunctions OutType;

	if (XorY == "x")
	{
		OutType = TrigFunctions::COS;

		// we know that we use cosine
		if (bIsPostiveDirectrix)
		{
			OutPoM = PlusOrMinusDenom::Plus;
		}
		else
		{
			OutPoM = PlusOrMinusDenom::Minus;
		}

	}
	else if (XorY == "y")
	{
		OutType = TrigFunctions::SIN;

		// we know that we use sine
		if (bIsPostiveDirectrix)
		{
			OutPoM = PlusOrMinusDenom::Plus;
		}
		else
		{
			OutPoM = PlusOrMinusDenom::Minus;
		}
	}
	else
	{
		throw std::exception("Error inside of FindPolarformGiveneccentri...etc");
	}

	//ConicType OutConicType = GetTypeFromEccentricity(Eccentricity);

	return PolarEquation(OutNumerator, OutDenomConst, Oute, OutType, OutPoM);
}