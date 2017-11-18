// MathCodeTestsTrig.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include <iostream>
#include <string>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <functional>
#include <map>

#include "TrigFunctions.h"
#include "LawOfCosines.h"
#include "ComplexNumber.h"
#include "MiscMathEquations.h"
#include "ParametricEquation.h"
#include "MPVector.h"
#include "MPMatrix.h"
#include "Ellipse.h"
#include "Hyperbola.h"
#include "Parabola.h"
#include "Conics.h"
#include "Sequence.h"
#include "CountingPrinciples.h"
#include "BionomialTheorm.h"
#include "Probability.h"
#include "CalculusFunction.h"
#include "Circle.h"

typedef std::pair<double, double> Point;
//#define NEGINFINITY   ((float)(_HUGE_ENUF * _HUGE_ENUF) *(-1))


const double my_gravityfeet(32);  // f/s^2


double GetSINEOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle);
double GetSINEOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle);

long mygcd(long a, long b);

enum class GivenSides
{
	AB,
	AC,
	BC
};

enum class GivenAngle
{
	Alpha,
	Beta,
	Gamma
};

enum class DoubleAngleReturnType
{
	SIN,
	COS,
	TAN
};

enum class HalfAngleFormula
{
	Sin,
	Cos,
	Tan
};

enum class GivenSide
{
	A,
	B,
	C
};




std::pair<double, double> OutputDecimalAsFraction(double input)
{
	double integral = std::floor(input);
	double frac = input - integral;

	const long precision = 100; // This is the accuracy. // was 1 million or so

	long gcd_ = mygcd(round(frac * precision), precision);

	long denominator = precision / gcd_;
	long numerator = round(frac * precision) / gcd_;
	
	/*
	if (integral != 0)
	{
		std::cout << integral << " + ";
	}*/


	//std::cout << numerator << "/" << denominator << std::endl;


	return std::pair<double, double>(numerator, denominator);
}

void PrintDecimalAsFraction(double input)
{
	double integral = std::floor(input);
	double frac = input - integral;

	const long precision = 100; // This is the accuracy. // was 1 million or so

	long gcd_ = mygcd(round(frac * precision), precision);

	long denominator = precision / gcd_;
	long numerator = round(frac * precision) / gcd_;

	/*
	if (integral != 0)
	{
	std::cout << integral << " + ";
	}*/

	std::cout << numerator << "/" << denominator << std::endl;
}

long mygcd(long a, long b)
{
	if (a == 0)
		return b;
	else if (b == 0)
		return a;

	if (a < b)
		return mygcd(a, b % a);
	else
		return mygcd(b, a % b);
}

// https://cnx.org/contents/E6wQevFf@6.38:aIPS8_HQ@8/Inverse-Trigonometric-Function
// Evaluating the Composition of a Sine with an Inverse Cosine
// We can use the Pythagorean identity to do this.
double CompositionofaSineWithAnInverseCosine(double InCos)
{
	// set up formula
	// sin^2(theta) + cos^2(theta) = 1
	// Use our known value for cosine
	InCos = std::pow(InCos, 2);
	// Solve for sine by subtracting our cosine from the other side(1)
	double RightHandSide = 1 - InCos;
	// Now need to take the square root to get sine off of a power.
	RightHandSide = std::sqrt(RightHandSide);

	// We know that the inverse cosine always gives an angle on the interval [0, π], 
	// so we know that the sine of that angle must be positive; 
	// therefore sin(cos^−1(4/5)) = sinθ = 3/5.
	return RightHandSide;
}


// TODO: This function needs to be set up to detect which quadrant the result is in
// Evaluating the Composition of a Sine with an Inverse Sine
// We can use the Pythagorean identity to do this.
double CompositionofaCosineWithAnInverseSine(double InSine)
{
	// set up formula
	// sin^2(theta) + cos^2(theta) = 1
	// Use our known value for Sine
	InSine = std::pow(InSine, 2);
	// Solve for sine by subtracting our cosine from the other side(1)
	double RightHandSide = 1 - InSine;
	// Now need to take the square root to get sine off of a power.
	RightHandSide = std::sqrt(RightHandSide);

	// We know that the inverse cosine always gives an angle on the interval [-π/2, π/2], 
	// so we know that the cosine of that angle must be positive; 
	// therefore cos(sin^−1(7/9)) = cosθ = 4*sqrt(2)/9
	return RightHandSide;
}


// Reciprocal Identities

double GetReciprocalIdentity(const double Angle)
{
	// sinθ = 1/cscθ
	// cosθ = 1/secθ
	// tanθ = 1/cotθ
	// csc θ = 1/sin θ
	// sec θ = 1/cos θ
	// cot  θ = 1/tan θ
	return 0.0;
}

// SUM AND DIFFERENCE FORMULAS FOR COSINE

// Cosine of the Difference of Two Angles
// Assumes the user entered the angles in radians
double GetCosineOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	// cos(alpha - beta) = cos(alpha)*cos(beta) + sin(alpha)*sin(beta)
	double OutResult = std::cos(AlphaAngle)*std::cos(BetaAngle) + 
						std::sin(AlphaAngle)*std::sin(BetaAngle);

	return OutResult;

}

// Cosine of the Sum of Two Angles
// Assumes the user entered the angles in radians
double GetCosineOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	// cos(alpha - beta) = cos(alpha)*cos(beta) + sin(alpha)*sin(beta)
	double OutResult = std::cos(AlphaAngle)*std::cos(BetaAngle) -
		std::sin(AlphaAngle)*std::sin(BetaAngle);

	return OutResult;

}


double GetProductOfSumCosACosB(const double& AlphaAngle, const double& BetaAngle)
{
	double OutResult{ 0.0 };

	OutResult = ((GetCosineOfTheDifferenceOfTwoAngles(AlphaAngle, BetaAngle) +
		GetCosineOfTheSumOfTwoAngles(AlphaAngle, BetaAngle)) / 2);

	return OutResult;

}

double GetProductOfSumSinACosB(const double& AlphaAngle, const double& BetaAngle)
{
	double OutResult{ 0.0 };

	OutResult = ((GetSINEOfTheSumOfTwoAngles(AlphaAngle, BetaAngle) +
		GetSINEOfTheDifferenceOfTwoAngles(AlphaAngle, BetaAngle)) / 2);

	return OutResult;
}

double GetProductOfSumSinASinB(const double& AlphaAngle, const double& BetaAngle)
{
	double OutResult{ 0.0 };

	OutResult = ((GetCosineOfTheDifferenceOfTwoAngles(AlphaAngle, BetaAngle) -
		GetCosineOfTheSumOfTwoAngles(AlphaAngle, BetaAngle)) / 2);

	return OutResult;
}

double GetProductOfSumCosASinB(const double& AlphaAngle, const double& BetaAngle)
{
	double OutResult{ 0.0 };

	OutResult = ((GetSINEOfTheSumOfTwoAngles(AlphaAngle, BetaAngle) -
		GetSINEOfTheDifferenceOfTwoAngles(AlphaAngle, BetaAngle)) / 2);

	return OutResult;
}

double SumToProductSinAddSin(const double& AlphaAngle, const double& BetaAngle)
{
	return 0.0;
}

double SumToProductCosSubtractCos(const double& AlphaAngle, const double& BetaAngle)
{
	// assumes the input angles are in radians

	double OutResult =
		(std::sin((AlphaAngle + BetaAngle) / 2.0) *
			std::sin((AlphaAngle - BetaAngle) / 2.0) * (-2.0));

	return OutResult;

}

double SumToProductCosSubtractCosDegrees(const double& AlphaAngle, const double& BetaAngle, const bool& ReturnRadians = true)
{
	double AlphaAngRad = DegreesToRadians(AlphaAngle);
	double BetaAngRad = DegreesToRadians(BetaAngle);

	double OutResult = 
		(std::sin((AlphaAngRad + BetaAngRad) / 2.0) *
			std::sin((AlphaAngRad - BetaAngRad) / 2.0) * (-2.0));

	if (ReturnRadians == false)
	{
		// return degrees change variable
		OutResult = RadianToDegrees(OutResult);

	}
	else
	{
		// nothing
	}

	return OutResult;

}


// SUM AND DIFFERENCE FORMULAS FOR SINE

// SINE of the Difference of Two Angles
// Assumes the user entered the angles in radians
double GetSINEOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	// cos(alpha - beta) = cos(alpha)*cos(beta) + sin(alpha)*sin(beta)
	double OutResult = std::sin(AlphaAngle)*std::cos(BetaAngle) -
						std::cos(AlphaAngle)*std::sin(BetaAngle);

	return OutResult;

}

// SINE of the Sum of Two Angles
// Assumes the user entered the angles in radians
double GetSINEOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	// cos(alpha - beta) = cos(alpha)*cos(beta) + sin(alpha)*sin(beta)
	double OutResult = std::sin(AlphaAngle)*std::cos(BetaAngle) +
						std::cos(AlphaAngle)*std::sin(BetaAngle);
		

	return OutResult;
}

// https://cnx.org/contents/E6wQevFf@6.38:LKVQKWmT@5/Sum-and-Difference-Identities
// The sum and difference formulas for tangent are:
// tan(alpha + beta) = tan(alpha) + tan(beta) / 1 - tan(alpha)*tan(beta)
double GetTanOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	// function assumes the user imputed the angles in radians

	double OutResult = 0.0;
	OutResult = (std::tan(AlphaAngle) + std::tan(BetaAngle)) / 
		(1.0 - (std::tan(AlphaAngle)*std::tan(BetaAngle)));

	return OutResult;
}

// tan(alpha - beta) = tan(alpha) - tan(beta) / 1 + tan(alpha)*tan(beta)
double GetTanOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle)
{
	double OutResult = 0.0;
	OutResult = (std::tan(AlphaAngle) - std::tan(BetaAngle)) /
		(1.0 + (std::tan(AlphaAngle)*std::tan(BetaAngle)));

	return OutResult;
}

// Finding Multiple Sums and Differences of Angles
/*
Given sin α = 3/5, 0<α<π2, cos β = −5/13, π<β<3π2,  find
	sin(α + β)sin(α + β)
	cos(α + β)cos(α + β)
	tan(α + β)tan(α + β)
	tan(α−β)
*/
void FindMultipleSumsAndDifferencesOfAngles(std::pair<double,double> InSineOfAlphaAngle,
						std::pair<double,double> InCosineOfBetaAngle)
{

	// alpha angle is in the first quadrant between 0 and pi/2
	// TODO: add a check for this?
	// triangle #1
	double SideOppositeAlpha = InSineOfAlphaAngle.first;
	double HypotenuseOfAlpha = InSineOfAlphaAngle.second;

	double OutBase = 0.0;
	FindMissingSidePythagoreanTheorm(OutBase, SideOppositeAlpha, HypotenuseOfAlpha);

	// triangle #2
	// beta angle is between pi and 3pi/2 which is in the third quadrant
	double SideAdjacentToBeta = InCosineOfBetaAngle.first;
	double HypotenuseOfBeta = InCosineOfBetaAngle.second;

	double OutHeight = 0.0;
	FindMissingSidePythagoreanTheorm(SideAdjacentToBeta, OutHeight, HypotenuseOfBeta);
	// Find missing side returns a positive result.
	// since this result is in the third quadrant make it negative (the square root is pos/neg)
	OutHeight = OutHeight*(-1);

	// Next step is finding the cosine of alpha and the sine of beta
	// cosine of alpha is the adjacent side over the hypotenuse
	// cos(alpha) = OutBase/HypotenuseOfAlpha
	double CosAngleOfAlpha = OutBase / HypotenuseOfAlpha;
	double SinAngleOfAlpha = SideOppositeAlpha / HypotenuseOfAlpha;

	// sine of beta is the opposite side over the hypotenuse
	// sin(beta) = OutHeight/HypotenuseOfBeta
	double  SinAngleOfBeta = OutHeight / HypotenuseOfBeta;
	double  CosAngleOfBeta = SideAdjacentToBeta / HypotenuseOfBeta;

	// Sin(AlphaAngle + BetaAngle) 
	double OutSinSumRes = SinAngleOfAlpha*CosAngleOfBeta + CosAngleOfAlpha*SinAngleOfBeta;
	std::cout << "Sum of Sin(Alpha + Beta) =  " << OutSinSumRes << std::endl;
	std::cout << std::endl;


	// Cos(AlphaAngle + BetaAngle) 
	double OutCosSumRes = CosAngleOfAlpha*CosAngleOfBeta - SinAngleOfAlpha*SinAngleOfBeta;
	std::cout << "Sum of Cos(Alpha + Beta) =  " << OutCosSumRes << std::endl;
	std::cout << std::endl;


	// tan(alpha + beta)
	double TanAngleAlpha = SinAngleOfAlpha / CosAngleOfAlpha;
	double TanAngleBeta = SinAngleOfBeta / CosAngleOfBeta;

	double TanAlphaPlusBeta = (TanAngleAlpha + TanAngleBeta) /
								(1 - (TanAngleAlpha * TanAngleBeta));

	std::cout << "Sum of Tan(Alpha + Beta) =  " << TanAlphaPlusBeta << std::endl;
	std::cout << std::endl;

	// To find tan(α−β), we have the values we need. We can substitute them in and evaluate.
	double TanAlphaMinusBeta = (TanAngleAlpha - TanAngleBeta) /
		(1 + (TanAngleAlpha * TanAngleBeta));

	std::cout << "Sum of Tan(Alpha - Beta) =  " << TanAlphaMinusBeta << std::endl;
	std::cout << std::endl;
}


// 68.622153857191046346385474553595
// 64.031242374328486864882176746218

// sinθ = 1/cscθ
// cosθ = 1/secθ
// tanθ = 1/cotθ
// csc θ = 1/sin θ
// sec θ = 1/cos θ
// cot  θ = 1/tan θ

double InSinOutCosCofunction(double InSinAngle)
{
	// input is a sin angle in radians

	// return the cosine cofunction angle

	double OutCosAngle = 0.0;
	OutCosAngle = (M_PI / 2) - InSinAngle;
	return OutCosAngle;
}

double InTanOutCotCofunction(double InTanAngle)
{
	// input is a sin angle in radians

	// return the angle in terms of cot

	double OutCotAngle = 0.0;
	OutCotAngle = (M_PI / 2) - InTanAngle;
	return OutCotAngle;
}

double InSecOutCscCofunction(double InSecAngle)
{
	// input is a sec angle in radians

	// return the csc cofunction angle

	double OutCscAngle = 0.0;
	OutCscAngle = (M_PI / 2) - InSecAngle;
	return OutCscAngle;
}

double InCosOutSinCofunction(double InCosAngle)
{
	// input is a cos angle in radians

	// return the angle in terms of sin

	double OutSinAngle = 0.0;
	OutSinAngle = (M_PI / 2) - InCosAngle;
	return OutSinAngle;
}

double InCotOutTanCofunction(double InCotAngle)
{
	// input is a cot angle in radians

	// return the tan cofunction angle

	double OutTanAngle = 0.0;
	OutTanAngle = (M_PI / 2) - InCotAngle;
	return OutTanAngle;
}

double InCscOutSecCofunction(double InCscAngle)
{
	// input is a csc angle in radians

	// return the angle in terms of sec

	double OutSecAngle = 0.0;
	OutSecAngle = (M_PI / 2) - InCscAngle;
	return OutSecAngle;
}



// use The double-angle formulas given tan of an angle
double InTanOfAnAngleOutDoubleAngleSin(std::pair<double,double> InTanOfAngle, 
	const DoubleAngleReturnType& TypeToUse)
{
	// quad 2 so cos is negative and sin is positive
	// so if input is -3/4
	// and tan is sin/cos
	

	// Functions Assumes tan is inputed in quadrant 2;
	double SideOppositeOfAngle = InTanOfAngle.first;
	double SideAdjacentAngle = InTanOfAngle.second;


	double HypotenuseOfAngle = 0.0;
	FindMissingSidePythagoreanTheorm(SideAdjacentAngle, SideOppositeOfAngle, HypotenuseOfAngle);
	
	// Let’s begin by writing the double-angle formula for sine.
	//sin(2θ) = 2 sin θ cos θ
	double SinOfAngle = SideOppositeOfAngle / HypotenuseOfAngle;
	double CosOfAngle = SideAdjacentAngle / HypotenuseOfAngle;

	if (TypeToUse == DoubleAngleReturnType::SIN)
	{
		return (2 * (SinOfAngle * CosOfAngle));
	}
	else if (TypeToUse == DoubleAngleReturnType::COS)
	{
		// cos(2θ)=cos^2(θ)−sin^2(θ)
		return (std::pow(CosOfAngle, 2) - (std::pow(SinOfAngle, 2)));

	}
	else // == tan
	{
		//tan(2θ) = 2 tan θ/ 1−tan^2(θ)
		// taken from the input
		double TanAngle = SideOppositeOfAngle / SideAdjacentAngle;

		return ((2 * TanAngle) / (1 - std::pow(TanAngle, 2)));
	}


}

/* for now.. these suck and dont work.
// The Trig Reduction Formulas
double SinSquaredReductionFormula(const double& InAngle)
{
	double OutReducedAngle{ 0.0 };

	// Input is a Sin^2(InAngle)
	// Output is the reduced version

	// cos(2θ)=cos^2θ−sin^2θ
	OutReducedAngle = ((1 - (std::pow(InAngle, 2) - std::pow(InAngle, 2))) / 2);
	std::cout << OutReducedAngle << std::endl;
	std::cout << std::sin(std::pow(InAngle, 2)) << std::endl;
	return OutReducedAngle;

}

double CosSquaredReductionFormula(const double InAngle)
{
	double OutReducedAngle{ 0.0 };
	return OutReducedAngle;

}

double TanSquaredReductionFormula(const double InAngle)
{
	double OutReducedAngle{ 0.0 };
	return OutReducedAngle;

}
*/

/*
double SinHalfAngleFormula(const double& AngleDegrees)
{
	// We do our calculations in Radians convert.
	double AngleRadians = DegreesToRadians(AngleDegrees);
	
	// sin(AngleDegrees/2) = +- sqrt(1-cos(a)/2)
	double OutVal{ 0.0 };
	OutVal = std::sqrt((1 - std::cos(AngleRadians)) / 2.0);


	if (AngleRadians < (M_PI/2) && AngleRadians > 0.0) // is in the first quadrant
	{
		// keep the positive result
		return OutVal;
	}
	
	// temp for later
	OutVal = OutVal * (-1);
	return OutVal;
}*/



double SinHalfAngleFormula2(const double& InSin, const std::pair<double, double>& InCos)
{
	// sin(AngleDegrees/2) = +- sqrt(1-cos(a)/2)
	double Cosine = InCos.first / InCos.second;
	double OutVal{ 0.0 };
	OutVal = std::sqrt((1 - Cosine) / 2.0);

	/*
	if (InSin < (M_PI / 2) && InSin > 0.0) // is in the first quadrant
	{
		// keep the positive result
		return OutVal;
	}*/

	// temp for later
	//OutVal = OutVal * (-1);
	return OutVal;
}

double CosHalfAngleFormula(const double& CosAngle , const std::pair<double, double>& InCos)
{
	// sin(AngleDegrees/2) = +- sqrt(1-cos(a)/2)
	double Cosine = InCos.first / InCos.second;
	double OutVal{ 0.0 };
	OutVal = std::sqrt((1 + Cosine) / 2.0);
	
	/*
	if (CosAngle < (M_PI / 2) && CosAngle > 0.0) // is in the first quadrant
	{
		// keep the positive result
		return OutVal;
	}*/

	// temp for later
	//OutVal = OutVal * (-1);
	return OutVal;
}

double TanHalfAngleFormula(const double& TanAngle, const std::pair<double, double>& InCos)
{
	// sin(AngleDegrees/2) = +- sqrt(1-cos(a)/2)
	double Cosine = InCos.first / InCos.second;
	double OutVal{ 0.0 };
	OutVal = std::sqrt((1 - Cosine) / (1 + Cosine));

	/*
	if ((TanAngle < (M_PI / 2) && TanAngle > 0.0)) // is in the first quadrant
	{
		// keep the positive result
		return OutVal;
	}*/

	// temp for later
	//OutVal = OutVal * (1.0);
	return OutVal;
}

// Before we start, we must remember that if  α is in quadrant III,
// then 180° < α < 270°,   so 180°/2 < α/2 < 270°/2. 
// This means that the terminal side of α/2   is in quadrant II, since 90°< α/2 < 135°.
double Problem1(std::pair<double, double> TanAngle, const Quadrants& InQuad, 
				const HalfAngleFormula& FormulaToUse)
{

	double TriangleHeight = TanAngle.first;
	double TriangleBase = TanAngle.second;

	double Hypotenuse{ 0.0 };
	// remember to set the missing part to 0.0 
	FindMissingSidePythagoreanTheorm(TriangleBase, TriangleHeight, Hypotenuse);

	
	double SinFuncOfAngle = (TriangleHeight / Hypotenuse);
	double CosFuncOfAngle = (TriangleBase / Hypotenuse);
	double TanFuncOfAngle = (TriangleHeight / TriangleBase);

	std::pair<double, double> TempCos{ TriangleBase, Hypotenuse };
	
	double HalfAngleResult{ 0.0 };

	if (FormulaToUse == HalfAngleFormula::Sin)
	{
		HalfAngleResult = SinHalfAngleFormula2(SinFuncOfAngle, TempCos);

		if (InQuad == Quadrants::Three)
		{
			// This means that the terminal side is in quadrant 2 since 90 < a/2 < 135
			if (HalfAngleResult > 0)
			{
				// keep the positive result since sin is positive in the second quadrant
				return HalfAngleResult;
			}
		}
	}
	else if (FormulaToUse == HalfAngleFormula::Cos)
	{
		HalfAngleResult = CosHalfAngleFormula(CosFuncOfAngle, TempCos);
		
		if (InQuad == Quadrants::Three)
		{
			// This means that the terminal side is in quadrant 2 since 90 < a/2 < 135
			// keep the negative result, since cos is negative in quadrant 2
			if (HalfAngleResult > 0)
			{
				return HalfAngleResult * (-1);
			}
		}
	}
	else if (FormulaToUse == HalfAngleFormula::Tan)
	{
		HalfAngleResult = TanHalfAngleFormula(TanFuncOfAngle, TempCos);

		if (InQuad == Quadrants::Three)
		{
			// This means that the terminal side is in quadrant 2 since 90 < a/2 < 135
			//  use negative result, since tan is negative in quadrant 2
			if (HalfAngleResult > 0)
			{
				return HalfAngleResult * (-1);
			}
		}
	}
	

	// TODO: temp
	return 0.0;
}



struct LawOfSinesInfoAAS
{
	LawOfSinesInfoAAS(const double Alpha,
					const double Beta,
					const double Gamma,
					const GivenSide& SideType,
					const double SideAmount)
	{
		m_Alpha = Alpha;
		m_Beta = Beta;
		m_Gamma = Gamma;
		

		switch (SideType)
		{
			case GivenSide::A:
			{
				SetSideA(SideAmount);
				break;
			}
			case GivenSide::B:
			{
				SetSideB(SideAmount);
				break;
			}
			case GivenSide::C:
			{
				SetSideC(SideAmount);
				break;
			}
		}
		
		if (m_Alpha == 0)
		{
			SetAlphaAngle();
		}
		else if (m_Beta == 0)
		{
			SetBetaAngle();
		}
		else
		{
			SetGammaAngle();
		}

		CalculateMissingSides(SideType);

	}

	void SetSideA(const double Length) { m_A = Length; }
	void SetSideB(const double Length) { m_B = Length; }
	void SetSideC(const double Length) { m_C = Length; }

	double GetAlphaAngle() const { return m_Alpha; }
	double GetGammaAngle() const { return m_Gamma; }
	double GetBetaAngle() const { return m_Beta; }

	double GetSideA() const { return m_A; }
	double GetSideB() const { return m_B; }
	double GetSideC() const { return m_C; }

	void PrintSideA() const
	{
		std::cout << std::setprecision(2) << "Side A: " << GetSideA() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideB() const
	{
		std::cout << std::setprecision(2) << "Side B: " << GetSideB() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideC() const
	{
		std::cout << std::setprecision(2) << "Side C: " << GetSideC() << std::endl;
		std::cout << std::fixed;
	}

	void PrintAllSidesNewLines() const
	{
		PrintSideA();
		PrintSideB();
		PrintSideC();
	}

	void PrintAngleAlpha() const
	{
		std::cout << "Alpha: " << GetAlphaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleBeta() const
	{
		std::cout << "Beta: " << GetBetaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleGamma() const
	{
		std::cout << "Gamma: " << GetGammaAngle() << " Degrees" << std::endl;
	}

	void PrintAllAnglesNewLines() const
	{
		PrintAngleAlpha();
		PrintAngleBeta();
		PrintAngleGamma();
	}

	void PrintAllSidesAndAngles() const
	{
		PrintAllSidesNewLines();
		PrintAllAnglesNewLines();
	}


	void CalculateMissingSides(const GivenSide& Side)
	{
		switch (Side)
		{
			case GivenSide::A:
			{
				// To find an unknown side, we need to know the corresponding angle and a known ratio. 
				// We know that Angle Alpha = 50 degrees and its side a = 10
				// we can now use the law of sines proprotion
				// Sin(AlphaAngle)/SideA = Sin(GammaAngle)/C
				// then solve for C
		
				SetSideC((std::sin(DegreesToRadians(GetGammaAngle()))*GetSideA()) /
						std::sin(DegreesToRadians(GetAlphaAngle())));

				SetSideB((std::sin(DegreesToRadians(GetBetaAngle()))*GetSideA()) /
					std::sin(DegreesToRadians(GetAlphaAngle())));

				break;
			}
			case GivenSide::B:
			{
				SetSideC((std::sin(DegreesToRadians(GetGammaAngle()))*GetSideB()) /
					std::sin(DegreesToRadians(GetBetaAngle())));

				SetSideA((std::sin(DegreesToRadians(GetAlphaAngle()))*GetSideB()) /
					std::sin(DegreesToRadians(GetBetaAngle())));

				break;
			}
			case GivenSide::C:
			{
				SetSideA((std::sin(DegreesToRadians(GetAlphaAngle()))*GetSideC()) /
					std::sin(DegreesToRadians(GetGammaAngle())));

				SetSideB((std::sin(DegreesToRadians(GetBetaAngle()))*GetSideC()) /
					std::sin(DegreesToRadians(GetGammaAngle())));

				break;
			}
		}
	}

private:

	double m_Alpha{ 0.0 };
	double m_Beta{ 0.0 };
	double m_Gamma{ 0.0 };

	double m_A{ 0.0 };
	double m_B{ 0.0 };
	double m_C{ 0.0 };


	void SetAlphaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Alpha = 180.0 - GetBetaAngle() - GetGammaAngle();
	}

	void SetBetaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Beta = 180.0 - GetAlphaAngle() - GetGammaAngle();
	}

	void SetGammaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Gamma = 180.0 - GetAlphaAngle() - GetBetaAngle();
	}
	
};


/*
// AAS (angle-angle-side) We know the measurements of two angles and
// a side that is not between the known angles. 
// Solving for Two Unknown Sides and Angle of an AAS Triangle
LawOfSinesInfoAAS LawOfSinesAngleAngleSide(const double& AlphaAngle, const double& BetaAngle, const double& GammaAngle,
								const GivenSide& Side, const double& SideAmount)
{
	LawOfSinesInfoAAS OutInfo(AlphaAngle, BetaAngle, GammaAngle, Side, SideAmount);

	return OutInfo;
}*/


struct SecondAnswerTriangleData
{
	typedef double Side;
	typedef double Angle;

	Side m_A;
	Side m_B;
	Side m_C;

	Angle m_Alpha;
	Angle m_Beta;
	Angle m_Gamma;

};


class LawOfSinesInfoSSA
{
public:
	typedef double Side;
	typedef double Angle;

	LawOfSinesInfoSSA(const Side& FirstSide, const Side& SecondSide,
		const GivenSides& InSides, const Angle& AngleAmount, const GivenAngle& InAngle)
	{

		switch (InSides)
		{
			case GivenSides::AB:
			{
				m_A = FirstSide;
				m_B = SecondSide;
				m_C = 0.0;

				SecondaryData.m_A = FirstSide;
				SecondaryData.m_B = SecondSide;
				SecondaryData.m_C = 0.0;
			

				if (m_B > m_A)
				{
					//m_IsObtuse = true;
				
				}
				else
				{
					m_IsObtuse = true;
					//throw std::exception("Error: You must enter an obtuse triangle to use this functionality");
				}

				break;
			}
			case GivenSides::AC:
			{
				m_A = FirstSide;
				m_C = SecondSide;
				m_B = 0.0;

				break;
			}
			case GivenSides::BC:
			{
				m_A = 0.0;
				m_B = FirstSide;
				m_C = SecondSide;

				SecondaryData.m_A = 0.0;
				SecondaryData.m_B = FirstSide;
				SecondaryData.m_C = SecondSide;

				if (m_C > m_B)
				{
					m_IsObtuse = true;
				}
				else
				{
					throw std::exception("Error: You must enter an obtuse triangle to use this functionality");
				}

				break;
			}
		}

		switch (InAngle)
		{
			case GivenAngle::Alpha:
			{
				m_Alpha = AngleAmount;
				SecondaryData.m_Alpha = AngleAmount;

				break;
			}
			case GivenAngle::Beta:
			{
				m_Beta = AngleAmount;
				break;
			}
			case GivenAngle::Gamma:
			{
				m_Gamma = AngleAmount;
				SecondaryData.m_Gamma = AngleAmount;
				break;
			}
		}



		CalculateMissingAngles(InAngle);
		CalculateMissingSide(InSides);

	}

	void CalculateMissingAngles(const GivenAngle& InAngle)
	{

		switch (InAngle)
		{
			case GivenAngle::Alpha:
			{
				double LHS = (m_B * std::sin(DegreesToRadians(m_Alpha))) / m_A;
				

				m_Beta = std::asin(LHS);
				m_Beta = RadianToDegrees(m_Beta);
				
				SecondaryData.m_Beta = m_Beta;
	
				
				if (m_IsObtuse == true)
				{
					m_Beta = 180.0 - m_Beta;
					m_Gamma = 180.0 - (m_Alpha)-(m_Beta);
				}

				SecondaryData.m_Gamma = 180.0 - SecondaryData.m_Alpha - SecondaryData.m_Beta;

				break;
			}
			case GivenAngle::Beta:
			{

				//m_Beta = AngleAmount;


				break;

			}
			case GivenAngle::Gamma:
			{

				double LHS = (m_B * std::sin(DegreesToRadians(m_Gamma))) / m_C;


				m_Beta = std::asin(LHS);
				m_Beta = RadianToDegrees(m_Beta);

				SecondaryData.m_Beta = m_Beta;


				if (m_IsObtuse == true)
				{
					m_Beta = 180.0 - m_Beta;
					m_Alpha = 180.0 - (m_Gamma)-(m_Beta);
				}

				SecondaryData.m_Alpha = 180.0 - SecondaryData.m_Gamma - SecondaryData.m_Beta;



				break;
			}
		}

	}

	void CalculateMissingSide(const GivenSides& InSides)
	{
		switch (InSides)
		{
			case GivenSides::AB:
			{

				m_C = (m_A * std::sin(DegreesToRadians(m_Gamma))) /
					std::sin(DegreesToRadians(m_Alpha));

				SecondaryData.m_C = (SecondaryData.m_A * std::sin(DegreesToRadians(SecondaryData.m_Gamma)) /
					std::sin(DegreesToRadians(SecondaryData.m_Alpha)));

				/*
				if (m_IsObtuse)
				{
					m_C = 
				}*/


				break;
			}
			case GivenSides::AC:
			{
				//m_A = FirstSide;
			//	m_C = SecondSide;
				//m_B = 0.0;

				break;
			}
			case GivenSides::BC:
			{

				m_A = (m_C * std::sin(DegreesToRadians(m_Alpha))) /
					std::sin(DegreesToRadians(m_Gamma));

				SecondaryData.m_A = (SecondaryData.m_C * std::sin(DegreesToRadians(SecondaryData.m_Alpha)) /
					std::sin(DegreesToRadians(SecondaryData.m_Gamma)));

				//m_A = 0.0;
				//m_B = FirstSide;
				//m_C = SecondSide;

				break;
			}
		}
	}

	double GetAlphaAngle() const { return m_Alpha; }
	double GetGammaAngle() const { return m_Gamma; }
	double GetBetaAngle() const { return m_Beta; }

	double GetSideA() const { return m_A; }
	double GetSideB() const { return m_B; }
	double GetSideC() const { return m_C; }

	void PrintSideA() const
	{
		std::cout << std::setprecision(2) << "Side A: " << GetSideA() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideB() const
	{
		std::cout << std::setprecision(2) << "Side B: " << GetSideB() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideC() const
	{
		std::cout << std::setprecision(2) << "Side C: " << GetSideC() << std::endl;
		std::cout << std::fixed;
	}

	void PrintAllSidesNewLines() const
	{
		PrintSideA();
		PrintSideB();
		PrintSideC();
	}

	void PrintAngleAlpha() const
	{
		std::cout << "Alpha: " << GetAlphaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleBeta() const
	{
		std::cout << "Beta: " << GetBetaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleGamma() const
	{
		std::cout << "Gamma: " << GetGammaAngle() << " Degrees" << std::endl;
	}

	void PrintAllAnglesNewLines() const
	{
		PrintAngleAlpha();
		PrintAngleBeta();
		PrintAngleGamma();
	}

	void PrintAllSidesAndAngles() const
	{
		std::cout << "If there is any negative angles in these triangles.\nThen that is a failed triangle (ignore the results for it) " << std::endl << std::endl;

		PrintAllSidesNewLines();
		PrintAllAnglesNewLines();
		PrintOtherAnswerTriangle();
	}

	void PrintOtherAnswerTriangle() const
	{
		if (m_IsObtuse)
		{
			// next
			std::cout << std::endl;
			std::cout << "Printing next possible triangle. " << std::endl << std::endl;

			std::cout << std::setprecision(2) << "Side A: " << SecondaryData.m_A << std::endl;
			std::cout << std::fixed;

			std::cout << std::setprecision(2) << "Side B: " << SecondaryData.m_B << std::endl;
			std::cout << std::fixed;

			std::cout << std::setprecision(2) << "Side C: " << SecondaryData.m_C << std::endl;
			std::cout << std::fixed;

			std::cout << "Alpha: " << SecondaryData.m_Alpha << " Degrees" << std::endl;

			std::cout << "Beta: " << SecondaryData.m_Beta << " Degrees" << std::endl;

			std::cout << "Gamma: " << SecondaryData.m_Gamma << " Degrees" << std::endl;
		}
	}

private:
	Side m_A;
	Side m_B;
	Side m_C;

	Angle m_Alpha;
	Angle m_Beta;
	Angle m_Gamma;

	bool m_IsObtuse = false;
	//bool m_NotPossibleResult = false;

	SecondAnswerTriangleData SecondaryData;
};

void FindAltitudeOfObject(const double& AngleA, const double& AngleB,
							const double& DistanceBetweenC)
{
	//To find the elevation of the aircraft, 
	// we first find the distance from one station to the aircraft,
	const double AngleC = 180.0 - AngleA - AngleB;
	
	const double A = (DistanceBetweenC * std::sin(DegreesToRadians(AngleA)) /
						std::sin(DegreesToRadians(AngleC)));

	std::cout << "The distance from one station to the aircraft is about " <<
			A << " units." << std::endl;

	// Now that we know a, we can use right triangle relationships to solve for  h
	const double Altitude = (A * (std::sin(DegreesToRadians(AngleB))));

	std::cout << "The aircraft is at an altitude of approximately " << Altitude << " units."
				<< std::endl;
}

double LawOfSinesFindAngleAlpha(const double& SideA, const double& BetaAngle, const double& SideB)
{
	double BetaAngRadians = DegreesToRadians(BetaAngle);

	double OutResult = (SideA * std::sin(BetaAngRadians)) / SideB;

	// find the inverse sine
	OutResult = std::asin(OutResult);

	// return the result as degrees
	OutResult = RadianToDegrees(OutResult);

	return OutResult;

}

double LawOfSinesFindAngleBeta(const double& SideA, const double& AlphaAngle, const double& SideB)
{
	double AlphaAngRadians = DegreesToRadians(AlphaAngle);

	double OutResult = (SideB * std::sin(AlphaAngRadians)) / SideA;

	// find the inverse sine
	OutResult = std::asin(OutResult);

	// return the result as degrees
	OutResult = RadianToDegrees(OutResult);

	return OutResult;

}

double LawOfSinesFindAngleGamma(const double& SideA, const double& AlphaAngle, const double& SideC)
{
	double AlphaAngRadians = DegreesToRadians(AlphaAngle);

	double OutResult = (SideC * std::sin(AlphaAngRadians)) / SideA;

	// find the inverse sine
	OutResult = std::asin(OutResult);

	// return the result as degrees
	OutResult = RadianToDegrees(OutResult);

	return OutResult;

}

void ConvertPolarToRectangle(const double& InRadius, const double InAngle, double& OutX, double& OutY)
{
	//double AngleToRadian = DegreesToRadians(InAngle);

	OutX = (std::cos(InAngle)) * InRadius;
	OutY = (std::sin(InAngle)) * InRadius;

	return;
	
}

void PrintCordinates(const double& x, const double& y)
{
	std::cout << std::fixed << "Rectangle Cords = (" << x << ',' << y << ")" << std::endl;
}

void PrintCordinatesSTDPair(const std::pair<double, double>& Cords)
{
	std::cout << std::fixed  << "Polar Cords = " <<
		"(" << Cords.first << "," << Cords.second << ")" << std::endl;
}

// Convert the rectangular coordinates (3, 3) to polar coordinates.
std::pair<double, double> ConvertToPolarFromRectangular(const double& x, const double &y)
{
	
	std::pair<double, double> OutPolarCords(0, 0);

	double AngleInRadians = std::atan(y / x);
	OutPolarCords.second = AngleInRadians;

	
	double Radius = std::sqrt(std::pow(x, 2) + std::pow(y, 2));

	if (AngleInRadians < 0)
	{
		Radius = Radius * (-1);
	}

	OutPolarCords.first = Radius;

	return OutPolarCords;
}

template <typename T>
void PrintAllPolarPoints(const T& CONT)
{
	for (const auto& pair : CONT)
	{
		std::cout << std::fixed << "(" << pair.first << ", " << pair.second << ")" << std::endl;
	}
}



struct PairSort
{
	template <typename T>
	bool operator()(T a, T b) const
	{
		return a.first > b.first;
	}
};


template <typename T>
void PrintMaxPairForPolarContainer(const T& CONT)
{
	PairSort sort;

	auto cont = CONT;
	std::sort(cont.begin(), cont.end(), sort);

	std::cout << std::fixed << "(" << cont[0].first << ", " << cont[0].second << ")" << std::endl;
}


template <typename T>
T GetPairsOfTheZerosForPolarContainer(T& CONT)
{
	T OutCONT;

	for (T::iterator it = CONT.begin(); it != CONT.end(); ++it)
	{
		if ((*it).first <= 0.1 && (*it).first >= -0.1)
		{
			OutCONT.push_back(*it);
		}
	}

	return OutCONT;

}

// Recall that, to find the zeros of polynomial functions,
// we set the equation equal to zero and then solve for x.
// We use the same process for polar equations. Set r=0, and solve for θ.

std::vector<std::pair<double,double>> ReturnMaxMinPointsSin(const double& PoleDist)
{
	// the maximum value of a sine function is when sin(pi/2) = 1 * DistanceBetweenPole
	std::vector<std::pair<double, double>> OutPairVector;
	
	// show the function we're using for clarification
	// Input = Theta  ||| r = PoleDist*Sin(Theta) ||| Output = r
	
	std::pair<double, double> PairOne(0, 0);
	// when theta = 0

	// r = PoleDist*Sin(0) NOTE: Sin(0) == 0 so this Input Output is (0,0)
	PairOne.first = 0;
	PairOne.second = 0;

	OutPairVector.push_back(PairOne);

	std::pair<double, double> PairTwo(0, 0);
	PairTwo.second = M_PI / 6.0;

	// r = PoleDist*Sin(M_PI/6.0)
	PairTwo.first = PoleDist * (std::sin(M_PI / 6.0));

	OutPairVector.push_back(PairTwo);

	std::pair<double, double> PairThree(0, M_PI / 3.0);
	PairThree.first = PoleDist * (std::sin(M_PI / 3.0));

	OutPairVector.push_back(PairThree);

	// This is the max input for sin(pi/2) = 1
	std::pair<double, double> PairFour(0, M_PI / 2.0);
	PairFour.first = PoleDist * (std::sin(M_PI / 2.0));

	OutPairVector.push_back(PairFour);

	std::pair<double, double> PairFive(0, TwoM_PI / 3.0);
	PairFive.first = PoleDist * (std::sin(TwoM_PI / 3.0));

	OutPairVector.push_back(PairFive);

	std::pair<double, double> PairSix(0, ((5 * M_PI) / 6.0));
	PairSix.first = PoleDist * (std::sin(((5 * M_PI) / 6.0)));

	OutPairVector.push_back(PairSix);

	std::pair<double, double> PairSeven(0, M_PI);
	PairSeven.first = PoleDist * (std::sin(M_PI));

	OutPairVector.push_back(PairSeven);

	return OutPairVector;
}

std::vector<std::pair<double, double>> ReturnMaxMinPointsCos(const double& PoleDist,
	const double& PlusMinusVar = 0)
{
	// the maximum value of a sine function is when sin(pi/2) = 1 * DistanceBetweenPole
	std::vector<std::pair<double, double>> OutPairVector;

	// show the function we're using for clarification
	// Input = Theta  ||| r = PoleDist*cos(Theta) ||| Output = r

	std::pair<double, double> PairOne(0, 0);
	// when theta = 0

	PairOne.first = (PoleDist * (std::cos(0))) + PlusMinusVar;
	PairOne.second = 0;

	OutPairVector.push_back(PairOne);

	std::pair<double, double> PairTwo(0, 0);
	PairTwo.second = M_PI / 6.0;

	// r = PoleDist*Sin(M_PI/6.0)
	PairTwo.first = (PoleDist * (std::cos(M_PI / 6.0))) + PlusMinusVar;

	OutPairVector.push_back(PairTwo);

	std::pair<double, double> PairThree(0, M_PI / 3.0);
	PairThree.first = (PoleDist * (std::cos(M_PI / 3.0))) + PlusMinusVar;

	OutPairVector.push_back(PairThree);

	// This is the max input for sin(pi/2) = 1
	std::pair<double, double> PairFour(0, M_PI / 2.0);
	PairFour.first = (PoleDist * (std::cos(M_PI / 2.0))) + PlusMinusVar;

	OutPairVector.push_back(PairFour);

	std::pair<double, double> PairFive(0, TwoM_PI / 3.0);
	PairFive.first = (PoleDist * (std::cos(TwoM_PI / 3.0))) + PlusMinusVar;

	OutPairVector.push_back(PairFive);

	std::pair<double, double> PairSix(0, ((5 * M_PI) / 6.0));
	PairSix.first = (PoleDist * (std::cos(((5 * M_PI) / 6.0)))) + PlusMinusVar;

	OutPairVector.push_back(PairSix);

	std::pair<double, double> PairSeven(0, M_PI);
	PairSeven.first = (PoleDist * (std::cos(M_PI))) + PlusMinusVar;

	OutPairVector.push_back(PairSeven);

	return OutPairVector;
}


class PolarFuncSymmetryResult
{
public:

	void PrintTestResultsNewLines() const
	{
		if (GetFirstTestResult())
		{
			std::cout << "Symmetry around pi/2 Test Passes" << std::endl;
		}
		if (GetSecondTestResult())
		{
			std::cout << "Polar Axis Symmetry Test Passes" << std::endl;
		}
		if (GetThirdTestResult())
		{
			std::cout << "Respect to pole Symmetry Test Passes" << std::endl;
		}
	}

	void PassFirstTest()
	{
		m_PI_2Symmetry = true;
	}
	void PassSecondTest()
	{
		m_PolarAxisSymmetry = true;
	}
	void PassThirdTest()
	{
		m_RespectToPoleSymmetry = true;
	}

	bool GetFirstTestResult() const { return m_PI_2Symmetry; }
	bool GetSecondTestResult() const { return m_PolarAxisSymmetry; }
	bool GetThirdTestResult() const { return m_RespectToPoleSymmetry; }


private:
	bool m_PI_2Symmetry = false;
	bool m_PolarAxisSymmetry = false;
	bool m_RespectToPoleSymmetry = false;
};

PolarFuncSymmetryResult PolarSymmentryTest(//const double& r,
	const double& Theta,
	const TrigFunctions& Formula,
	const double& FuncMultiplier,
	const double& PlusMinusVar = 0)
{
	PolarFuncSymmetryResult OutResult;

	if (Formula == TrigFunctions::SIN)
	{
		// Replacing (r,θ) with (−r,−θ)
		double r = (FuncMultiplier * (std::sin(Theta)) + PlusMinusVar);
		double NegR = (-r);
		if (NegR == (FuncMultiplier * (std::sin(-Theta))) + PlusMinusVar)
		{
			OutResult.PassFirstTest();
		}


		//  (r,−θ)  for polar axis symmetry
		if (r == (FuncMultiplier * (std::sin(-Theta))) + PlusMinusVar)
		{
			OutResult.PassSecondTest();
		}

		//  (−r,θ) for symmetry with respect to the pole.
		if (NegR == (FuncMultiplier * (std::sin(Theta))) + PlusMinusVar)
		{
			OutResult.PassThirdTest();
		}
	}

	// do first test using cos TODO: add others 
	if (Formula == TrigFunctions::COS)
	{
		// Replacing (r,θ) with (−r,−θ)
		double r = (FuncMultiplier * (std::cos(Theta)) + PlusMinusVar);
		double NegR = (-r);
		if (NegR == (FuncMultiplier * (std::cos(-Theta))) + PlusMinusVar)
		{
			OutResult.PassFirstTest();
		}


		//  (r,−θ)  for polar axis symmetry
		if (r == (FuncMultiplier * (std::cos(-Theta))) + PlusMinusVar)
		{
			OutResult.PassSecondTest();
		}

		//  (−r,θ) for symmetry with respect to the pole.
		if (NegR == (FuncMultiplier * (std::cos(Theta))) + PlusMinusVar)
		{
			OutResult.PassThirdTest();
		}
	}

	return OutResult;
}

double XAsAFunctionOfTime(const double& timeT)
{
	return timeT;
}

double YAsAFunctionOfTime(const double& timeT)
{
	return std::pow(timeT, 2) - 1.0;
}

double XAsAFunctionOfTimev2(const double& timeT)
{
	return timeT - 3;
}

double YAsAFunctionOfTimev2(const double& timeT)
{
	return (2 * timeT) + 4;
}

double XAsAFunctionOfTimev3(const double& timeT)
{
	return std::pow(timeT, 2) + 1;
}

double YAsAFunctionOfTimev3(const double& timeT)
{
	return 2 + timeT;
}

// An example to print numbers counting down :
void print(int p)
{
	if (p == 0)
		return;
	std::cout << p;
	print(p - 1);
	return;
}

//  An example to print counting up:
void print(double p)
{
	if (p == 0)
		return;
	print(p - 1);
	std::cout << p;
	return;
}


// An example to produce the fibonacci number for a given index in the series:
int Fibonacci(int n)
{
	if (n == 0)
		return 0;
	if (n == 1)
		return 1;
	
	//std::cout << (Fibonacci(n - 2) + Fibonacci(n - 1)) << std::endl;
	return (Fibonacci(n - 2) + Fibonacci(n - 1));
}

int RecursiveTest(int n)
{
	if (n == 1)
		return 1;

	return RecursiveTest((2 * n) - 1);
}

//double TestFunc(const double& x)
//{
//	return (std::pow(x, 2) * (-16)) + 64;
//}

//double TestFunc(const double& x)
//{
//	return std::pow(x, 2) + 1;
//}

//double TestFunc(const double& x)
//{
//	return std::sin(x) / x;
//}

//double TestFunc(const double& x)
//{
//	double Numerator = std::sqrt(x) - 2;
//	double Denominator = x - 4;
//
//	return Numerator / Denominator;
//}

//double TestFunc(const double& x)
//{
//	double Numerator = (1 / x) - 1;
//	double Denominator = x - 1;
//
//	return Numerator / Denominator;
//}

//double TestFunc(const double& x)
//{
//	double Numerator = std::abs(std::pow(x, 2) - 4);
//	double Denominator = x - 2;
//
//	return Numerator / Denominator;
//}

//double TestFunc(const double& x)
//{
//	double Numerator = 1;
//	double Denominator = x;
//
//	return std::sin(Numerator / Denominator);
//}
//
//double TestFunc(const double& x)
//{
//	if (x < 2)
//	{
//		return x + 1;
//	}
//	
//	if (x >= 2)
//	{
//		return std::pow(x, 2) - 4;
//	}
//}


//double TestFunc(const double& x)
//{
//	if (x < 2)
//	{
//		return x + 1;
//	}
//
//	if (x >= 2)
//	{
//		return std::pow(x, 2) - 4;
//	}
//}

//double TestFunc(const double& x)
//{
//	double Numerator = std::abs(std::pow(x, 2) - 4);
//	double Denominator = x - 2;
//
//	return Numerator / Denominator;
//
//}

//double TestFunc(const double& x)
//{
//	return 1 / x;
//}

//double TestFunc(const double& x)
//{
//	return 1 / std::pow(x,2);
//}

double TestFunc(const double& x)
{
	return std::sqrt(x - 3.0);
}

int main()
{
	try
	{


		//Point FirstPoint(1.0, 1.0);
		//Point SecondPoint(2.0, 4.0);
		//Point ThirdPoint(5.0 / 4.0, 25.0 / 16.0);


		////std::cout << GetSlopeOfSecantLineTwoPoints(FirstPoint, SecondPoint) << std::endl;
		//std::cout << GetSlopeOfSecantLineTwoPoints(FirstPoint, ThirdPoint) << std::endl;

		//std::cout << GetSlopeOfSecantLineTwoPoints(SecondPoint, ThirdPoint) << std::endl;

		//double DesiredTime = 0.5;

		//std::cout << GetAverageVelocity(0.49, DesiredTime,TestFunc) << std::endl;
		//std::cout << GetAverageVelocity(0.51, DesiredTime, TestFunc) << std::endl;
		// The instantaneous velocity is somewhere between these two

		//double DesiredTime = 2.0;


		//std::cout << GetAverageVelocity(2.001, DesiredTime, TestFunc) << std::endl;
		//std::cout << GetAverageVelocity(1.999, DesiredTime, TestFunc) << std::endl;


		//double underres = GetAreaUnderCurve(0, 3, TestFunc, true);
		//double overres = GetAreaUnderCurve(0, 3, TestFunc, false);

		//
		//std::cout << underres << std::endl;
		//std::cout << overres << std::endl;

		//Limit TestLimit(TestFunc, 2);




		//Limit TestLimit(TestFunc, 0);

		//LinearFunction TestLinearFunc(4, 2);

		//Limit LinearFuncLimit(TestLinearFunc, -3);

		//QuadraticFunction TestQuad(2.0, -3.0, 1.0);
		//CubicFunction TestCubeic(1.0, 0.0, 0.0, 4.0);

		//RationalFunction TestRational(TestQuad, TestCubeic);

		//Limit TestRationalLimit(TestRational, 2.0);


		//QuadraticFunction TestQuad(2.0, -3.0, 1.0);
		//LinearFunction TestLinear(5.0, 4.0);

		//RationalFunction TestRational(TestQuad, TestLinear);

		//Limit TestRationalLimit(TestRational, 3.0);



		//CubicFunction TestCubic(3.0,0.0, -2.0, 7.0);

		//Limit TestCubicLimit(TestCubic, -2.0);


		//QuadraticFunction TestQuad1(1.0, -3.0);
		//std::vector<double> RealZeros1 = GetZerosQuadraticFormula(TestQuad1);
		//TestQuad1.PrintBasicFunctionInfo();
		//
		//std::cout << std::endl;
		//
		//
		//QuadraticFunction TestQuad2(2.0, -5.0, -3.0);
		//std::vector<double> RealZeros2 = GetZerosQuadraticFormula(TestQuad2);
		//TestQuad2.PrintBasicFunctionInfo();


		//RationalFunction TestRational(TestQuad1, TestQuad2);
		//Limit TestRationalLimit(TestRational, 3.0);

		//QuadraticFunction TestQuad1(1.0, 4.0, 3.0);
		//QuadraticFunction TestQuad2(1.0, 0, -9.0);

		//RationalFunction TestRational(TestQuad1, TestQuad2);
		//Limit TestRationalLimit(TestRational, -3.0);


		//RootFunction RootFuncTest(2.0, 1.0, -2.0, -1.0);
		//LinearFunction LinearFuncTest(1.0, 1.0);

		//RationalFunction RationalFuncTest(RootFuncTest, LinearFuncTest);

		//Limit LimitRootLinearTest(RationalFuncTest, -1.0);

		//RootFunction RootFuncTest(2.0, 1.0, 1.0, -2.0);
		//LinearFunction LinearFuncTest(1.0, -5.0);

		//RationalFunction RationalFuncTest(RootFuncTest, LinearFuncTest);

		//Limit LimitRootLinearTest(RationalFuncTest, 5.0);

		//auto func = TestFunc;

		//Limit TestLimit(TestFunc, 3.0);
		
		QuadraticFunction TestQuadStrRead("(x-3)^2+0");
		Limit TakeLimitQUad(TestQuadStrRead, 2.0);


	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

