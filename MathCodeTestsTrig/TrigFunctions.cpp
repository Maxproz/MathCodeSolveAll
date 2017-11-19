#include "stdafx.h"
#include "TrigFunctions.h"


#include <iostream>
#include <string>
#include <utility>
#include <tuple>
#include <stdexcept>

#include "Conics.h"

double GetReferenceAngle(const Quadrants& Quad, const AngleType& Type,
	const double InAngle)
{
	double OutAngleDeg = 0.0;
	double OutAngleRad = 0.0;

	if (Type == AngleType::Degrees)
	{
		switch (Quad)
		{
			case Quadrants::One:
			{
				// an angle in the first quadrant is own ref angle
				OutAngleDeg = InAngle;
				break;
			}
			case Quadrants::Two:
			{
				OutAngleDeg = (abs(180 - InAngle));
				break;
			}
			case Quadrants::Three:
			{
				OutAngleDeg = (abs(180 - InAngle));
				break;
			}
			case Quadrants::Four:
			{
				OutAngleDeg = (360 - InAngle);
				break;
			}
		}

		while (OutAngleDeg < 0)
		{
			OutAngleDeg = (OutAngleDeg + (360));
		}
		while (OutAngleDeg > 360)
		{
			OutAngleDeg = (OutAngleRad - 360);
		}

		return OutAngleDeg;
	}
	else
	{
		switch (Quad)
		{
			case Quadrants::One:
			{
				// an angle in the first quadrant is own ref angle
				OutAngleRad = InAngle;
				break;
			}
			case Quadrants::Two:
			{
				OutAngleRad = (abs(M_PI - InAngle));
				break;
			}
			case Quadrants::Three:
			{
				OutAngleRad = (abs(M_PI - InAngle));
				break;
			}
			case Quadrants::Four:
			{
				OutAngleRad = ((2 * M_PI) - InAngle);
				break;
			}
		}

		while (OutAngleRad < 0)
		{
			OutAngleRad = (OutAngleRad + (2 * M_PI));
		}
		while (OutAngleRad >(2 * M_PI))
		{
			OutAngleRad = (OutAngleRad - (2 * M_PI));
		}

		return OutAngleRad;
	}
}

Degrees RadianToDegrees(const Radians InRad)
{
	Degrees OutDeg{ 0.0 };
	OutDeg = ((InRad * (Degrees)180) / M_PI);
	return OutDeg;
}

Radians DegreesToRadians(const Degrees InDeg)
{
	Radians OutRad{ 0.0 };
	OutRad = ((InDeg * M_PI) / 180);
	return OutRad;
}

// the division problem s/r -- how many r's are there in s?
Radians RadianMeasureOfArcInUnitCircle(const double ArcLength, const double Radius)
{
	Radians OutRad{ 0.0 };
	OutRad = (ArcLength / Radius);
	return OutRad;
}


// circumference of circle
// C = 2(pi)(r) = pi(d)
double FindCircumferenceOfCircle(const double& Radius)
{
	return (2 * M_PI*Radius);
}

// area of circle
// A = (pi)(radius^2)
double FindAreaOfCircle(const double& Radius)
{
	return M_PI*std::pow(Radius, 2);
}

// Area of triangle
// A = (1/2)(b)(h) or (1/2)(b)(c)(sinA) or (1/2)(a)(c)(sinB) or (1/2)(a)(b)(sinC)
double FindAreaOfTriangle(const double& Base, const double& Height)
{
	return ((0.5)*(Base)*(Height));
}

// the following three functions assume the user inputed the angle in degrees
// it is then converted to radians 
double FindAreaOfTriangleAlpha(const double & B, const double & C, const double & AlphaAngle)
{
	double AlphaInRadians(DegreesToRadians(AlphaAngle));
	return (0.5*(B*(C*(std::sin(AlphaInRadians)))));
	
}

double FindAreaOfTriangleBeta(const double & A, const double & C, const double & BetaAngle)
{

	double BetaInRadians(DegreesToRadians(BetaAngle));
	return (0.5*(A*(C*(std::sin(BetaInRadians)))));

}

double FindAreaOfTriangleGamma(const double & A, const double & B, const double & GammaAngle)
{

	double GammaInRadians(DegreesToRadians(GammaAngle));
	return (0.5*(A*(B*(std::sin(GammaInRadians)))));
}



// Pythagorean Theorm
// a squared + bsquared = c squared
// user supplies 2 sides with a value and one side that they set to 0
// return true if success
bool FindMissingSidePythagoreanTheorm(double& base, double& height, double& hypotenuse)
{
	if (base == 0)
	{
		base = std::pow(hypotenuse, 2) - std::pow(height, 2);
		base = sqrt(base);
		return true;
	}
	else if (height == 0)
	{
		height = std::pow(hypotenuse, 2) - std::pow(base, 2);
		height = sqrt(height);
		return true;
	}
	else // c == 0
	{
		hypotenuse = std::pow(base, 2) + std::pow(height, 2);
		hypotenuse = std::sqrt(hypotenuse);
		return true;
	}

	// failed to use algorithm
	return false;
}

// distance formula
// input is 4 points (2 cord pairs) 
// output is a double as the distance between them
double FindDistanceBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2)
{
	// what will be returned after doing calculations
	double OutDistance = 0;

	OutDistance = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
	return OutDistance;
}

// midpoint formula
// input is 4 points (2 cord pairs) 
// output is a double pair as the midpoint cord's between them
std::pair<double, double>
FindMidpointBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2)
{
	std::pair<double, double> MidPoint(0, 0);

	MidPoint.first = ((x1 + x2) / 2);
	MidPoint.second = ((y1 + y2) / 2);

	return MidPoint;
}

double FindSlopeWithTwoPoints(const double& x1, const double& y1, const double& x2, const double& y2)
{
	double OutSlope = 0.0;

	OutSlope = ((y2 - y1) / x2 - x1);
	return OutSlope;
}

void OutputSignsOfTrigFunctionQuadrants(const Quadrants& InQuad)
{
	switch (InQuad)
	{
		case Quadrants::One:
		{
			std::cout << "All functions are positive here.";
			break;
		}
		case Quadrants::Two:
		{
			std::cout << "Sine and its reciprocial cosecant are positive here.";
			break;
		}
		case Quadrants::Three:
		{
			std::cout << "Tangent and its reciprocial cotangent are positive here.";
			break;
		}
		case Quadrants::Four:
		{
			std::cout << "Cosine and its reciprocial secant are positive here.";
			break;
		}
		default:
			break;
	}
}

void OutputTrigEquivelantEquations()
{
	std::cout << "Tan(t) == sin(t)/cos(t)" << std::endl;
	std::cout << "Sec(t) == 1/cos(t)" << std::endl;
	std::cout << "Csc(t) == 1/sin(t)" << std::endl;
	std::cout << "Cot(t) == 1/tan(t) == cos(t)/sin(t) " << std::endl;
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

double SumToProductCosSubtractCosDegrees(const double& AlphaAngle, const double& BetaAngle, const bool& ReturnRadians)
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
void FindMultipleSumsAndDifferencesOfAngles(std::pair<double, double> InSineOfAlphaAngle,
	std::pair<double, double> InCosineOfBetaAngle)
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
double InTanOfAnAngleOutDoubleAngleSin(std::pair<double, double> InTanOfAngle,
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
	std::cout << std::fixed << "Polar Cords = " <<
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





// Recall that, to find the zeros of polynomial functions,
// we set the equation equal to zero and then solve for x.
// We use the same process for polar equations. Set r=0, and solve for θ.

std::vector<std::pair<double, double>> ReturnMaxMinPointsSin(const double& PoleDist)
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
	const double& PlusMinusVar)
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


PolarFuncSymmetryResult PolarSymmentryTest(//const double& r,
	const double& Theta,
	const TrigFunctions& Formula,
	const double& FuncMultiplier,
	const double& PlusMinusVar)
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

double CosHalfAngleFormula(const double& CosAngle, const std::pair<double, double>& InCos)
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
