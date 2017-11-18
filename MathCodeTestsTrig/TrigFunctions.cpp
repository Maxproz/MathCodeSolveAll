#include "stdafx.h"
#include "TrigFunctions.h"


#include <iostream>
#include <string>
#include <utility>
#include <tuple>
#include <stdexcept>

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


