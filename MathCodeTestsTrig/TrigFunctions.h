#pragma once

#ifndef TRIGFUNC
#define TRIGFUNC

#include <utility>

typedef double Radians;
typedef double Degrees;

constexpr double M_PI = 3.14159;
constexpr auto TwoM_PI = M_PI * 2;



//constexpr auto M_PI = 3.14159265358979323846;

enum class Quadrants
{
	One,
	Two,
	Three,
	Four
};

enum class AngleType
{
	Radians,
	Degrees
};

double GetReferenceAngle(const Quadrants& Quad, const AngleType& Type,
	const double InAngle);

Degrees RadianToDegrees(const Radians InRad);

Radians DegreesToRadians(const Degrees InDeg);

// the division problem s/r -- how many r's are there in s?
Radians RadianMeasureOfArcInUnitCircle(const double ArcLength, const double Radius);

// circumference of circle
// C = 2(pi)(r) = pi(d)
double FindCircumferenceOfCircle(const double& Radius);
// area of circle
// A = (pi)(radius^2)
double FindAreaOfCircle(const double& Radius);

// Area of triangle ( these return square units )
// A = (1/2)(b)(h) or (1/2)(b)(c)(sinA) or (1/2)(a)(c)(sinB) or (1/2)(a)(b)(sinC)
double FindAreaOfTriangle(const double& Base, const double& Height);
double FindAreaOfTriangleAlpha(const double& B, const double& C, const double& AlphaAngle);
double FindAreaOfTriangleBeta(const double& A, const double& C, const double& BetaAngle);
double FindAreaOfTriangleGamma(const double& A, const double& B, const double& GammaAngle);

// Pythagorean Theorm
// a squared + bsquared = c squared
// user supplies 2 sides with a value and one side that they set to 0
// return true if success
bool FindMissingSidePythagoreanTheorm(double& base, double& height, double& hypotenuse);

// distance formula
// input is 4 points (2 cord pairs) 
// output is a double as the distance between them
double FindDistanceBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2);
// midpoint formula
// input is 4 points (2 cord pairs) 
// output is a double pair as the midpoint cord's between them
std::pair<double, double>
FindMidpointBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2);

double FindSlopeWithTwoPoints(const double& x1, const double& y1, const double& x2, const double& y2);

void OutputSignsOfTrigFunctionQuadrants(const Quadrants& InQuad);

void OutputTrigEquivelantEquations();


#endif TRIGFUNC