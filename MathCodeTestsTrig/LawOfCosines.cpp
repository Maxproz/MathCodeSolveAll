#include "stdafx.h"
#include "LawOfCosines.h"

#include "TrigFunctions.h"
#include <cmath>

double LawOfCosinesFindSideA(const double & SideB, const double & SideC, const double & AngleAlpha)
{
	double OutASquared{ 0.0 };
	double AlphaAngRadians = DegreesToRadians(AngleAlpha);
	
	OutASquared = std::pow(SideB, 2) + std::pow(SideC, 2) - (2.0 * (SideB * (SideC * (std::cos(AlphaAngRadians)))));
	double OutSideA = std::sqrt(OutASquared);

	return OutSideA;
}

double LawOfCosinesFindSideB(const double & SideA, const double & SideC, const double & AngleBeta)
{
	double OutBSquared{ 0.0 };
	double BetaAngRadians = DegreesToRadians(AngleBeta);

	OutBSquared = std::pow(SideA, 2) + std::pow(SideC, 2) - (2.0 * (SideA * (SideC * (std::cos(BetaAngRadians)))));
	double OutSideB = std::sqrt(OutBSquared);

	return OutSideB;
}

double LawOfCosinesFindSideC(const double & SideA, const double & SideB, const double & AngleGamma)
{
	double OutCSquared{ 0.0 };
	double GammaAngRadians = DegreesToRadians(AngleGamma);

	OutCSquared = std::pow(SideA, 2) + std::pow(SideB, 2) - (2.0 * (SideA * (SideB * (std::cos(GammaAngRadians)))));
	double OutSideC = std::sqrt(OutCSquared);

	return OutSideC;
}
#include <iostream>

// TODO: These currently have an issue, need fixed later
// the following functions return the cos of an angle using the given sides
double LawOfCosinesFindAngleAlpha(const double & SideA, const double & SideB, const double & SideC)
{
	double FirstHalf = ((std::pow(SideA, 2) - std::pow(SideB, 2) - std::pow(SideC, 2)) /
		((-2.0) * SideB * SideC));
	std::cout << FirstHalf << std::endl;
	// return the angle
	FirstHalf = std::acos(FirstHalf);
	FirstHalf = RadianToDegrees(FirstHalf);

	return FirstHalf;
}

double LawOfCosinesFindDistanceSAS(const double & SideA, const double & SideB, const double& AngleGamma)
{
	double AngleToRad = DegreesToRadians(AngleGamma);

	double LHS = (std::pow(SideA, 2) + std::pow(SideB, 2));
	double LHS2 = 2 * SideA *SideB * std::cos(AngleToRad);
	double LHSSum = LHS - LHS2;

	LHSSum = sqrt(LHSSum);
	
	// return the distance 
	return LHSSum;
}

double LawOfCosinesFindAngleBeta(const double & SideA, const double & SideB, const double & SideC)
{


	return (std::pow(SideA, 2) + std::pow(SideC, 2) - std::pow(SideB, 2)) / ((2.0) * SideA * SideC);
}

double LawOfCosinesFindAngleGamma(const double & SideA, const double & SideB, const double & SideC)
{


	return (std::pow(SideA, 2) + std::pow(SideB, 2) - std::pow(SideC, 2)) / ((2.0) * SideA * SideB);
}
