#pragma once


double LawOfCosinesFindSideA(const double& SideB, const double& SideC, 
								const double& AngleAlpha);

double LawOfCosinesFindSideB(const double& SideA, const double& SideC,
	const double& AngleBeta);

double LawOfCosinesFindSideC(const double& SideA, const double& SideB,
	const double& AngleGamma);

double LawOfCosinesFindAngleAlpha(const double& SideA, const double& SideB,
	const double& SideC);

double LawOfCosinesFindAngleBeta(const double& SideA, const double& SideB,
	const double& SideC);

double LawOfCosinesFindAngleGamma(const double& SideA, const double& SideB,
	const double& SideC);


double LawOfCosinesFindDistanceSAS(const double & SideA, const double & SideB, const double& AngleGamma);
