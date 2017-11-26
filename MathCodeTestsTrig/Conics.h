#pragma once


#include "MiscMathEquations.h"

#include <utility>
#include <string>
#include <map>
#include <sstream>  
#include <iostream>

template<typename Map>
void print_mapp(Map& m);

void PrintDecimalAsFraction(double input);

enum class ConicType
{
	Ellipse,
	Circle,
	Hyperbola,
	Parabola
};

enum class TrigFunctions
{
	NONE,
	SIN,
	COS,
	TAN
};

enum class PlusOrMinusDenom
{
	Plus,
	Minus
};

struct PolarEquation
{
	double m_Numerator{ 0.0 };
	
	double m_DenominatorConstant{ 0.0 };
	double m_DenominatorFuncMultiple{ 0.0 };
	TrigFunctions m_Function = TrigFunctions::NONE;
	PlusOrMinusDenom m_PoM;

	// Variables found after evaluation
	double m_Eccentricity{ 0.0 };
	std::string m_DirectrixString;

	explicit PolarEquation(
		double Numerator,
		double DenomConst,
		double DenomFuncMult,
		TrigFunctions FuncType,
		PlusOrMinusDenom PlusOrMinus = PlusOrMinusDenom::Plus) :
		m_Numerator(Numerator),
		m_DenominatorConstant(DenomConst),
		m_DenominatorFuncMultiple(DenomFuncMult),
		m_Function(FuncType),
		m_PoM(PlusOrMinus) {}

	//PolarEquation(const PolarEquation&) = delete;
	//PolarEquation(PolarEquation&&) = delete;
};

// where B=0 (not rotation for this function),
// and A and C  are nonzero real numbers.
// If B = 0, the conic section will have a vertical and/or horizontal axes
ConicType IdentifyConicType(
	const double& A,
	const double& C,
	const double& D,
	const double& E,
	const double& F,
	const double& B = 0);

void PrintConicType(const ConicType& InType);


ConicType IdentifyConicType(
	const double& A,
	const double& C,
	const double& D,
	const double& E,
	const double& F,
	const double& RotationAngleRadians,
	const double& B = 0);

void WriteEquationOfRatatedConicStandardForm(
	const double& A,
	const double& B,
	const double& C,
	const double& Equals);

//std::pair<double, double> OutputDecimalAsFract(double input);

ConicType UseDescrimentToFindConicType(
	const double& B,
	const double& A,
	const double& C);

ConicType GetTypeFromEccentricity(const double& Eccentricity);

PolarEquation EvaluatePolarEquation(PolarEquation& InEquation);

void PrintPolarEquationInfo(
	const std::string& InDirectrixStr,
	const double& InEccentricity,
	const ConicType& InType);

std::map<std::string, double> GetPolarGraphInfo(PolarEquation& InEquation);

double RunPolarFunction(const PolarEquation& InFunction, const double& AngleInput);

// TODO: remake this equation to deal with eccentricity inputs that 
// are not whole numbers (fractions)
// may need to do the solution involving the convert to fraction function 
// and diving the two parts
// http://mathforum.org/dr.math/gifs/todd5.6.03.gif
PolarEquation FindPolarFormGivenEccentricityAndDirectrix(
	const double& Eccentricity,
	const std::string& Directrix);

// TODO: make an equation that converts a given polar equation to rectangular form
// Use the two identities.. r = std::sqrt(x^2 + y^2) 
// and x rcos(theta), and y = rsin(theta)