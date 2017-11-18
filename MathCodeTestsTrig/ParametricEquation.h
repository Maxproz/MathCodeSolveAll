#pragma once



#include <vector>
#include <functional>
#include <utility>



class ParametricEquation
{
public:
	ParametricEquation();
	~ParametricEquation();

};


void FindingParametricEquationsThatModelGivenCriteria(
	std::pair<double, double> StartLoc,
	std::pair<double, double> EndLoc,
	const int& timeT
	/* const std::string MeasuredInString */);

std::vector<std::pair<double, double >> ParameterizeCurveFunctionTable(
	std::function<double(const double&)> XFunc,
	std::function<double(const double&)> YFunc,
	const int& StartNum,
	const int& EndNum);

void PrintParametricCurvePairTable(const std::vector<std::pair<double, double>>& data,
	const int& start,
	const int& end);


// example: // x = acos(t) // y = bsin(t)
void EliminatingtheParameterfromPairofTrigonometricParametricEquations(
	const double& a,
	const double& b);


// TODO: make a function to evaluate trig parametric equations sin cos


void ProjectileMotionProblemSolver(
	const double& InitialSpeed,
	const double& AngleDegrees,
	const double& HeightPropelled,
	const double& TimeSincePropelled,
	const double& Gravity = 9.8);


// returns only the positive result
double UseQuadraticFormulaToFindHowLongInAir(
	const double& a,
	const double& b,
	const double&c);

bool FindOutIfHitIsHomerun(
	const double& DistanceFromHomePlate,
	const double& initalspeed,
	const double& AngleInRadians,
	const double& a, // variables we used from the quadratic formula
	const double& b,
	const double& c,
	double& OutHeightWhenCleared);
