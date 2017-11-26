#include "ParametricEquation.h"

#include <iostream>
#include <string>
#include "MathConstants.h"


using std::vector;
using std::pair;
using std::function;
using std::cout;
using std::endl;

//const double my_gravityfeet(32);  // f/s^2
const double pi(3.14159);


ParametricEquation::ParametricEquation()
{

}


ParametricEquation::~ParametricEquation()
{

}


vector<pair<double, double>> ParameterizeCurveFunctionTable(
	function<double(const double&)> XFunc,
	function<double(const double&)> YFunc,
	const int& StartNum,
	const int& EndNum)
{
	// store a free function
	std::function<double(const double&)> xFunc = XFunc;
	std::function<double(const double&)> yFunc = YFunc;

	vector<pair<double, double>> OutVecPairs; // 0 - 8 default slots
	const int VecSize = EndNum - StartNum;
	int CurrentNumber = StartNum;

	for ( ; CurrentNumber <= EndNum; ++CurrentNumber)
	{
		std::pair<double, double> NextPair;
		NextPair.first = xFunc(CurrentNumber);
		NextPair.second = yFunc(CurrentNumber);

		OutVecPairs.push_back(NextPair);
	}

	return OutVecPairs;

}

void FindingParametricEquationsThatModelGivenCriteria(
	std::pair<double, double> StartLoc,
	std::pair<double, double> EndLoc,
	const int& timeT
	/* const std::string MeasuredInString */)
{
	constexpr double StartTime = 0.0;
	const double EndingTime = timeT;

	const double StartX = StartLoc.first;
	const double StartY = StartLoc.second;

	const double EndX = EndLoc.first;
	const double EndY = EndLoc.second;

	double TotalXDistTraveled = EndX - StartX;
	double XSlope = TotalXDistTraveled / timeT; // measured in meters/second

	auto XFunc = [=](const double& timeT)
	{
		return (XSlope * timeT) + StartX; // y = mx + b
	};

	double TotalYDistTraveled = EndY - StartY;
	double YSlope = TotalYDistTraveled / timeT; // measured in meters/second

	auto YFunc = [=](const double& timeT)
	{
		return (YSlope * timeT) + StartY; // y = mx + b
	};

	std::vector<std::pair<double, double >> CurveTable = ParameterizeCurveFunctionTable(
		XFunc,
		YFunc,
		StartTime,
		timeT);

	PrintParametricCurvePairTable(CurveTable, StartTime, EndingTime);

}


void PrintParametricCurvePairTable(const std::vector<std::pair<double, double>>& data,
	const int& start,
	const int& end)
{
	cout << "Printing data pairs from time " << start << " to " << end << endl;

	for (const auto& datapair : data)
	{
		cout << "(" << datapair.first << "," << datapair.second << ")" << endl;
	}
}

void EliminatingtheParameterfromPairofTrigonometricParametricEquations(const double & a, const double & b)
{
	// x / a = cost
	// y / b = sint

	// cos^2(t) + sin^2(t) = 1

	std::cout << "(x/" << a << ")^2 + " << "(y/" << b << ")^2 == 1 " << std::endl;

}

void ProjectileMotionProblemSolver(
	const double & InitialSpeed,
	const double & AngleDegrees,
	const double & HeightPropelled, 
	const double & TimeSincePropelled, 
	const double & Gravity)
{
	// first translate my inputed angle into radians
	const double AngleInRadians = (pi * AngleDegrees) / 180.0;

	// The horizontal position is found using the parametric equation for x
	const double HorizontalPosX = (InitialSpeed * (std::cos(AngleInRadians)) * TimeSincePropelled);


	const double GravityG = my_gravityfeet;
	const double TimeSquared = std::pow(TimeSincePropelled, 2);
	const double NegOneHalf = -0.5;
	const double VerticalPosY =
		(NegOneHalf*GravityG*TimeSquared) +
		((InitialSpeed * std::sin(AngleInRadians)) * 
			TimeSincePropelled) + HeightPropelled;

	// TODO: After 2 seconds, the ball is 198 feet away from the batter’s box and 137 feet above the ground.
	// output string for that data
	std::cout << "After " << TimeSincePropelled << " seconds " <<
		" the ball is " << HorizontalPosX << " feet away from the " <<
		"batters box, and " << VerticalPosY << " feet above the ground"
		<< std::endl;




	// TODO: Solve the quadratic assign variables etc...
	// find out how long in air
	// set y = 0 solve 
	const double a = NegOneHalf*GravityG;
	const double b = (InitialSpeed * std::sin(AngleInRadians));
	const double c = HeightPropelled;
	
	const double TimeUntilBallHitsGround =
		UseQuadraticFormulaToFindHowLongInAir(a, b, c);

	std::cout << std::endl;

	std::cout << "The ball will hit the ground after "
		<< TimeUntilBallHitsGround << " seconds. " << std::endl;

	const double WallDistanceFromHome = 400.0;

	double HeightWhenClearedFence = 0.0;
	const bool bIsHomeRun = FindOutIfHitIsHomerun(
		WallDistanceFromHome,
		InitialSpeed,
		AngleInRadians,
		a,
		b,
		c,
		HeightWhenClearedFence);


	if (bIsHomeRun)
	{
		std::cout << "It was a homerun " << std::endl;
		std::cout << "The ball is " << HeightWhenClearedFence << " feet in the air "
			<< "when it soars out of the ballpark." << std::endl;
	}
	else
	{
		std::cout << "It was not a homerun" << std::endl;
	}

}

// returns only the positive result
double UseQuadraticFormulaToFindHowLongInAir(
	const double& a,
	const double& b,
	const double&c)
{
	const double NegativeB = b * (-1);
	const double BSquared = std::pow(b, 2);
	const double TwoAC = 2 * a * c;
	const double TwoA = 2 * a;

	
	double res1 = (NegativeB + (std::sqrt(BSquared - TwoAC))) / TwoA;
	double res2 = (NegativeB - (std::sqrt(BSquared - TwoAC))) / TwoA;


	// returns only the positive result
	if (res1 > 0)
	{
		return res1;
	}
	else
	{
		return res2;
	}



}

bool FindOutIfHitIsHomerun(const double & DistanceFromHomePlate,
	const double & initalspeed,
	const double & AngleInRadians,
	const double& a, // variables we used from the quadratic formula
	const double& b, 
	const double& c,
	double& OutHeightWhenCleared)

{
	// when the wall is 400 feet away and 10 feet high

	// set x = 400 and solve for t
	const double t = DistanceFromHomePlate / (initalspeed * (std::cos(AngleInRadians)));

	// input t into y
	const double y = (a * (std::pow(t, 2)) + (b * t) + c);

	if (y > 10.0)
	{
		OutHeightWhenCleared = y;
		return true;
	}
	else
	{
		OutHeightWhenCleared = 0;
		return false;
	}
}
