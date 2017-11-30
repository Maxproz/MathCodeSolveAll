#pragma once



#include <functional>
#include <algorithm>    // std::find_if

/* old includes before refactor */
#include "Circle.h"
/* newest includes */

#include "TranscendentalFunction.h"
#include "PolynomialFunction.h"
#include "CubicFunction.h"
#include "QuadraticFunction.h"
#include "LinearFunction.h"
#include "RootFunction.h"
#include "MathConstants.h" // for euluere
#include "RationalFunction.h"
#include "PiecewiseFunction.h"
#include "MiscMathEquations.h"
#include "Limit.h"
#include "Derivative.h"

//typedef pair<double, double> Point;


// TODO: Move the important constants into one header
// TODO: Calculate the domain and ranges of each individual function after splitting up files (and other data)
// TODO: Set the domain and range variables in the constructors of all the different function types.

// Basic Important Notes on Functions
/* 1. (algebraic functions can only have powers that are rational numbers.)
*/


void PrintPointSlopeForm(const double& Slope, const Point& InPoint);

double GetSlope(const Point& FirstPoint, const Point& SecondPoint);

void PrintSlopeInterceptForm(const Point& Point, const double& Slope);


// Find the rate of change of profit when 10,000 games are produced.
/*
The profit P(x) earned by producing xx gaming systems is R(x)−C(x), 
where R(x) is the revenue obtained from the sale of x games. 
Since the company can sell x games at p=−0.01x+400 per game,

R(x) = x*(−0.01x+400)


@Input 1: Quadratic Reveune Function
@input 2: Linear Cost Function

@Output: Bool for weather you should increase production based on profit calculations with derivative
*/

inline bool ShouldProductionBeIncreasedUsingRateOfChange(const QuadraticFunction& Revenue, const LinearFunction& Cost, const int& AmountProduced)
{
	// Calculate a new Profit Quadratic Function P(x) = R(x) - C(x) Inside
	// Take the derivative of that new function for 10000 games produced

	double NewB{ 0 };
	double NewC{ 0 };

	// A stays the same when subtracting a linear from a quadratic (what this function is specific for)
	double NewA = std::get<0>(Revenue.GetABC());
	NewB = std::get<1>(Revenue.GetABC()) - (Cost.GetA());
	NewC = std::get<2>(Revenue.GetABC()) - (Cost.GetB());

	QuadraticFunction Profit(NewA, NewB, NewC);
	double ProfitOfAmountProduced = Profit(AmountProduced);

	Derivative<QuadraticFunction, LinearFunction> DerivativeOfProfit(Profit);
	double RateOfChangeOfAmountProduced = DerivativeOfProfit.EstimateDerivative(AmountProduced);


	// Make an if else for the sentence below
	// Since the rate of change of profit P′(10,000) > 0 AND P(10,000) > 0, the company should increase production.
	// TODO: When should decrease production?
	
	if (ProfitOfAmountProduced > 0 && RateOfChangeOfAmountProduced > 0)
	{
		// increase production
		return true;
	}
	else
	{
		// decrease production
		// NOTE: cant help but think i am missing more cases here
		// what about Profit > 0, Rate of change < 0 ?
		// I will probably end up learning later once i learn more about statistics.
		return false;
	}
}

inline bool ShouldProductionBeIncreasedUsingRateOfChange(const QuadraticFunction& Profit, const double& PriceOfItem)
{

	double ProfitOfAmountProduced = Profit(PriceOfItem);

	QuadraticFunction LocalProfitFunc = Profit;

	Derivative<QuadraticFunction, LinearFunction> DerivativeOfProfit(LocalProfitFunc);
	double RateOfChangeOfAmountProduced = DerivativeOfProfit.EstimateDerivative(PriceOfItem);


	// Make an if else for the sentence below
	// Since the rate of change of profit P′(10,000) > 0 AND P(10,000) > 0, the company should increase production.
	// TODO: When should decrease production?

	if (ProfitOfAmountProduced > 0 && RateOfChangeOfAmountProduced > 0)
	{
		// increase production
		return true;
	}
	else
	{
		// decrease production
		// NOTE: cant help but think i am missing more cases here
		// what about Profit > 0, Rate of change < 0 ?
		// I will probably end up learning later once i learn more about statistics.
		return false;
	}
}







// Start of calculus chapter 
// TODO: Sort/organize this stuff later

// Slope of a secant line formula
// through Points (a,f(a)) and (x,f(x))
inline double GetSlopeOfSecantLineTwoPoints(const Point& FirstPoint, const Point& SecondPoint)
{
	double a = FirstPoint.first;
	double FOfa = FirstPoint.second;

	double x = SecondPoint.first;
	double FOfx = SecondPoint.second;

	// use formula return result
	double Numerator = FOfx - FOfa;
	double Denominator = x - a;

	double SlopeOfSecantLine = Numerator / Denominator;

	return SlopeOfSecantLine;
}

//template <typename Func>
//inline double GetSlopeOfSecantLine(const double& x, const double& a, const Func& F)
//{
//	// use formula return result
//	double Numerator = F. - F(a);
//	double Denominator = x - a;
//
//	double SlopeOfSecantLine = Numerator / Denominator;
//
//	return SlopeOfSecantLine;
//}

//bool IsOdd(int i) 
//{
//	return ((i % 2) == 1);
//}



LinearFunction FactorNumeratorRationalFunc(const QuadraticFunction& InNumerator, const LinearFunction& InDenominator);

inline double GetSlopeOfSecantLine(const double& a, const QuadraticFunction& F)
{
	// use formula return result

	auto TupleABC = F.GetABC();
	double QuadA = std::get<0>(TupleABC);
	double QuadB = std::get<1>(TupleABC);
	double QuadC = std::get<2>(TupleABC);

	double FofA = F(a);
	std::cout << "FofA = 9 == " << FofA << endl;

	QuadC = QuadC - FofA;

	QuadraticFunction Numerator = QuadraticFunction(QuadA, QuadB, QuadC);
	//double InA = a*-1;
	LinearFunction Denominator(1, a * (-1));

	LinearFunction FactoredFunc = FactorNumeratorRationalFunc(Numerator, Denominator);

	FactoredFunc.PrintLinearFunctionInfo();

	// Take the limit of the factored function to find the slope of the secant line
	Limit<LinearFunction> LocalLimit(FactoredFunc, a);
	double SlopeOfSecantLine = LocalLimit.GetLimitResult();

	//// Next, find a point on the tangent line.

	//// Since the line is tangent to the graph of f(x)x=3,
	//// it passes through the point (3,f(3)).
	//// We have f(3)=9, so the tangent line passes through the point (3,9).
	//
	//// tangent line at x = a
	//LinearFunction TangetLineAtXIsEqualToA(1, 0);
	//ConvertFromPointSlopeFromToLinearForm(Point(a, FofA), SlopeOfSecantLine, TangetLineAtXIsEqualToA);



	// TODO: Finish next steps

	return SlopeOfSecantLine;
}

// Assumes the zeros are not fractions
inline LinearFunction FactorNumeratorRationalFunc(const QuadraticFunction& InNumerator, const LinearFunction& InDenominator)
{
	// assumes a simple factor with - QuadraticFunction m_a == 1;
	auto NumABTemp = InNumerator.GetAB();
	double ExceptionCheckA = std::get<0>(NumABTemp);
	if (ExceptionCheckA != 1)
	{
		throw std::logic_error("Simple factor is not able to handle this quadratic");
	}


	auto QuadNumVec = InNumerator.GetRealNumberZerosVec();
	auto LinearNumVec = InDenominator.GetRealNumberZerosVec();

	cout << QuadNumVec.size() << endl;
	cout << LinearNumVec.size() << endl;
	
	//for (int i = 0; i < QuadNumVec.size(); ++i)
	//{
	//	cout << QuadNumVec[i] << endl;
	//}

	// Have both zeros vectors 
	// assign a variable to the denominator zero
	// check to see if this variable is in the numerator
	// if it is "erase" it from numerator and denominator
	// and return a linear function that displays the factored form

	//int i = 0;

	if (LinearNumVec.size() >= 1)
	{
		double NumeratorZero = LinearNumVec.front();

		auto TryFindZero = std::find(QuadNumVec.begin(), QuadNumVec.end(), NumeratorZero);
		
		if (TryFindZero != QuadNumVec.end())
		{
			//std::cout << "DEBUG 1" << endl;
			
			// Found a zero
			// "erase" it by reformatting a quadratic / a linear
			double MatchingZero = *TryFindZero;
			std::cout << "MatchingZero: " << *TryFindZero << endl;

			QuadNumVec.erase(TryFindZero);
			QuadNumVec.shrink_to_fit();
			LinearNumVec.pop_back();

			//std::sort(QuadNumVec.begin(), QuadNumVec.end(), std::greater<double>());
		
			//LinearNumVec.pop_back();

			//for (int i = 0; i < QuadNumVec.size(); ++i)
			//{
			//	cout << QuadNumVec[i] << endl;
			//}

			//auto NumRemove = std::remove(QuadNumVec.begin(), QuadNumVec.end(), MatchingZero);
			//auto DenomRemove = std::remove(LinearNumVec.begin(), LinearNumVec.end(), MatchingZero);
		}
		else
		{
			throw std::logic_error("Zero found in denominator not found in numerator");
		}
	}
	else
	{
		throw std::logic_error("Assuming the wrong func form FactorNumeratorRationalFunc");
	}


	// here is where I would add checks for fractional numbers


	// Set up the remaining linear function to be returned after factoring
	double NewAVar{ 0.0 };
	double NewBVar{ 0.0 };

	// assumes a simple factor with - QuadraticFunction m_a == 1;
	NewAVar = 1;

	cout << QuadNumVec.size() << endl;
	cout << LinearNumVec.size() << endl;

	double OtherBZero = QuadNumVec.front();
	// current form x = OtherBZero


	cout << "LastZero: " << OtherBZero << endl;

	//double LHSBForm{ 0.0 };

	// Flip the sign so the function we return is in the correct form.
	FlipSign<double>(OtherBZero);

	//if (OtherBZero < 0)
	//{
	//	LHSBForm = LHSBForm + (OtherBZero * (-1));
	//}
	//else if (OtherBZero > 0)
	//{
	//	LHSBForm = LHSBForm + (OtherBZero * (-1));
	//}
	//else
	//{
	//	// == 0 ?? I will see the reaction to this later
	//	throw std::logic_error("Undefined Logic Error");
	//}
	
	LinearFunction OutFunc(NewAVar, OtherBZero);

	return OutFunc;

}



inline void FindEquationOfTangentLineToAFunction(const double& a, const QuadraticFunction& F)
{
	// use formula return result

	auto TupleABC = F.GetABC();
	double QuadA = std::get<0>(TupleABC);
	double QuadB = std::get<1>(TupleABC);
	double QuadC = std::get<2>(TupleABC);

	double FofA = F(a);
	std::cout << "FofA = 9 == " << FofA << endl;

	QuadC = QuadC - FofA;

	QuadraticFunction Numerator = QuadraticFunction(QuadA, QuadB, QuadC);
	//double InA = a*-1;
	LinearFunction Denominator(1, a * (-1));

	LinearFunction FactoredFunc = FactorNumeratorRationalFunc(Numerator, Denominator);

	FactoredFunc.PrintLinearFunctionInfo();

	// Take the limit of the factored function to find the slope of the secant line
	Limit<LinearFunction> LocalLimit(FactoredFunc, a);
	double SlopeOfSecantLine = LocalLimit.GetLimitResult();

	// Next, find a point on the tangent line.

	// Since the line is tangent to the graph of f(x)x=3,
	// it passes through the point (3,f(3)).
	// We have f(3)=9, so the tangent line passes through the point (3,9).

	// tangent line at x = a
	LinearFunction TangentLineAtXIsEqualToA(1, 0);
	ConvertFromPointSlopeFromToLinearForm(Point(a, FofA), SlopeOfSecantLine, TangentLineAtXIsEqualToA);

	std::cout << "Original Function: " << endl;
	F.PrintFunction();
	cout << endl;
	cout << "Equation of line tangent to origional function at x == " << a << endl;
	TangentLineAtXIsEqualToA.PrintLinearFunctionInfo();

	return;

	//return SlopeOfSecantLine;
}

//  if h≠0 is chosen so that a+h is in the interval
template <typename Func>
inline double GetSlopeOfSecantLineTwoPointsWithIncrementH(const double& a, const double& h, const Func& F)
{
	//double a = FirstPoint.first;
	//double FOfa = FirstPoint.second;

	//double a = SecondPoint.first;
	//double FOfx = SecondPoint.second;

	//// use formula return result
	//double Numerator = FOfx - FOfa;
	//double Denominator = x - a;

	//double SlopeOfSecantLine = Numerator / Denominator;

	double SlopeOfSecantLine{ 0 };

	double Numerator = (F(a + h)) - (F(a));
	double Denominator = h;

	SlopeOfSecantLine = Numerator / Denominator;
	//std::cout << Numerator << std::endl;
	//std::cout << Denominator << std::endl;

	return SlopeOfSecantLine;
}


// TODO: This functions variables need to be better understood
//  average velocity of an object over a time period to be the change in its position divided by the length of the time period.
inline double GetAverageVelocity(const double& a, const double& t,
	std::function<double(const double&)> s)
{
	// x == t
	/*
	Let s(t)s(t) be the position of an object moving along a coordinate axis at time t. 
	The average velocity of the object over a time interval [a,t] where a<t
	(or [t,a] if t<a)t<a) 
	is
	
	s(t) - s(a)
	/
	t - a

	*/
	if (a < t)
	{
		std::cout << "a < t using iterval [a, t]\n";

		// use formula return result
		double Numerator = s(t) - s(a);
		double Denominator = t - a;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}

	if (a > t)
	{

		std::cout << "t < a using iterval [t, a]\n";
		// use formula return result
		double Numerator = s(a) - s(t);
		double Denominator = a - t;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}
}

// Interval [0,3] under estimate means
// evaluate the function for 0,1,2 and add the results
// underestimate = false means evaluate for 1 2 3
inline double GetAreaUnderCurve(const int& IntervalStart, const int& IntervalEnd,
	std::function<double(const double&)> func,
	const bool& bIsUnderEstimate = true)
{
	double OutResult{ 0.0 };

	if (bIsUnderEstimate)
	{
		for (int i = IntervalStart; i < IntervalEnd; i++)
		{
			OutResult = OutResult + func(i);
		}
	}

	// TODO: is this bottom part correct?
	if (bIsUnderEstimate == false)
	{
		for (int i = IntervalStart + 1; i <= IntervalEnd; ++i)
		{
			OutResult = OutResult + func(i);
		}
	}

	return OutResult;
}

//template <typename FirstFunc, typename SecondFunc>
//bool DetermineContinunityAtAPoint(RationalFunction<FirstFunc,SecondFunc>& InFunc, const int& InPoint);

bool DetermineContinunityAtAPoint(PiecewiseFunction<QuadraticFunction, LinearFunction>& InFunc, const int& InPoint);


// TODO: How I did this is super sloppy and specific, should alter it somehow.
template <typename FirstFunc, typename SecondFunc, int ThirdFunctionConstant>
inline bool DetermineContinunityAtAPoint(const PiecewiseFunctionThreeFunctions<FirstFunc, SecondFunc, ThirdFunctionConstant>& InFunc, const int& InPoint)
{
	//QuadraticFunction FirstFuntion = InFunc.GetFirstFunction();
	//LinearFunction SecondFunction = InFunc.GetSecondFunction();

	// Step 1: check to see if f(a) is defined
	double TestOne = InFunc(InPoint);

	if (std::isnan(TestOne))
	{
		// failed 
		std::cout << "TestOneFailed: The function is not continuous at " << InPoint << "\n";
		return false;
	}
	else
	{
		//  If f(a) is defined, continue to step 2.

		// Step 2: Compute Limit from both sides
		// If Limit does not exist (that is, it is not a real number),
		// then the function is not continuous at a and the problem is solved.



		Limit TestTwoLimit(InFunc, InPoint);
		double TestTwo = TestTwoLimit.GetLimitResult();

		// if Limit Exists go to step 3
		// TODO: Need a rational function check to return bool to see if it exists
		// TODO: Need more example data in order to fix

		if (TestOne != TestTwo)
		{
			std::cout << "TestThreeFailed: The function is not continuous at " << InPoint << "\n";
			std::cout << "Because f(" << InPoint << ") " << " = " << TestOne << " != " << TestTwo << " = Limit fof(x)_x->m_a \n";
			return false;
		}
		else
		{
			std::cout << "TestThreePassed: The function is continuous at " << InPoint << "\n";
			return true;
		}

	}
}


/* Types of discontinunities 
Intuitively, a removable discontinuity is a discontinuity for which there is a hole in the graph,
a jump discontinuity is a noninfinite discontinuity for which the sections of the function do not meet up,
and an infinite discontinuity is a discontinuity located at a vertical asymptote.
*/

// If a function is not continuous then you can classify what type its discontinunity is.

// 1. Removeable if the limit exists

// 2. Jump if limit exists from both sides but they are not equal.

// 3. Infinity if limit is +- infinity on both sides

inline bool ApplyIntermediateValueTherom(const CubicFunction& InCubicFunc, const int& ClosedIntervalStart, const int& ClosedIntervalEnd)
{
	auto FullFunctionForm = InCubicFunc.GetABCD();

	//auto a = std::get<0>(FullFunctionForm);
	//auto b = std::get<1>(FullFunctionForm);
	//auto c = std::get<2>(FullFunctionForm);
	//auto d = std::get<3>(FullFunctionForm);

	double StartIntervalRes = InCubicFunc(ClosedIntervalStart);
	double EndIntervalRes = InCubicFunc(ClosedIntervalEnd);

	bool OppositeSigns = (StartIntervalRes > 0 && EndIntervalRes < 0) || (StartIntervalRes < 0 && EndIntervalRes > 0);

	if (OppositeSigns)
	{
		// There exists at least one zero in the interval
		return true;
	}
}


//inline void ProveLinearFunctionLimitEpsilonDelta(const Limit& LinearFunctionLimit)
//{
//	LinearFunction LinearFunc = LinearFunctionLimit.GetLinearFunctionIfExists();
//
//	double a = LinearFunc.GetA();
//	double b = LinearFunc.GetB();
//
//	// Let  epsilon > 0
//
//	double Epsilon = 1;
//	double Delta = 1;
//	// This means we must prove that whatever follows is true no matter what positive value of ε is chosen.
//
//	auto L = LinearFunctionLimit.GetLimitResult();
//
//	// Factor the form
//	// (ax + b) = L
//	if (L > 0)
//	{
//		b = b - L;
//	}
//	else if (L < 0)
//	{
//		b = b + L;
//	}
//	else
//	{
//		// shouldnt reach here
//		throw std::logic_error("reached invalid location Prove epsilon delta linear func");
//	}
//
//	// current form is |ax + b| < epsilon
//	// check if needs factoring
//	
//	// Should only need this variable sometimes?
//	double Divisor = 0;
//	
//	if (a == b)
//	{
//		Divisor = a;
//
//		a = a / Divisor;
//		b = b / Divisor;
//		// Epsilion = Epsilion / Divisior;
//
//		// current form
//		// std::abs(LinearFunction NewForm(a, b)) <  Epsilion / Divisior;
//		
//		// Thus, it would seem that delta = Epsilion / Divisior is appropriate.
//		// Delta is std::abs(LinearFunction NewForm(a, b)) 
//		// Delta = Epsilion / Divisior;
//
//
//	}
//
//	// Now Assume 0 < | x − 1 | < Delta
//	// When Delta has been chosen
//
//	// Then 0<|x−1|< delta ,  then |(2x+1)−3|< epsilion.
//
//	// |2| |x-1|
//	// 2 *|x-1|
//	// < 2 * Delta ~~~~~~ here’s where we use the assumption that 0<|x−1|< delta
//	// Let our choice of delta = epsilion / divisior
//	// = 2*epsilion/divisor  // if divisior == 2  = epsilion
//
//
//}
