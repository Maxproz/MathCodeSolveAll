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

const double my_gravityfeet(32);  // f/s^2


long mygcd(long a, long b);


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

double TestFunc(const double& x)
{
	return std::sqrt(x - 3.0);
}

int main()
{
	try
	{
		//auto func = TestFunc;

		//Limit TestLimit(TestFunc, 3.0);
		
		//QuadraticFunction TestQuadStrRead("(x-3)^2+0");
		//Limit TakeLimitQUad(TestQuadStrRead, 2.0);

		//std::pair<double, LinearFunction> NumeratorA(1, LinearFunction(1,1));
		//std::pair<int, int> NumeratorConstantFraction(-1,2);

		//LinearFunction Denominator(1,-1);

		//ComplexFraction TestComplexFract(NumeratorA, NumeratorConstantFraction,
		//	Denominator);

		//Limit TestComplexFractLimit(TestComplexFract, 1);

		//RootFunction TestRootFunc(2, 1, 3, 0);
		//Limit TestRootLimit(TestRootFunc, 3);

		//LinearFunction FirstFunc(4, -3);
		//QuadraticFunction SecondFunc("(x-3)^2 + 0");

		//PiecewiseFunction<LinearFunction, QuadraticFunction> 
		//	PiecewiseFunc(FirstFunc, "<", SecondFunc, ">=", 2);

		//Limit PiecewiseLimitTest(PiecewiseFunc, 2);

		//QuadraticFunction QuadTestOne(1, 0, -4);
		//LinearFunction LinearTestOne(1, -2);

		//RationalFunction RationalTestOne(QuadTestOne, LinearTestOne);
		//
		//const int PointToCheck = 2;

		//DetermineContinunityAtAPoint(RationalTestOne, PointToCheck);



	//	RationalFunction RationalTestOne(QuadTestOne, LinearTestOne);

		


		//QuadraticFunction QuadTestOne(-1, 0, 4);
		//LinearFunction LinearTestOne(4, -8);

		//const int PointToCheck = 3;


		//PiecewiseFunction<QuadraticFunction, LinearFunction>
		//	PiecewiseFunc(QuadTestOne, "<=", LinearTestOne, ">", 3);

		//DetermineContinunityAtAPoint(PiecewiseFunc, PointToCheck);

		LinearFunction TestFuncOne(2, 1);
		const int ConstFuncReturn = 2;
		LinearFunction TestFuncTwo(-1, 4);


		PiecewiseFunctionThreeFunctions<LinearFunction, LinearFunction, ConstFuncReturn>
			PiecewiseFuncOne(TestFuncOne, "<", TestFuncTwo, ">",1, "==");

		DetermineContinunityAtAPoint(PiecewiseFuncOne, 1);


	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

