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
		
		QuadraticFunction TestQuadStrRead("(x-3)^2+0");
		Limit TakeLimitQUad(TestQuadStrRead, 2.0);


	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

