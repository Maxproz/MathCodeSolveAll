﻿// MathCodeTestsTrig.cpp : Defines the entry point for the console application.


#include "stdafx.h"

#include "CalculusFunction.h"
#include "ConstantFunction.h"
#include "TrigonometricFunction.h"

// TODO: Should I create a class for typedefs?
typedef std::pair<double, double> Point;


int main()
{
	try
	{


		//LinearFunction Numer(1, 2);
		//LinearFunction Denom(1, 1);

		//RationalFunction RatFunc(Numer, Denom);
		//const int PointToCheck = -1;

		//RatFunc(PointToCheck);

		//if (RatFunc.GetAmountOfDiscontinunitiesFound() > 0)
		//{

		//	bool bIsContinous = DetermineContinunityAtAPoint(RatFunc, PointToCheck);
		//}

		//auto ResPtr = RatFunc.GetCurrentDiscontinunityPtrInfo();

		//std::cout << "\nPrinting First Discontinunity Info\n";
		//PrintDiscontinunityType(ResPtr.first);
		//std::cout << "At x = " << ResPtr.second << std::endl;


		//CubicFunction TestCubic(1, -1, -3, 1);

		//int* IntervalStart = new int(0);
		//int* IntervalEnd = new int(1);
	
		//bool CheckIfAZeroExists = ApplyIntermediateValueTherom(
		//	TestCubic, 
		//	*IntervalStart,
		//	*IntervalEnd);

		//if (CheckIfAZeroExists)
		//{
		//	std::cout << "We found at least one zero" << std::endl;
		//}
		// 
		//delete IntervalStart;
		//delete IntervalEnd;


		//LinearFunction TestLinear(2, 1);

		//const int x = 1;

		//Limit TestLimit(TestLinear, x);
		//
		////std::cout << Res << std::endl;

		//// prove using epsilon delta

		//ProveLinearFunctionLimitEpsilonDelta(TestLimit);


		//QuadraticFunction TestQuad(2, -3, 1);
		//LinearFunction TestLinear(5, 4);

		//RationalFunction<QuadraticFunction, LinearFunction> TestRational(TestQuad, TestLinear);
		//

		//Limit<decltype(TestRational)> TestLimit(TestRational, 3);

		//Limit<RationalFunction<QuadraticFunction, LinearFunction>> TestLimit;
		

		//QuadraticFunction TestQuad1(1, 0, -4);
		//LinearFunction TestLinear1(1, -2);

		//RationalFunction<QuadraticFunction, LinearFunction> TestRational1(TestQuad1, TestLinear1);

		//Limit<decltype(TestRational1)> TestLimit1(TestRational1, 2);

		//ConstantFunction ConstFunc1(1);
		//QuadraticFunction TestQuad1("(x-2)^2 + 0");


		//RationalFunction<ConstantFunction, QuadraticFunction> TestRational1(ConstFunc1, TestQuad1);

		//Limit<decltype(TestRational1)> TestLimit1(TestRational1, 2);

		TrigometricFunction<MPSIN> SinFunctionTest(1, 1, 0, 0);
		auto Result = SinFunctionTest(0.1, true);
		std::cout << Result << std::endl;

	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

