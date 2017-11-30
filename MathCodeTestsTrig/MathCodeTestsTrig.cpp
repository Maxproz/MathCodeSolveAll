// MathCodeTestsTrig.cpp : Defines the entry point for the console application.


#include "stdafx.h"

#include "CalculusFunction.h"
#include "ConstantFunction.h"
#include "TrigonometricFunction.h"
#include "Derivative.h"


// TODO: Should I create a class for typedefs?
//typedef std::pair<double, double> Point;


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

		//TrigometricFunction<MPSIN> SinFunctionTest(1, 1, 0, 0);
		//auto Result = SinFunctionTest(0.1, true);
		//std::cout << Result << std::endl;


		//QuadraticFunction TestQuad(1, -2, 0);
		//QuadraticFunction TestQuad(1, 0, 0);
		//Derivative<QuadraticFunction, LinearFunction> TestDerivative(TestQuad);

		//LinearFunction LinearFunctionTest(0, 0);
		//LinearFunctionTest = TestDerivative.GetDerivativeFunction();
		//LinearFunctionTest.PrintLinearFunctionInfo();
			

		//LinearFunction TestLinear(-3, 2);
		//ConstantFunction TestConst(0);

		//Derivative<LinearFunction, ConstantFunction> TestDerivative(TestLinear);

		//TestConst = TestDerivative.GetDerivativeFunction();
		//TestConst.PrintConstantFunctionInfo();
		
		//QuadraticFunction TestQuad(1, 0, 0);
		//FindEquationOfTangentLineToAFunction(3, TestQuad);

		//std::cout << Result << std::endl;

		//QuadraticFunction TestQuad(1, 0, 0);
		//Derivative<QuadraticFunction, LinearFunction> DerivativeTest(TestQuad);
		//auto DerivativeResult = DerivativeTest.EstimateDerivative(3);
		
		//QuadraticFunction TestQuad(3, -4, 1);
		//Derivative<QuadraticFunction, LinearFunction> DerivativeTest(TestQuad);
		//auto DerivativeResult = DerivativeTest.EstimateDerivative(2);

		//std::cout << "DerivativeResult: = " << DerivativeResult << endl;


		//QuadraticFunction TestQuad(1, 3, 2);
		//Derivative<QuadraticFunction, LinearFunction> DerivativeTest(TestQuad);
		//auto DerivativeResult = DerivativeTest.EstimateDerivative(1);

		//std::cout << "DerivativeResult: = " << DerivativeResult << endl;

		//TrigometricFunction<MPSIN> TestTrig(1, 1, 0, 0);
		//Derivative<decltype(TestTrig), TrigometricFunction<MPCOS>> TrigDerivTest(TestTrig);

		////TrigometricFunction<MPCOS> OutDeriv = TrigDerivTest.GetDerivativeFunction();
		//TrigDerivTest.EstimateDerivative(0);

		//TrigometricFunction<MPSIN> TestTrig(1, 1, 0, 0);
		//Derivative<decltype(TestTrig), TrigometricFunction<MPCOS>> TrigDerivTest(TestTrig);

		////TrigometricFunction<MPCOS> OutDeriv = TrigDerivTest.GetDerivativeFunction();
		//TrigDerivTest.EstimateDerivative(0);
		
		//QuadraticFunction TestQuad(-16, 0, 64);

		//Derivative<QuadraticFunction, LinearFunction> TestDeriv(TestQuad);
		//TestDeriv.EstimateDerivative(1);

		//// Revenue = PriceFunction * x
		//QuadraticFunction Revenue(-0.01, 400, 0);
		//LinearFunction Cost(100,10000);

		//int* AmountProduced = new int(10000);
		//bool ShouldIncreaseProduction = ShouldProductionBeIncreasedUsingRateOfChange(Revenue, Cost, *AmountProduced);

		//if (ShouldIncreaseProduction)
		//{
		//	cout << "Increase Production" << endl;
		//}
		//else
		//{
		//	cout << "Decrease Production" << endl;
		//}

		//delete AmountProduced;

		// Revenue = PriceFunction * x


		//QuadraticFunction Profit(-20, 150,-10);
		//std::unique_ptr<double> PriceOfItem = std::make_unique<double>(3.25);

		//bool ShouldIncreaseProduction = ShouldProductionBeIncreasedUsingRateOfChange(Profit, *PriceOfItem.get());

		//if (ShouldIncreaseProduction)
		//{
		//	cout << "Increase Production" << endl;
		//}
		//else
		//{
		//	cout << "Decrease Production" << endl;
		//}


		RootFunction<2> RootFuncTest(1, 0, 0);
		Derivative <RootFunction<2>, RootFunction<-2>> RootDeriv(RootFuncTest);
		RootFunction<-2> OutFunc = RootDeriv.GetDerivativeFunction();
		OutFunc.PrintFunction();




	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

