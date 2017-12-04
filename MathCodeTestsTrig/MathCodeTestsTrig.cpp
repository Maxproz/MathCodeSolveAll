// MathCodeTestsTrig.cpp : Defines the entry point for the console application.


#include "stdafx.h"

#include "CalculusFunction.h"
#include "ConstantFunction.h"
#include "TrigonometricFunction.h"
#include "Derivative.h"



// TODO: Should I create a class for typedefs? (they are inside MathConstants for now)
// TODO: Create the 5th/6th/7th degree polynomial classses
// TODO: Create the base functionality for the quotient rule for derivatives
// TODO: Setup power functions to take the derivative of a negative number exponent,
// base functionality is maybe there? Needs confirmed.
// TODO: I removed the std::unique_ptr from the rational function class to simplify copy/moving.
// Need to remember that I have that variable doing nothing in some limit functions etc...
// TODO: Create a function that uses the abs(velocity) = speed function to better understand it


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


		//RootFunction<2> RootFuncTest(1, 0, 0);
		//Derivative <RootFunction<2>, RootFunction<-2>> RootDeriv(RootFuncTest);
		//RootFunction<-2> OutFunc = RootDeriv.GetDerivativeFunction();
		//OutFunc.PrintFunction();
			
		//QuadraticFunction TestQuadratic(2, -3, 1);
		//ConstantFunction TestConstFunc =
		//	GetSecondDerivativeFunction<QuadraticFunction, LinearFunction, ConstantFunction>(TestQuadratic);

		//TestConstFunc.PrintConstantFunctionInfo();

		//PowerFunction<4> TestPowerFunc(3, 4, 5, 6);
		//Derivative<PowerFunction<4>, PowerFunction<3>> TestPowerDerivative(TestPowerFunc);
		//
		//auto DerivativeFunc = TestPowerDerivative.GetDerivativeFunction();
		//DerivativeFunc.PrintFunction();
			
		//PowerFunction<7> TestPowerFunc(1, 1, 0, 0);
		//Derivative<PowerFunction<7>, PowerFunction<6>> TestPowerDerivative(TestPowerFunc);

		//auto DerivativeFunc = TestPowerDerivative.GetDerivativeFunction();
		//DerivativeFunc.PrintFunction();
		//
		//PowerFunction<0> TestFunc1(0, 0, 0, 7);
		//PowerFunction<5> TestFunc2(2, 1, 0, 0);

		////Derivative<PowerFunction<5>, PowerFunction<4>> FirstDeriv(TestFunc2);
		////Derivative<PowerFunction<0>, PowerFunction<0>> SecondDeriv(TestFunc1);
		////
		////auto First = FirstDeriv.GetDerivativeFunction();
		////auto Second = SecondDeriv.GetDerivativeFunction();

		////First.PrintFunction();
		////Second.PrintFunction();


		//auto GetRes1 = std::get<4>(TestFunc1.GetNAKDC());
		//GetRes1 = ApplyDerivativeConstantRule(GetRes1);
		//
		//auto a = std::get<1>(TestFunc2.GetNAKDC());
		//auto n = std::get<0>(TestFunc2.GetNAKDC());

		//double newa;
		//double newn;
		//ApplyDerivativePowerRules(a, n, newa, newn);

		//std::cout << newa << "x^" << newn << "+" << GetRes1;

		//cout << endl;

		//QuadraticFunction InputFunction1(1, 0, 2);
		//CubicFunction InputFunction2(3, 0, -5, 0);

		//QuarticFunction TestFunc = ApplyDerivativeProductRule(InputFunction1, InputFunction2);

		//TestFunc.PrintFunction();
		//cout << endl;
		

		//QuadraticFunction TestNumer(5, 0, 0);
		//LinearFunction TestDenom(4, 3);

		//RationalFunction<QuadraticFunction, QuadraticFunction> TestRatRes =	ApplyDerivativeQuotientRule(TestNumer, TestDenom);
		//TestRatRes.PrintFunction();
		//cout << endl;

		//LinearFunction TestNumer(3, 1);
		//LinearFunction TestDenom(4, -3);

		//RationalFunction<LinearFunction, QuadraticFunction> TestRatRes = ApplyDerivativeQuotientRule(TestNumer, TestDenom);
		//TestRatRes.PrintFunction();
		//cout << endl;

		//PowerFunction<-7> TestPowerFunc(1, 1, 0, 0);
		//Derivative<PowerFunction<-7>, PowerFunction<-8>> Deriv(TestPowerFunc);

		//auto DerivativeFunc = Deriv.GetDerivativeFunction();

		//DerivativeFunc.PrintFunction();
		//cout << endl;

		//PowerFunction<-2> TestPowerFunc(6, 1, 0, 0);
		//Derivative<PowerFunction<-2>, PowerFunction<-3>> Deriv(TestPowerFunc);

		//auto DerivativeFunc = Deriv.GetDerivativeFunction();

		//DerivativeFunc.PrintFunction();
		//cout << endl;

		//CubicFunction TestCubic(1, -7, 8, 1);
		//QuadraticFunction DerivFunc = TestCubic.GetDerivativeFunction();
		////DerivFunc.PrintFunction();
		////cout << endl;

		//TestCubic.PrintHorizontalTangetLineXValues();
		//cout << endl;
		//C:\Users\Maxpro\Downloads

		
		//PrintAllFilesInDirectory("C:/Users/Maxpro/Downloads");

		//double VelocityAtOne{ 0 };
		//double AccelerationAtOne{ 0 };

		//CubicFunction PositionFunction(1, 0, -4, 2);
		//
		//QuadraticFunction VelocityFunction = PositionFunction.GetDerivativeFunction();
		//VelocityAtOne = VelocityFunction(1);

		//LinearFunction AccelerationFunction = VelocityFunction.GetDerivativeFunction();
		//AccelerationAtOne = AccelerationFunction(1);

		//std::cout << "Velocity at time t = 1: is " << VelocityAtOne << endl;
		//cout << "Acceleration at time t = 1: is " << AccelerationAtOne << endl;


		// f(t) =
		// Where t >= 0
		//const int TimeTEqualZero = 0;


		//CubicFunction PositionFunction(1, -9, 24, 4);
		//DisplayParticlePositionAndVelocityMovingAlongAxis<CubicFunction>(PositionFunction);
		
		//const int TimeT = 3;
		//QuadraticFunction PositionQuad(1, -5, 1);
		//FindParticleMovementAtTimeT(PositionQuad, TimeT);

		
		//EstimatePopulationSizeAfterTDays(3000, 100, 3);

		//LinearFunction PriceFunction(-0.03, 9);
		//const int StartInterval = 0;
		//const int EndInterval = 300;
		//const int Items = 101;
		//EstimateTheRevenueObtainedFromSellingXItems(StartInterval, EndInterval, Items, PriceFunction);

		//QuadraticFunction ProfitFunction(-0.03, 8, -50);
		//UseMarginalProfitToEstimateSaleOfXthItem(101, ProfitFunction);
		//cout << endl;
		
		//MPSIN<1> SinTestFunc1(1, 1, 0, 0);
		//cout << SinTestFunc1(2) << endl;

		// Only <MPSIN, 1> stores the derivative at the moment
		//TrigometricFunction<MPSIN, 2> SinTestFunc2(1, 1, 0, 0);

		
		//cout << SinTestFunc2(2) << endl;

		// TODO: Add the Derivative function setup to these classes so I can grab it like the Quadratic Function class
		//MPCOS<1> CosDerivTestFunc1 = SinTestFunc1.GetDerivativeFunction();
		//cout << CosDerivTestFunc1(2) << endl;
		//MPNEGSIN<1> TestFunc = CosDerivTestFunc1.GetDerivativeFunction();
		//cout << TestFunc(2) << endl;


		//SinTestFunc1(2)
		//cout << endl;
		//cout << CosDerivTestFunc1(2) << endl;

		//TrigometricFunction<MPSIN, 2> SinTestFunc2(1, 1, 0, 0);


		//MPCOT<1> COTTESTFUNC(1, 1, 0, 0);
		//cout << COTTESTFUNC(2) << endl;
		//
		//MPNEGCSC<2> NEGCSCSQUAREDDERIVFUNC = COTTESTFUNC.GetDerivativeFunction();
		//cout << NEGCSCSQUAREDDERIVFUNC(2) << endl;

		//MPSEC<1> SECTESTFUNC(1, 1, 0, 0);
		//cout << SECTESTFUNC(2) << endl;

		//MPSECTAN<1> SECXTANXDERIVFUNC = SECTESTFUNC.GetDerivativeFunction();
		//cout << SECXTANXDERIVFUNC(2) << endl;

		//MPCSC<1> CSCTESTFUNC(1, 1, 0, 0);
		//cout << CSCTESTFUNC(2) << endl;

		//MPNEGCSCCOT<1> NEGCSCXCOTXDERIVFUNC = CSCTESTFUNC.GetDerivativeFunction();
		//cout << NEGCSCXCOTXDERIVFUNC(2) << endl;

		
		// TODO: Test out the "Start New Feature" button in sourcetree

		MPCOT<1> TestFunc(1, 1, 0, 0);
		std::string TestStr = TestFunc.GetEquationForATangentLineAtInputAngle(M_PIOverFour);
		cout << TestStr << endl << endl;


	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

