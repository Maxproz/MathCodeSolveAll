// MathCodeTestsTrig.cpp : Defines the entry point for the console application.




#include "stdafx.h"

#include "CalculusFunction.h"


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


		LinearFunction TestLinear(2, 1);

		const int x = 1;

		Limit TestLimit(TestLinear, x);
		
		//std::cout << Res << std::endl;

		// prove using epsilon delta

		ProveLinearFunctionLimitEpsilonDelta(TestLimit);

	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
		std::cout << std::endl;
	}


    return 0;
}

