#include "CalculusFunction.h"


void PrintPointSlopeForm(const double& Slope, const Point& InPoint)
{
	double X1 = InPoint.first;
	double Y1 = InPoint.second;

	std::cout << "y - " << Y1 << " = " << Slope << "(x - " << X1 << ")" << std::endl;

}

double GetSlope(const Point& FirstPoint, const Point& SecondPoint)
{
	double X1, X2, Y1, Y2;

	X1 = FirstPoint.first;
	Y1 = FirstPoint.second;
	X2 = SecondPoint.first;
	Y2 = SecondPoint.second;

	double Numerator = Y2 - Y1;
	double Denominator = X2 - X1;

	return Numerator / Denominator;
}

void PrintSlopeInterceptForm(const Point& Point, const double& Slope)
{
	double RHS = Point.first;
	double LHS = Point.second;

	if (LHS > 0)
	{
		LHS = LHS*-1;
	}

	if (RHS > 0)
	{
		RHS = RHS*-1;
	}

	RHS = Slope*RHS;


	if (LHS < 0)
	{
		RHS = RHS + LHS;
	}
	else
	{
		// LHS number is needs subtracted ( > 0)
		RHS = RHS - LHS;
	}

	double B = RHS;

	std::cout << "y = " << Slope << "x + " << B << std::endl;
}


//bool DetermineContinunityAtAPoint(RationalFunction& InFunc, const int& InPoint)
//{
//	//QuadraticFunction Numerator = InFunc.GetNumeratorQuadratic();
//	//LinearFunction Denominator = InFunc.GetDenominatorLinear();
//
//	// Step 1: check to see if f(a) is defined
//	double TestOne = InFunc.GetLastCalculatedRes();
//
//	std::cout << "TESTONE" <<  TestOne << std::endl;
//
//	if (std::isnan(TestOne))// || std::isinf(TestOne))
//	{
//		// failed 
//		std::cout << "TestOneFailed: The function is not continuous at " << InPoint << "\n";
//		return false;
//	}
//	else
//	{
//		//  If f(a) is defined, continue to step 2.
//
//		// Step 2: Compute Limit from both sides
//		// If Limit does not exist (that is, it is not a real number),
//		// then the function is not continuous at a and the problem is solved.
//
//
//
//		//Limit TestTwoLimit(InFunc, InPoint);
//		//double TestTwo = TestTwoLimit.GetLimitResult();
//
//		// if Limit Exists go to step 3
//		// TODO: Need a rational function check to return bool to see if it exists
//		// TODO: Need more example data in order to fix
//
//		Limit TestTwoLimit(InFunc, InPoint);
//		double TestTwoLimitNegDir = TestTwoLimit.GetLimitFromNegDir();
//		double TestTwoLimitPosDir = TestTwoLimit.GetLimitFromPosDir();
//
//		bool NegDirIsPosOrNegInf = (std::isinf(TestTwoLimitNegDir));
//		bool PosDirIsPosOrNegInf = (std::isinf(TestTwoLimitPosDir));
//
//		if (NegDirIsPosOrNegInf &&  PosDirIsPosOrNegInf)
//		{
//			TestTwoLimit.CheckAndSetRationalFuncDiscontinunities(InFunc, InPoint);
//
//			std::cout << "TestTwoFailed: The function is not continuous at " << InPoint << "\n";
//			return false;
//
//		}
//
//		//if (TestOne != TestTwo)
//		//{
//		//	std::cout << "TestThreeFailed: The function is not continuous at " << InPoint << "\n";
//		//	return false;
//		//}
//		//else
//		//{
//		//	std::cout << "TestThreePassed: The function is continuous at " << InPoint << "\n";
//		//	return true;
//		//}
//
//	}
//}


// TODO: Function needs edited after template changes
bool DetermineContinunityAtAPoint(PiecewiseFunction<QuadraticFunction, LinearFunction>& InFunc, const int& InPoint)
{
	////QuadraticFunction FirstFuntion = InFunc.GetFirstFunction();
	////LinearFunction SecondFunction = InFunc.GetSecondFunction();

	//// Step 1: check to see if f(a) is defined
	//double TestOne = InFunc.GetLastEvaluatedResult();

	//if (std::isnan(TestOne))
	//{
	//	// failed 
	//	std::cout << "TestOneFailed: The function is not continuous at " << InPoint << "\n";
	//	return false;
	//}
	//else
	//{
	//	//  If f(a) is defined, continue to step 2.

	//	// Step 2: Compute Limit from both sides
	//	 //If Limit does not exist (that is, it is not a real number),
	//	// then the function is not continuous at a and the problem is solved.



	//	//Limit TestTwoLimit(InFunc, InPoint);
	//	

	//	Limit<PiecewiseFunction<QuadraticFunction, LinearFunction>> TestTwoLimit(InFunc, InPoint);
	//	double TestTwoLimitNegDir = TestTwoLimit.GetLimitFromNegDir();
	//	double TestTwoLimitPosDir = TestTwoLimit.GetLimitFromPosDir();

	//	bool BothSidedLimitsAreEqual = (TestTwoLimitNegDir == TestTwoLimitPosDir);
	//		//TestTwoLimit.GetLimitResult();

	//	//std::cout << "TestTwo: = " << TestTwo << std::endl;

	//	double TestTwo = 0;

	//	if (BothSidedLimitsAreEqual == false)
	//	{
	//		TestTwoLimit.CheckAndSetPiecewiseFuncQuadLinearDiscontinunities(InFunc, InPoint);

	//		std::cout << "TestTwoFailed: The function is not continuous at " << InPoint << "\n";
	//		return false;

	//	}
	//	else
	//	{
	//		TestTwo = TestTwoLimit.GetLimitResult();
	//	}
	//	
	//	// if Limit Exists go to step 3
	//	// TODO: Need a rational function check to return bool to see if it exists
	//	// TODO: Need more example data in order to fix

	//	if (TestOne != TestTwo)
	//	{
	//		std::cout << "TestThreeFailed: The function is not continuous at " << InPoint << "\n";

	//		return false;
	//	}
	//	else
	//	{
	//		std::cout << "TestThreePassed: The function is continuous at " << InPoint << "\n";
	//		return true;
	//	}

	//}
	
	return false;
}



