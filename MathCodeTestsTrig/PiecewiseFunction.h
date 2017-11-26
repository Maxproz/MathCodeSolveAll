#pragma once

// TODO: Need to go through these operator() functions and fix the logic in them. [look at lines (128 - 143) to refresh memory]

#ifndef PIECEWISEFUNCTION_H
#define PIECEWISEFUNCTION_H

#include "FunctionEnums.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::unique_ptr;
using std::vector;
using std::pair;
using std::isnan;


template <typename FirstFunc, typename SecondFunc>
class PiecewiseFunction
{
private:
	FirstFunc m_FirstFunc;
	SecondFunc m_SecondFunc;

	std::string m_FirstFuncRangeOp;
	std::string m_SecondFuncRangeOp;

	int m_RangeVariable{ 0 };


	std::vector<int> m_CurrentDiscontinousLocations;
	static int m_AmountOfDiscontinunitiesFound;

	// Pointer to a type and a point it is discontinous at
	std::unique_ptr<std::pair<DiscontinunityType, int>> m_Discontinuity = nullptr;

	void IncreaseDiscontinunitiesFound() const { m_AmountOfDiscontinunitiesFound++; }

	mutable double m_LastEvaluatedResult = 0.0;

public:


	PiecewiseFunction(const FirstFunc& InFirstFunc, const std::string& FirstFuncRangeOp,
		const SecondFunc& InSecondFunc, const std::string& SecondFuncRangeOp, const int& RangeVariable)
		: m_FirstFunc(InFirstFunc), m_FirstFuncRangeOp(FirstFuncRangeOp), m_SecondFunc(InSecondFunc),
		m_SecondFuncRangeOp(SecondFuncRangeOp), m_RangeVariable(RangeVariable)
	{

		// Temp Testing



	}

	//PiecewiseFunction(const PiecewiseFunction&) = default;

	//PiecewiseFunction<QuadraticFunction, LinearFunction>


	//PiecewiseFunction(const PiecewiseFunction<QuadraticFunction, LinearFunction>&) = default;

	double GetLastEvaluatedResult() const { return m_LastEvaluatedResult; }

	const int GetAmountOfDiscontinunitiesFound() const { return m_AmountOfDiscontinunitiesFound; }
	std::vector<int> GetFoundDiscontinunityLocations() const { return m_CurrentDiscontinousLocations; }

	void SetDiscontinunityPtr(const DiscontinunityType& Type, const int AtXEqualTo)

	{
		// Implicit move operation into the variable that stores the result.
		m_Discontinuity = std::make_unique<std::pair<DiscontinunityType, int>>(Type, AtXEqualTo);
	}

	std::pair<DiscontinunityType, int> GetCurrentDiscontinunityPtrInfo()
	{
		if (m_Discontinuity != nullptr)
		{
			DiscontinunityType Type = m_Discontinuity->first;
			int XLoc = m_Discontinuity->second;

			std::pair<DiscontinunityType, int> OutRes(Type, XLoc);

			return OutRes;
		}
	}

	// TODO: add functionality for a piecewise function that takes 3 functions


	double operator()(const double& x) const
	{

		double Res{ 0.0 };

		// Get the ranges with if else evaluate need to check all possiblites
		const double RangeVar = m_RangeVariable;

		// Check the x value options for the first function
		if (m_FirstFuncRangeOp == "<")
		{
			if (x < RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_FirstFuncRangeOp == ">")
		{
			if (x > RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_FirstFuncRangeOp == "<=")
		{

			// TODO: I probably want to make the rest of them like this too huh..
			if (x <= RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
			else
			{
				//return m_FirstFunc(x);
				Res = m_SecondFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_FirstFuncRangeOp == ">=")
		{
			if (x >= RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_FirstFuncRangeOp == "==")
		{
			if (x == RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_FirstFuncRangeOp == "!=")
		{
			if (x != RangeVar)
			{
				//return m_FirstFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}

		// Second Function options
		if (m_SecondFuncRangeOp == "<")
		{
			if (x < RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_SecondFuncRangeOp == ">")
		{
			if (x > RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_SecondFuncRangeOp == "<=")
		{
			if (x <= RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_SecondFuncRangeOp == ">=")
		{
			if (x >= RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);

				return Res;
			}

		}
		else if (m_SecondFuncRangeOp == "==")
		{
			if (x == RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}
		else if (m_SecondFuncRangeOp == "!=")
		{
			if (x != RangeVar)
			{
				//return m_SecondFunc(x);
				Res = m_FirstFunc(x);
				m_LastEvaluatedResult = Res;
				TestAndPushResForDiscontinunity(Res, x);
				return Res;
			}
		}


		return Res;

	}

	inline void TestAndPushResForDiscontinunity(const double& Res, const int& x) const
	{
		if (!(std::isnan(Res)))
		{
			GetFoundDiscontinunityLocations().push_back(x);
			IncreaseDiscontinunitiesFound();
		}
	}

	////template <typename FirstFunc>
	//FirstFunc GetFirstFunction() const { return m_FirstFunc; }
	//
	////template <typename SecondFunc>
	//SecondFunc GetSecondFunction() const { return m_SecondFunc; }

};



template <typename FirstFunc, typename SecondFunc, int ThirdFunctionConstant>
class PiecewiseFunctionThreeFunctions
{
private:
	FirstFunc m_FirstFunc;
	SecondFunc m_SecondFunc;
	//int m_ThirdFunctionReturn = ThirdFunctionConstant;

	std::string m_FirstFuncRangeOp;
	std::string m_SecondFuncRangeOp;
	std::string m_ThirdFuncRangeOp;

	int m_RangeVariable = { 0 };

public:



	explicit PiecewiseFunctionThreeFunctions(const FirstFunc& InFirstFunc, const std::string& FirstFuncRangeOp,
		const SecondFunc& InSecondFunc, const std::string& SecondFuncRangeOp, const int& RangeVariable,
		const std::string& ThirdFunctionRangeOp)
		: m_FirstFunc(InFirstFunc), m_FirstFuncRangeOp(FirstFuncRangeOp), m_SecondFunc(InSecondFunc),
		m_SecondFuncRangeOp(SecondFuncRangeOp), m_RangeVariable(RangeVariable), m_ThirdFuncRangeOp(ThirdFunctionRangeOp)
	{

		// Temp Testing



	}


	// TODO: add functionality for a piecewise function that takes 3 functions


	double operator()(const double& x) const
	{

		// Get the ranges with if else evaluate need to check all possiblites
		const double RangeVar = m_RangeVariable;

		// Check the x value options for the first function
		if (m_FirstFuncRangeOp == "<")
		{
			if (x < RangeVar)
			{
				return m_FirstFunc(x);
			}
		}
		else if (m_FirstFuncRangeOp == ">")
		{
			if (x > RangeVar)
			{
				return m_FirstFunc(x);
			}
		}
		else if (m_FirstFuncRangeOp == "<=")
		{
			if (x <= RangeVar)
			{
				return m_FirstFunc(x);
			}
		}
		else if (m_FirstFuncRangeOp == ">=")
		{
			if (x >= RangeVar)
			{
				return m_FirstFunc(x);
			}
		}
		else if (m_FirstFuncRangeOp == "==")
		{
			if (x == RangeVar)
			{
				return m_FirstFunc(x);
			}
		}
		else if (m_FirstFuncRangeOp == "!=")
		{
			if (x != RangeVar)
			{
				return m_FirstFunc(x);
			}
		}

		// Second Function options
		if (m_SecondFuncRangeOp == "<")
		{
			if (x < RangeVar)
			{
				return m_SecondFunc(x);
			}
		}
		else if (m_SecondFuncRangeOp == ">")
		{
			if (x > RangeVar)
			{
				return m_SecondFunc(x);
			}
		}
		else if (m_SecondFuncRangeOp == "<=")
		{
			if (x <= RangeVar)
			{
				return m_SecondFunc(x);
			}
		}
		else if (m_SecondFuncRangeOp == ">=")
		{
			if (x >= RangeVar)
			{
				return m_SecondFunc(x);
			}

		}
		else if (m_SecondFuncRangeOp == "==")
		{
			if (x == RangeVar)
			{
				return m_SecondFunc(x);
			}
		}
		else if (m_SecondFuncRangeOp == "!=")
		{
			if (x != RangeVar)
			{
				return m_SecondFunc(x);
			}
		}

		// Third Function options
		if (m_ThirdFuncRangeOp == "<")
		{
			if (x < RangeVar)
			{
				return ThirdFunctionConstant;
			}
		}
		else if (m_ThirdFuncRangeOp == ">")
		{
			if (x > RangeVar)
			{
				return ThirdFunctionConstant;
			}
		}
		else if (m_ThirdFuncRangeOp == "<=")
		{
			if (x <= RangeVar)
			{
				return ThirdFunctionConstant;
			}
		}
		else if (m_ThirdFuncRangeOp == ">=")
		{
			if (x >= RangeVar)
			{
				return ThirdFunctionConstant;
			}

		}
		else if (m_ThirdFuncRangeOp == "==")
		{
			if (x == RangeVar)
			{
				return ThirdFunctionConstant;
			}
		}
		else if (m_ThirdFuncRangeOp == "!=")
		{
			if (x != RangeVar)
			{
				return ThirdFunctionConstant;
			}
		}

	}

	////template <typename FirstFunc>
	//FirstFunc GetFirstFunction() const { return m_FirstFunc; }
	//
	////template <typename SecondFunc>
	//SecondFunc GetSecondFunction() const { return m_SecondFunc; }

};



#endif
