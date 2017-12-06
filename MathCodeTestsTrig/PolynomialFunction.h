#pragma once


#ifndef POLYNOMIALFUNCTION_H
#define POLYNOMIALFUNCTION_H

#include "FunctionEnums.h"
#include "MiscMathEquations.h"

//#include "PowerFunction.h"
//
//#include <list>
//
//using std::list;




enum class PolynomialFunctionType
{
	LINEAR,
	RATIONAL,
	CUBIC,
	QUADRATIC,
	CONSTANT,
	POWER,
	QUARTIC,
	//ROOT, // not polynomial
};

class PolynomialFunction
{
private:

	//std::list<double> m_Coefficents;




protected:
	int m_Degree = 0;
	bool m_bIsEvenFunction = false;
	bool m_bIsLeadingCoefficentPositive = false;

	// // https://www.varsitytutors.com/hotmath/hotmath_help/topics/end-behavior-of-a-function
	// The end behavior of a polynomial function is the behavior of the graph of f(x)f(x) as xx approaches positive infinity or negative infinity.
	EndBehavior m_EndBehavior;

	FuncEndBehavior m_XGoesNegDir = std::make_pair(0.f, 0.f); // = std::make_pair(Interval(0, 0, IntervalType::IT_UNASSIGNED), Interval(0, 0, IntervalType::IT_UNASSIGNED));
	FuncEndBehavior m_XGoesPosDir = std::make_pair(0.f, 0.f);

	Domain m_Domain;
	Range m_Range;

	PolynomialFunctionType m_PolyFunctionType;

	double m_LeadingCoefficent = 0;

	// Might use these later
	inline void SetEndBehaviorNegDir(const FuncEndBehavior& InBehav) { m_XGoesNegDir = InBehav; }
	inline void SetEndBehaviorPosDir(const FuncEndBehavior& InBehav) { m_XGoesPosDir = InBehav; }

	void SetDegree(const double& InDegree) { m_Degree = InDegree; }
	void SetLeadingCoefficent(const double& InVariable) { m_LeadingCoefficent = InVariable; }

	// This function needs to have the degree and leading coefficent set before its called on a polynomial function
	inline void SetPolynomialEndBehaviours()
	{
		// As x goes to the negative direction is the first interval.
		//Interval AsXGoesToNegativeDirectionInterval = m_EndBehaviorIntervals.first;
		//float AsXGoesNegStartInterval = 0; 
		//float AsXGoesNegEndInterval = 0; // std::get<1>(AsXGoesToNegativeDirectionInterval);
		//IntervalType AsXGoesNegIntervalType = IntervalType::IT_UNASSIGNED; // std::get<2>(AsXGoesToNegativeDirectionInterval);

		//// As x goes to the positive direction is the Second interval.
		////Interval AsXGoesToPositiveDirectionInterval = m_EndBehaviorIntervals.second;
		//float AsXGoesPosStartInterval = 0;// std::get<0>(AsXGoesToPositiveDirectionInterval);
		//float AsXGoesPosEndInterval = 0; // std::get<1>(AsXGoesToPositiveDirectionInterval);
		//IntervalType AsXGoesPosIntervalType = IntervalType::IT_UNASSIGNED;// std::get<2>(AsXGoesToPositiveDirectionInterval);

		bool bIsLeadingCoefficentPositive = (m_LeadingCoefficent > 0);
		bool bIsFunctionDegreeEven = isEven(m_Degree);

		if (bIsLeadingCoefficentPositive)
		{
			// Leading coefficent is positive
			if (bIsFunctionDegreeEven)
			{
				//// leading coefficent positive and function degree even
				//float AsXGoesNegDirectionToStartInterval = 0; 
				//float FofXGoesToNegDir = 0; 

				m_XGoesNegDir = std::make_pair(NEGINFINITY, INFINITY);
				m_XGoesPosDir = std::make_pair(INFINITY, INFINITY);

				//IntervalType AsXGoesNegIntervalType = IntervalType::IT_OPEN;

				//float AsXGoesPosDirectionToEndInterval = 0;
				//float FofXGoesPosDir = 0;
				//IntervalType AsXGoesPosIntervalType = IntervalType::IT_UNASSIGNED;

			}
			else
			{
				// leading coefficent positive and function degree odd
				m_XGoesNegDir = std::make_pair(NEGINFINITY, NEGINFINITY);
				m_XGoesPosDir = std::make_pair(INFINITY, INFINITY);

			}
		}
		else
		{
			// Leading Coefficent is negative
			if (bIsFunctionDegreeEven)
			{
				// leading coefficent negative and function degree even

				m_XGoesNegDir = std::make_pair(NEGINFINITY, NEGINFINITY);
				m_XGoesPosDir = std::make_pair(INFINITY, NEGINFINITY);
			}
			else
			{
				// leading coefficent negative and function degree odd

				m_XGoesNegDir = std::make_pair(NEGINFINITY, INFINITY);
				m_XGoesPosDir = std::make_pair(INFINITY, NEGINFINITY);
			}
		}
	}

public:
	explicit PolynomialFunction() = default;

	PolynomialFunction(const PolynomialFunction&) = default;

	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }

	//explicit PolynomialFunction()
	//{
	//	int i = 0;
	//	while (HighestDegree - i != 0)
	//	{
	//		

	//		i = i + 1;
	//	}
	//}

	FuncEndBehavior GetEndBehaviourNegDir() const { return m_XGoesNegDir; }
	FuncEndBehavior GetEndBehaviourPosDir() const { return m_XGoesPosDir; }


	void PrintEndBehaviours() const;


};


#endif