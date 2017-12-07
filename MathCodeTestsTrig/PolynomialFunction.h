#pragma once


#ifndef POLYNOMIALFUNCTION_H
#define POLYNOMIALFUNCTION_H


#include "FunctionEnums.h"
#include "MiscMathEquations.h"


#include <tuple>
#include <vector>

using std::tuple;
using std::vector;

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
	bool m_bIsContinuousFunction = false;
	bool m_bIsDifferientableEveryWhere = true; // all polynomials are
	bool m_bHasDomainBeenRestricted = false;



	// // https://www.varsitytutors.com/hotmath/hotmath_help/topics/end-behavior-of-a-function
	// The end behavior of a polynomial function is the behavior of the graph of f(x)f(x) as xx approaches positive infinity or negative infinity.
	EndBehavior m_EndBehavior; // not really using this one anymore

	FuncEndBehavior m_XGoesNegDir = std::make_pair(0.f, 0.f); // = std::make_pair(Interval(0, 0, IntervalType::IT_UNASSIGNED), Interval(0, 0, IntervalType::IT_UNASSIGNED));
	FuncEndBehavior m_XGoesPosDir = std::make_pair(0.f, 0.f);

	Domain m_Domain; // not really using this one anymore
	Range m_Range; // not really using this one anymore

	Interval m_DomainInterval = std::make_tuple(0.0f, 0.0f, IntervalType::IT_UNASSIGNED);
	Interval m_RangeInterval = std::make_tuple(0.0f, 0.0f, IntervalType::IT_UNASSIGNED);

	PolynomialFunctionType m_PolyFunctionType;

	double m_LeadingCoefficent = 0;

	// any value in its domain where its derivative is 0
	std::vector<double> m_CriticalPoints;
	// any value in its domain where its derivative is 0
	std::vector<Point> m_CriticalValueCordPoints;

	std::vector<Point> m_LocalMaximumPoints;
	std::vector<Point> m_LocalMinimumPoints;

	// unnassigned at the moment
	double m_AbsoluteMaximum = 0;
	double m_AbsoluteMinimum = 0;
	bool bHasAbsoluteMax = false;
	bool bHasAbsoluteMin = false;

	// Might use these later
	inline void SetEndBehaviorNegDir(const FuncEndBehavior& InBehav) { m_XGoesNegDir = InBehav; }
	inline void SetEndBehaviorPosDir(const FuncEndBehavior& InBehav) { m_XGoesPosDir = InBehav; }

	void SetDegree(const double& InDegree) { m_Degree = InDegree; }
	void SetLeadingCoefficent(const double& InVariable) { m_LeadingCoefficent = InVariable; }

	// This function needs to have the degree and leading coefficent set before its called on a polynomial function
	void SetPolynomialEndBehaviours();
	

	inline void SetIsEvenFunction(const bool& Input) { m_bIsEvenFunction = Input; }
	inline void SetIsContinuousFunction(const bool& Input) { m_bIsContinuousFunction = Input; }


	void SetAbsoluteMaximum(const double& InVal) { m_AbsoluteMaximum = InVal; }
	void SetAbsoluteMinimum(const double& InVal) { m_AbsoluteMinimum = InVal; }
	void SetHasAbsoluteMaximum(const bool& InRes) { bHasAbsoluteMax = InRes; }
	void SetHasAbsoluteMinimum(const bool& InRes) { bHasAbsoluteMin = InRes; }
	

	// any value in its domain where its derivative is 0
	virtual void FindCriticalPoints() = 0;



public:
	explicit PolynomialFunction() = default;

	PolynomialFunction(const PolynomialFunction&) = default;

	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }

	FuncEndBehavior GetEndBehaviourNegDir() const { return m_XGoesNegDir; }
	FuncEndBehavior GetEndBehaviourPosDir() const { return m_XGoesPosDir; }


	inline Interval GetDomainInterval() const { return m_DomainInterval; }
	inline Interval GetRangeInterval() const { return m_RangeInterval; }


	void PrintEndBehaviours() const;
	void PrintDomain() const;
	void PrintRange() const;

	inline void RestrictDomain(const Interval& NewDomainInterval) { m_bHasDomainBeenRestricted = true;  m_DomainInterval = NewDomainInterval; }
};


#endif