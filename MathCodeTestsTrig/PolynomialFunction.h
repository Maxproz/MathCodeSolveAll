#pragma once


#ifndef POLYNOMIALFUNCTION_H
#define POLYNOMIALFUNCTION_H


#include "FunctionEnums.h"
#include "MiscMathEquations.h"


#include <tuple>
#include <vector>

using std::tuple;
using std::vector;


// if a function has a local extremum at a critical point, then the sign of f′ switches as x increases through that point.


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

	

protected:
	int m_Degree = 0;
	double m_LeadingCoefficent = 0;
	double m_AbsoluteMaximum = 0;
	double m_AbsoluteMinimum = 0;
	bool m_bHasAbsoluteMax = false;
	bool m_bHasAbsoluteMin = false;
	bool m_bIsEvenFunction = false;					// set during SetPolynomialEndBehaviours()
	bool m_bIsLeadingCoefficentPositive = false;	// set during SetPolynomialEndBehaviours()
	bool m_bIsContinuousFunction = true;
	bool m_bIsDifferientableEveryWhere = true; // all polynomials are
	bool m_bHasDomainBeenRestricted = false;
	

	// TODO: Sloppy should fix this later, whats the point in two varibles
	// any value in its domain where its derivative is 0
	std::vector<double> m_AllCriticalPoints;
	

	std::vector<Point> m_AllPotentialExtrenums;

	std::vector<Point> m_CriticalPointsNotLocalMaxOrMins;
	std::vector<Point> m_LocalMaximumPoints;
	std::vector<Point> m_LocalMinimumPoints;


	// // https://www.varsitytutors.com/hotmath/hotmath_help/topics/end-behavior-of-a-function
	// The end behavior of a polynomial function is the behavior of the graph of f(x) as x approaches positive infinity or negative infinity.
	EndBehavior m_XGoesNegDir = std::make_pair(0.f, 0.f);
	EndBehavior m_XGoesPosDir = std::make_pair(0.f, 0.f);

	Interval m_DomainInterval = std::make_tuple(0.0f, 0.0f, IntervalType::IT_UNASSIGNED);
	Interval m_RangeInterval = std::make_tuple(0.0f, 0.0f, IntervalType::IT_UNASSIGNED);

	PolynomialFunctionType m_PolyFunctionType;

	std::vector<std::pair<double, double>> m_IncreasingIntervals;
	std::vector<std::pair<double, double>> m_DecreasingIntervals;


	// So these two functions are not used at the momemnt, but I might use these to set the endbehavior of a polynomial after I restriect the domain of it.
	inline void SetEndBehaviorNegDir(const EndBehavior& InBehav) { m_XGoesNegDir = InBehav; }
	inline void SetEndBehaviorPosDir(const EndBehavior& InBehav) { m_XGoesPosDir = InBehav; }


	inline void SetDegree(const double& InDegree) { m_Degree = InDegree; }
	inline void SetLeadingCoefficent(const double& InVariable) { m_LeadingCoefficent = InVariable; }
	// This function needs to have the degree and leading coefficent set before its called on a polynomial function
	
	// All polynomails are continuous over their domains, will need to change for rational functions
	inline void SetIsContinuousFunction(const bool& Input) { m_bIsContinuousFunction = Input; }
	
	inline void SetAbsoluteMaximum(const double& InVal) { m_AbsoluteMaximum = InVal; }
	inline void SetAbsoluteMinimum(const double& InVal) { m_AbsoluteMinimum = InVal; }
	inline void SetHasAbsoluteMaximum(const bool& InRes) { m_bHasAbsoluteMax = InRes; }
	inline void SetHasAbsoluteMinimum(const bool& InRes) { m_bHasAbsoluteMin = InRes; }
	
	void SetPolynomialEndBehaviours();

	// any value in its domain where its derivative is 0 (or undefined?)
	virtual void FindCriticalPoints() = 0;
	// These should probably be converted to virtual functions
	// THIS FUNCTION CAN ONLY BE CALLED AFTER THE DEGREE AND LEADING COEFFICENT HAVE BEEN SET
	
	//inline void SetIsEvenFunction(const bool& Input) { m_bIsEvenFunction = Input; }
	virtual void SetDefaultDomainInterval() = 0;
	virtual void SetDefaultRangeInterval() = 0;

	virtual void SetIncreasingDecreasingIntervals() = 0; 

public:
	explicit PolynomialFunction() = default;
	explicit PolynomialFunction(const PolynomialFunction&) = default;


	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }


	// End behaviors asssigned by SetPolynomialEndBehaviors()
	EndBehavior GetEndBehaviourNegDir() const { return m_XGoesNegDir; }
	EndBehavior GetEndBehaviourPosDir() const { return m_XGoesPosDir; }


	inline std::vector<Point> GetLocalMinimumPoints() const { return m_LocalMinimumPoints; }
	inline std::vector<Point> GetLocalMaximumPoints() const { return m_LocalMaximumPoints; }
	inline std::vector<Point> GetAllPotentialExtrenums() const { return m_AllPotentialExtrenums; }
	inline std::vector<double> GetAllCriticalPoints() const { return m_AllCriticalPoints; }


	inline Interval GetDomainInterval() const { return m_DomainInterval; }
	inline Interval GetRangeInterval() const { return m_RangeInterval; }


	void PrintEndBehaviours() const;
	void PrintDomain() const;
	void PrintRange() const;
	void PrintIncreasingDecreasingIntervals() const;

	inline void RestrictDomain(const Interval& NewDomainInterval) { m_bHasDomainBeenRestricted = true;  m_DomainInterval = NewDomainInterval; }

	void PrintRelatedCriticalPointData();

};


#endif