#pragma once


#ifndef RATIONALFUNCTION_H
#define RATIONALFUNCTION_H


#include "PolynomialFunction.h"
#include "QuadraticFunction.h"
#include "LinearFunction.h"
#include "CubicFunction.h"
#include "RootFunction.h"

#include <memory>

using std::unique_ptr;
using std::vector;
using std::pair;
using std::make_unique;



class RationalFunction : PolynomialFunction
{
private:

	// Possible Numerator Functions
	QuadraticFunction m_NumeratorQuadratic;
	RootFunction m_NumeratorRoot; // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	LinearFunction m_NumeratorLinear; // denominator func cannot be ==  to 0


									  // Possible Denominator Functions
	CubicFunction m_DenominatorCubic; // denominator func cannot be == to 0
	LinearFunction m_DenominatorLinear; // denominator func cannot be ==  to 0
	QuadraticFunction m_DenominatorQuadratic; // denominator func cannot be ==  to 0
	RootFunction m_DenominatorRoot; // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH

	PolynomialFunctionType m_NumeratorFuncType;
	PolynomialFunctionType m_DenominatorFuncType;


	std::vector<int> m_CurrentDiscontinousLocations;
	static int m_AmountOfDiscontinunitiesFound;

	// Pointer to a type and a point it is discontinous at
	std::unique_ptr<std::pair<DiscontinunityType, int>> m_Discontinuity = nullptr;

	void IncreaseDiscontinunitiesFound() const { m_AmountOfDiscontinunitiesFound++; }

	mutable double m_LastCalculatedRes = 0;

public:

	RationalFunction() = default;
	RationalFunction(const RationalFunction&) = default;
	RationalFunction(RationalFunction&&) = default;

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





	explicit RationalFunction(const QuadraticFunction& Numerator, const CubicFunction& Denominator)
		: m_NumeratorQuadratic(Numerator), m_DenominatorCubic(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	explicit RationalFunction(const QuadraticFunction& Numerator, const LinearFunction& Denominator)
		: m_NumeratorQuadratic(Numerator), m_DenominatorLinear(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	explicit RationalFunction(const QuadraticFunction& Numerator, const QuadraticFunction& Denominator)
		: m_NumeratorQuadratic(Numerator), m_DenominatorQuadratic(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	// Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	explicit RationalFunction(const RootFunction& Numerator, const LinearFunction& Denominator)
		: m_NumeratorRoot(Numerator), m_DenominatorLinear(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	// Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	explicit RationalFunction(const LinearFunction& Numerator, const RootFunction& Denominator)
		: m_NumeratorLinear(Numerator), m_DenominatorRoot(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	explicit RationalFunction(const LinearFunction& Numerator, const LinearFunction& Denominator)
		: m_NumeratorLinear(Numerator), m_DenominatorLinear(Denominator)
	{
		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	}

	inline PolynomialFunctionType GetNumeratorFunctionType() const
	{
		return m_NumeratorFuncType;
	}

	inline PolynomialFunctionType GetDenominatorFunctionType() const
	{
		return m_DenominatorFuncType;
	}

	double CallingRationalFunctionFuncOperator(const double& x)
	{
		return operator()(x);
	}

	double operator()(const double& x)
	{
		double NumeratorRes{ 0.0 };
		double DenominatorRes{ 0.0 };

		switch (GetNumeratorFunctionType())
		{
			case PolynomialFunctionType::QUADRATIC:
			{
				NumeratorRes = m_NumeratorQuadratic(x);
				break;
			}
			case PolynomialFunctionType::LINEAR:
			{
				NumeratorRes = m_NumeratorLinear(x);
				break;
			}
		}

		switch (GetDenominatorFunctionType())
		{
			case PolynomialFunctionType::LINEAR:
			{
				DenominatorRes = m_DenominatorLinear(x);
				break;
			}
		}



		m_LastCalculatedRes = NumeratorRes / DenominatorRes;

		if (!(std::isnan(m_LastCalculatedRes)))
		{
			GetFoundDiscontinunityLocations().push_back(x);
			IncreaseDiscontinunitiesFound();
		}


		return NumeratorRes / DenominatorRes;
	}

public:



	double GetLastCalculatedRes() const { return m_LastCalculatedRes; }

	//friend double CallMemberFunctionFuncOperator(const double& x);

	// Possible NumeratorGetters
	QuadraticFunction GetNumeratorQuadratic() const { return m_NumeratorQuadratic; }
	RootFunction GetNumeratorRoot() const { return m_NumeratorRoot; } // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	LinearFunction GetNumeratorLinear() const { return m_NumeratorLinear; }

	// Possible Denominator Getters
	CubicFunction GetDenominatorCubic() const { return m_DenominatorCubic; }
	LinearFunction GetDenominatorLinear() const { return m_DenominatorLinear; }
	QuadraticFunction GetDenominatorQuadratic() const { return m_DenominatorQuadratic; }
	RootFunction GetDenominatorRoot() const { return m_DenominatorRoot; } // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH

};



#endif
