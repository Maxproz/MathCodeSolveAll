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

//template< class T >
//std::unique_ptr<T> copy_unique(const std::unique_ptr<T>& source)
//{
//	return source ? std::make_unique<T>(*source) : nullptr;
//}

template <typename NumFunc, typename DenomFunc>
class RationalFunction : public PolynomialFunction
{
private:


	NumFunc m_NumeratorFunction;
	DenomFunc m_DenominatorFunction;

	//// Possible Numerator Functions
	//QuadraticFunction m_NumeratorQuadratic;
	//RootFunction m_NumeratorRoot; // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	//LinearFunction m_NumeratorLinear; // denominator func cannot be ==  to 0


	//								  // Possible Denominator Functions
	//CubicFunction m_DenominatorCubic; // denominator func cannot be == to 0
	//LinearFunction m_DenominatorLinear; // denominator func cannot be ==  to 0
	//QuadraticFunction m_DenominatorQuadratic; // denominator func cannot be ==  to 0
	//RootFunction m_DenominatorRoot; // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH

	PolynomialFunctionType m_NumeratorFuncType;
	PolynomialFunctionType m_DenominatorFuncType;


	std::vector<int> m_CurrentDiscontinousLocations;
	int m_AmountOfDiscontinunitiesFound = 0;

	//// Pointer to a type and a point it is discontinous at
	//std::unique_ptr<std::pair<DiscontinunityType, int>> m_Discontinuity = nullptr;

	void IncreaseDiscontinunitiesFound() { m_AmountOfDiscontinunitiesFound++; }

	mutable double m_LastCalculatedRes = 0;

public:

	////template<typename NumFunc, typename DenomFunc>
	RationalFunction<NumFunc,DenomFunc>() = default;

	////template<typename NumFunc, typename DenomFunc>
	//RationalFunction& operator=(const RationalFunction<NumFunc, DenomFunc> &) = default;

	//template <typename NumFunc, typename DenomFunc>
	RationalFunction<NumFunc, DenomFunc>(const RationalFunction<NumFunc, DenomFunc>&) = default;
	
	//RationalFunction() : a_(new int(33)) {}

	//RationalFunction<NumFunc, DenomFunc>(RationalFunction<NumFunc, DenomFunc> &&data)
	//	: m_Discontinuity(std::move(data.m_Discontinuity = nullptr;))
	//{
	//	m_NumeratorFunction = data.m_NumeratorFunction;
	//	m_DenominatorFunction = data.m_DenominatorFunction;
	//	m_NumeratorFuncType = data.m_NumeratorFuncType;
	//	m_DenominatorFuncType = data.m_DenominatorFuncType;
	//	m_CurrentDiscontinousLocations = data.m_CurrentDiscontinousLocations;
	//	// TODO: how do i move the amount of discontinuties found?
	//	m_AmountOfDiscontinunitiesFound = data.m_AmountOfDiscontinunitiesFound;
	//	m_LastCalculatedRes = data.m_LastCalculatedRes;
	//}

	//RationalFunction<NumFunc, DenomFunc>& operator=(RationalFunction<NumFunc,DenomFunc> &&data)
	//{
	//	m_Discontinuity = std::move(data.m_Discontinuity);
	//	m_NumeratorFunction = data.m_NumeratorFunction;
	//	m_DenominatorFunction = data.m_DenominatorFunction;
	//	m_NumeratorFuncType = data.m_NumeratorFuncType;
	//	m_DenominatorFuncType = data.m_DenominatorFuncType;
	//	m_CurrentDiscontinousLocations = data.m_CurrentDiscontinousLocations;
	//	m_AmountOfDiscontinunitiesFound = data.m_AmountOfDiscontinunitiesFound;
	//	m_LastCalculatedRes = data.m_LastCalculatedRes;

	//	return *this;
	//}



	////template<typename NumFunc, typename DenomFunc>
	//RationalFunction(RationalFunction<NumFunc, DenomFunc>&&) = default;

	const int GetAmountOfDiscontinunitiesFound() const { return m_AmountOfDiscontinunitiesFound; }
	std::vector<int> GetFoundDiscontinunityLocations() const { return m_CurrentDiscontinousLocations; }

	void SetDiscontinunityPtr(const DiscontinunityType& Type, const int AtXEqualTo)

	{
		//// Implicit move operation into the variable that stores the result.
		//m_Discontinuity = std::make_unique<std::pair<DiscontinunityType, int>>(Type, AtXEqualTo);
	}

	std::pair<DiscontinunityType, int> GetCurrentDiscontinunityPtrInfo()
	{
		//if (m_Discontinuity != nullptr)
		//{
		//	DiscontinunityType Type = m_Discontinuity->first;
		//	int XLoc = m_Discontinuity->second;

			std::pair<DiscontinunityType, int> OutRes(Type, XLoc);

			return OutRes;
		
	}


	//RationalFunction<NumFunc, DenomFunc>(const RationalFunction<NumFunc, DenomFunc>&) = default;
	//RationalFunction<NumFunc,DenomFunc>& operator=(const RationalFunction<NumFunc, DenomFunc>&) = default;
	//RationalFunction(const RationalFunction<NumFunc, DenomFunc>&) = default;

	//RationalFunction(const RationalFunction<NumFunc, DenomFunc> &) = default;

	explicit RationalFunction(const NumFunc& Numerator, const DenomFunc& Denominator)
		: m_NumeratorFunction(Numerator), m_DenominatorFunction(Denominator)
	{

		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();

	}

	//RationalFunction(const RationalFunction&) = delete;
	//RationalFunction& operator =(const RationalFunction&) = delete;

	//explicit RationalFunction(const QuadraticFunction& Numerator, const CubicFunction& Denominator)
	//	: m_NumeratorQuadratic(Numerator), m_DenominatorCubic(Denominator)
	//{


	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();


	//}

	//explicit RationalFunction(const QuadraticFunction& Numerator, const LinearFunction& Denominator)
	//	: m_NumeratorQuadratic(Numerator), m_DenominatorLinear(Denominator)
	//{
	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	//}

	//explicit RationalFunction(const QuadraticFunction& Numerator, const QuadraticFunction& Denominator)
	//	: m_NumeratorQuadratic(Numerator), m_DenominatorQuadratic(Denominator)
	//{
	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	//}

	//// Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	//explicit RationalFunction(const RootFunction& Numerator, const LinearFunction& Denominator)
	//	: m_NumeratorRoot(Numerator), m_DenominatorLinear(Denominator)
	//{
	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	//}

	//// Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	//explicit RationalFunction(const LinearFunction& Numerator, const RootFunction& Denominator)
	//	: m_NumeratorLinear(Numerator), m_DenominatorRoot(Denominator)
	//{
	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	//}

	//explicit RationalFunction(const LinearFunction& Numerator, const LinearFunction& Denominator)
	//	: m_NumeratorLinear(Numerator), m_DenominatorLinear(Denominator)
	//{
	//	m_PolyFunctionType = PolynomialFunctionType::RATIONAL;

	//	m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
	//	m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
	//}

	inline PolynomialFunctionType GetNumeratorFunctionType() const
	{
		return m_NumeratorFuncType;
	}

	inline PolynomialFunctionType GetDenominatorFunctionType() const
	{
		return m_DenominatorFuncType;
	}

	//double CallingRationalFunctionFuncOperator(const double& x)
	//{
	//	return operator()(x);
	//}

	double operator()(const double& x)
	{
		double NumeratorRes{ 0.0 };
		double DenominatorRes{ 0.0 };

		NumeratorRes = m_NumeratorFunction(x);
		DenominatorRes = m_DenominatorFunction(x);

		//switch (GetNumeratorFunctionType())
		//{
		//	case PolynomialFunctionType::QUADRATIC:
		//	{
		//		QuadraticFunction Numer = 

		//		NumeratorRes = m_NumeratorQuadratic(x);
		//		break;
		//	}
		//	case PolynomialFunctionType::LINEAR:
		//	{
		//		NumeratorRes = m_NumeratorLinear(x);
		//		break;
		//	}
		//}

		//switch (GetDenominatorFunctionType())
		//{
		//	case PolynomialFunctionType::LINEAR:
		//	{
		//		DenominatorRes = m_DenominatorLinear(x);
		//		break;
		//	}
		//}



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

	NumFunc GetNumeratorFunction() const { return m_NumeratorFunction; }
	DenomFunc GetDenominatorFunction() const { return m_DenominatorFunction; }

	//// Possible NumeratorGetters
	//QuadraticFunction GetNumeratorQuadratic() const { return m_NumeratorQuadratic; }
	//RootFunction GetNumeratorRoot() const { return m_NumeratorRoot; } // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH
	//LinearFunction GetNumeratorLinear() const { return m_NumeratorLinear; }

	//// Possible Denominator Getters
	//CubicFunction GetDenominatorCubic() const { return m_DenominatorCubic; }
	//LinearFunction GetDenominatorLinear() const { return m_DenominatorLinear; }
	//QuadraticFunction GetDenominatorQuadratic() const { return m_DenominatorQuadratic; }
	//RootFunction GetDenominatorRoot() const { return m_DenominatorRoot; } // Rational FUNCTIONS CANNOT HAVE ROOT FUNCTIONS AHHHHHHHHHHH


	inline void PrintFunction() const
	{
		m_NumeratorFunction.PrintFunction();
		cout << " / ";
		m_DenominatorFunction.PrintFunction();
	}

};

//template<>
//class RationalFunction<QuadraticFunction, LinearFunction> : PolynomialFunction
//{
//private:
//
//
//public:
//
//	explicit RationalFunction(const QuadraticFunction& Numerator, const LinearFunction& Denominator)
//		: m_NumeratorFunction(Numerator), m_DenominatorFunction(Denominator)
//	{
//
//		m_PolyFunctionType = PolynomialFunctionType::RATIONAL;
//
//		m_NumeratorFuncType = Numerator.GetCurrentFunctionType();
//		m_DenominatorFuncType = Denominator.GetCurrentFunctionType();
//
//	}
//
//};


#endif
