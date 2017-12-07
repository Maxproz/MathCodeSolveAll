#pragma once

#ifndef QUADRATICFUNCTION_H
#define QUADRATICFUNCTION_H


#include "PolynomialFunction.h"

#include "LinearFunction.h"

//#include "Derivative.h"


#include <vector>
#include <cctype> // isdigit()
#include <tuple>
#include <sstream> // istringstream
#include <iostream>


using std::vector;
//using std::isdigit;
using std::tuple;
using std::istringstream;
using std::cout;
using std::endl;


class QuarticFunction;
class LinearFunction;
class QuadraticFunction;
//class MPSIN;
//class MPCOS;






void SetZerosQuadraticFormula(QuadraticFunction& QuadraticFunc);

std::vector<double> GetZerosQuadraticFormula(const double& a, const double& b, const double& c);

void AutoSetDerivativeFunction(QuadraticFunction& InFunc);

// quadratic function has the form f(x)=ax^2+bx+c,where a≠0
class QuadraticFunction : public PolynomialFunction
{
private:

	friend class ConstantFunction;

protected:

	double m_a;
	double m_b;
	double m_c;

	double m_VertexX;
	double m_VertexY;

private:

	ParabolaOpen m_ParabolaOpens;

	//vector<double> m_RealNumberZeros;
	int m_AmountOfRealNumberZeros;
	//AmountOfRealNumberSoltuions m_RealSolutionAmount;

	vector<double> m_AllZeros;

	double m_SymmetryLine;

	// no idea what these variables are tracking
	double m_MaxValueAtXIsEqualTo;
	double m_MinValueAtXIsEqualTo;

	void AutoSetHowManyRealZeroVariables();

	bool m_bIsFunctionGeneralForm = false;
	bool m_bIsFunctionVertexForm = false;


	void PrintFunctionEndBehavior() const;
	void PrintParabolaOpensDirection() const;


	LinearFunction m_DerivativeFunction;// = LinearFunction(1, 0);

	// Open Down means function is upside down bowl shape
	inline void SetParabolaOpeningDirection(const double& LeadingCoefficent);

	void AutoSetDomainInterval();
	void AutoSetRangeInterval();

	void AutoSetVertexFromGeneralForm();
	void AutoSetVertexFromVertexForm(const double& h, const double& k);


	//virtual void FindCriticalPoints() override;

	//void FindGlobalExtremums();

public:

	virtual void FindCriticalPoints() override;

	void FindGlobalExtremums();

	 QuadraticFunction() = default;
	 QuadraticFunction(const QuadraticFunction&) = default;
	//QuadraticFunction(QuadraticFunction&&) = default;
	//QuadraticFunction& operator =(const QuadraticFunction &) = default;

	explicit QuadraticFunction(double a, double b, double c = 0) : m_a(a), m_b(b), m_c(c)
	{
		m_PolyFunctionType = PolynomialFunctionType::QUADRATIC;

		if (a == 0)
			throw std::exception("a cannot == 0 for quadratic func initalization (general form constructor)");

		AutoSetVertexFromGeneralForm();
		AutoSetDomainInterval();
	

		// Qudratic degree is 2 and also even 
		SetDegree(2); 
		SetLeadingCoefficent(m_a);
		SetPolynomialEndBehaviours();
		SetIsEvenFunction(true);
		SetParabolaOpeningDirection(a);
		m_bIsFunctionGeneralForm = true;

		// This needs to be called after ParabolaOpensDirection and After the VertexYCord has been set
		AutoSetRangeInterval();

		// Normal polynomials are continuous 
		SetIsContinuousFunction(true);

		AutoSetHowManyRealZeroVariables();

		SetZerosQuadraticFormula(*this);

		AutoSetDerivativeFunction(*this);
		

		//FindCriticalPoints();


	}

	// assumes user has put in Vertex Form example y = a*(x - h)^2 + k
	// TODO: This is an awesome function, polish it and finish it up
	// TODO: need to add additional functionality for the different operator options
	// TODO: need to try different values for the "a" input and see what happens and add fixes 
	explicit QuadraticFunction(const std::string& FuncForm);


	// multiplication
	QuarticFunction operator*(QuadraticFunction const& rhs) const;
	QuadraticFunction operator-(QuadraticFunction const& rhs) const;
	
	// The operator overloaded directly below didnt work
	//QuadraticTrigometric operator*(const TrigometricFunction<MPSIN>& rhs) const;
	/*QuarticFunction& operator*=(CubicFunction const& rhs);*/

	double operator()(const double& x) const
	{
		double TermOne, TermTwo, TermThree;
		TermOne = (m_a * (std::pow(x, 2)));
		TermTwo = x*m_b;
		TermThree = m_c;
		return TermOne + TermTwo + TermThree;
	}

	inline std::tuple<double, double, double> GetABC() const { return std::tuple<double, double, double>(m_a, m_b, m_c); }
	inline std::tuple<double, double> GetAC() const { return std::tuple<double, double>(m_a, m_c); }

	void SetAllZeroVec(std::vector<double> InVec) { m_AllZeros = InVec; }
	//void SetRealNumberZeroVec(std::vector<double> InVec) { m_RealNumberZeros = InVec; }
	void SetLineOfSymmetry(double InNum) { m_SymmetryLine = InNum; }
	void SetTheMaxMinValue(double InNum);


	void PrintAllZeros() const;
	//void PrintAllRealNumberZeros() const;
	void PrintNumberOfRealNumberSoltions() const;

	void PrintBasicFunctionInfo() const;
	void PrintFunction() const;
	std::string GetFunctionString() const;

	

	//std::vector<double> GetRealNumberZerosVec() const { return m_RealNumberZeros; }
	//inline int GetAmountOfRealZeros() const { return m_AmountOfRealNumberZeros; }

	std::vector<double> GetAllZerosVec() const { return m_AllZeros; }
	

	LinearFunction GetDerivativeFunction() const;
	void SetDerivativeFunction(const LinearFunction& InFunc);

};

template <typename T>
bool is_close_to_zero(T x)
{
	return std::abs(x) < std::numeric_limits<T>::epsilon();
}




////inline bool isEven(int n);
//// A QuadraticFunction in the form
//// ax^2TrigA(x) + bxTrigB(x) + TrigC(x)
//class QuadraticTrigometric : public QuadraticFunction
//{
//private:
//
//	TrigometricFunction<MPSIN> m_TrigA = TrigometricFunction<MPSIN>(1, 1, 0, 0);
//	//TrigometricFunction<BType> m_TrigB;
//	//TrigometricFunction<CType> m_TrigC;
//
//public:
//
//	explicit QuadraticTrigometric(const double& a, const TrigometricFunction<MPSIN>& triga, const double& b, const double& c)
//		: m_TrigA(triga)/*, m_TrigB(trigb), m_TrigC(trigc)*/
//	{
//		m_a = a;
//		m_b = b;
//		m_c = c;
//
//	}
//
//	QuadraticTrigometric(const QuadraticTrigometric&) = default;
//};
//


#endif

