#pragma once


#ifndef ROOTFUNCTION_H
#define ROOTFUNCTION_H

// TODO: THIS FUNCTION NEEDS REFACTORED, A ROOT FUNCTION IS NOT A POLYNOMIAL FUNCTION!!!!!

#include "PolynomialFunction.h"

#include <exception>
#include <stdexcept>
#include <iostream>
#include <tuple>

using std::exception;	 // exception
using std::domain_error; // domain_error
using std::cout;
using std::endl;
using std::tuple;

inline bool isEven(int n) // add to helper function?
{
	if (n % 2 == 0)
		return true;
	else
		return false;
}


/*
You also learned that radical equations are equations with variables under square roots. 
Radical equations can be solved by isolating the square root and squaring both sides.
Sometimes radical equations will produce extraneous solutions, which are not really solutions,
so it is important to always check your answers to radical equations.

*/

////template <typename Root>
//class RootFunction : public PolynomialFunction
//{
//private:
//	double m_n;
//
//	double m_a;
//	double m_b;
//	double m_c;
//
//	// TODO: Set domain and range in constructor
//
//	// The domain goes to pos infinity in an even function 
//	// find the starting number by setting the expression inside of the root >= 0
//	double m_StartingDomainNum = 0;
//
//
//
//public:
//	RootFunction() = default;
//
//	explicit RootFunction(const double& n, const double& a,
//		const double& b, const double& c) : m_n(n), m_a(a), m_b(b), m_c(c)
//	{
//
//		// ROOT FUNCTIONS ARE NOT POLYNOMIAL FUNCTIONS AHHHHHHHHHHH
//		//m_PolyFunctionType = PolynomialFunctionType::ROOT;
//
//		if (n == 0)
//			throw exception("Pretty sure this is invalid algebra");
//
//		if (n >= 1 && !isEven(n))
//		{
//			m_Domain = Domain::NegInfinityToPosInfinity;
//			m_bIsEvenFunction = false;
//
//		}
//
//		if (n >= 2 && isEven(n))
//		{
//			m_Domain = Domain::InclusiveZeroToPosInfinity;
//		}
//
//		if (m_b < 0)
//		{
//			// in the form of x + m_b >= 0 ----- [-m_b, PosInfinity)
//			m_StartingDomainNum = m_b*(-1);
//		}
//		else if (m_b > 0)
//		{
//			// in the form x - m_B >= 0 ----- [m_b, PosInfinity)
//			m_StartingDomainNum = m_b;
//		}
//	}
//
//	double operator()(const double x) const
//	{
//		if (x <= 0)
//		{
//			throw domain_error("x has to be >= 0");
//		}
//
//		double Power = 1.0 / m_n;
//		double Base = x - m_b;
//		double RootRes = std::pow(Base, Power);
//
//
//		cout << "Power: " << Power << endl;
//		cout << "Base: " << Base << endl;
//		cout << "RootRes: " << RootRes << endl;
//		cout << "A: " << m_a << endl;
//		cout << "C: " << m_c << endl;
//		return (m_a * (RootRes)) + m_c;
//	}
//
//	inline tuple<double, double, double, double> GetNABC() const
//	{
//		return tuple<double, double, double, double>(m_n, m_a, m_b, m_c);
//	}
//
//	double GetStartingDomainNum() const { return m_StartingDomainNum; }
//
//};

// f(x)= a*root(x - b) + c
template <int Root>
class RootFunction : public PolynomialFunction
{
private:
	double m_n;

	double m_a;
	double m_b;
	double m_c;

	// TODO: Set domain and range in constructor

	// The domain goes to pos infinity in an even function 
	// find the starting number by setting the expression inside of the root >= 0
	double m_StartingDomainNum = 0;
	


public:
	RootFunction() = default;

	RootFunction(const RootFunction&) = default;

	explicit RootFunction(const double& a,const double& b, const double& c) 
		: m_n(Root), m_a(a), m_b(b), m_c(c)
	{

		// ROOT FUNCTIONS ARE NOT POLYNOMIAL FUNCTIONS AHHHHHHHHHHH
		//m_PolyFunctionType = PolynomialFunctionType::ROOT;

		if (m_n == 0)
			throw exception("Pretty sure this is invalid algebra");

		if (m_n >= 1 && !isEven(m_n))
		{
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_bIsEvenFunction = false;

		}

		if (m_n >= 2 && isEven(m_n))
		{
			m_Domain = Domain::InclusiveZeroToPosInfinity;
		}

		if (m_b < 0)
		{
			// in the form of x + m_b >= 0 ----- [-m_b, PosInfinity)
			m_StartingDomainNum = m_b*(-1);
		}
		else if (m_b > 0)
		{
			// in the form x - m_B >= 0 ----- [m_b, PosInfinity)
			m_StartingDomainNum = m_b;
		}
	}

	double operator()(const double x) const
	{
		if (x <= 0)
		{
			throw domain_error("x has to be >= 0");
		}

		double Power = 1.0 / m_n;
		double Base = x - m_b;
		double RootRes = std::pow(Base, Power);


		cout << "Power: " << Power << endl;
		cout << "Base: " << Base << endl;
		cout << "RootRes: " << RootRes << endl;
		cout << "A: " << m_a << endl;
		cout << "C: " << m_c << endl;
		return (m_a * (RootRes)) + m_c;
	}

	inline tuple<double, double, double, double> GetNABC() const
	{
		return tuple<double, double, double, double>(m_n, m_a, m_b, m_c);
	}

	double GetStartingDomainNum() const { return m_StartingDomainNum; }

	
	inline void PrintFunction() const
	{
		// TODO:
		// Edit this to display a better formatted text
		cout << "f(x) = " << m_a << "(x - " << m_b << ")^1/" << Root << " + " << m_c << endl;
	}

};



#endif

