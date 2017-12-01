#pragma once


#ifndef POWERFUNCTION_H
#define POWERFUNCTION_H

#include "PolynomialFunction.h"

#include <cmath>
#include <tuple>
#include <iostream>

using std::pow;
using std::tuple;
using std::cout;
using std::endl;

// Hmm how do I want to handle this class....

// The power function f(x)=x^n is an even function if n is even and n≠0,
// and it is an odd function if  n  is odd

// y = a[k(x – d)]^n + c

// I already have functions defined for power functions with Exponents == 2 and 3
// So as of right now this is for polynomial functions with exponents > 3 

//
template <int Exponent>
class PowerFunction : PolynomialFunction
{
private:
	// assigned on Initialization
	double m_n;

	double m_a;
	double m_k;
	double m_d;
	double m_c;

public:

	PowerFunction() = default;

	PowerFunction(const PowerFunction&) = default;

	explicit PowerFunction(const double& a, const double& k, const double& d, const double& c)
		: m_a(a), m_k(k), m_d(d), m_c(c), m_n(Exponent)
	{
		m_PolyFunctionType = PolynomialFunctionType::POWER;



	}

	double operator()(const double& x) const
	{

		double First = x - m_d;
		double Second = First * m_k;
		double Third = std::pow(Second, m_n);
		double Fourth = Third * m_a;
		double Fifth = Fourth + m_c;

		return Fifth;
	}

	inline tuple<double, double, double, double, double> GetNAKDC() const
	{
		return tuple<double, double, double, double, double>(m_n, m_a, m_k, m_d, m_c);
	}

	inline void PrintFunction() const
	{
		// TODO: Edit this to display a better formatted text
		/*cout << "f(x) = ";*/
		bool bAddSecondBracket = false;

		bool bShouldShowParenthesis = true;

		if (m_d == 0)
		{
			bShouldShowParenthesis = false;
		}

		if (!(m_a == 1))
		{
			cout << m_a;
		}

		if (!(m_k == 1))
		{
			cout << "[" << m_k;
			bAddSecondBracket = true;
		}

		if (bShouldShowParenthesis)
		{
			cout << "(x";
		}
		else
		{
			cout << "x";
		}


		if (m_d == 0)
		{

		}
		else
		{
			cout << " - " << m_d;
		}

		if (bShouldShowParenthesis)
			cout << ")";

		if (bAddSecondBracket == true)
			cout << "]";

		cout << "^" << m_n;

		if (!(m_c == 0))
		{
			cout << " + " << m_c;
		}

	}

};


//template <>
//class PowerFunction<4> : PolynomialFunction
//{
//private:
//	// assigned on Initialization
//	double m_n;
//
//	double m_a;
//	double m_k;
//	double m_d;
//	double m_c;
//
//public:
//
//	PowerFunction() = default;
//
//	PowerFunction(const PowerFunction&) = default;
//
//	explicit PowerFunction(const double& a, const double& k, const double& d, const double& c)
//		: m_a(a), m_k(k), m_d(d), m_c(c), m_n(4)
//	{
//		m_PolyFunctionType = PolynomialFunctionType::POWER;
//
//
//
//	}
//
//	double operator()(const double& x) const
//	{
//
//		double First = x - m_d;
//		double Second = First * m_k;
//		double Third = std::pow(Second, m_n);
//		double Fourth = Third * m_a;
//		double Fifth = Fourth + m_c;
//
//		return Fifth;
//	}
//
//	inline tuple<double, double, double, double, double> GetNAKDC() const
//	{
//		return tuple<double, double, double, double, double>(m_n, m_a, m_k, m_d, m_c);
//	}
//
//	inline void PrintFunction() const
//	{
//		// TODO: Edit this to display a better formatted text
//		/*cout << "f(x) = ";*/
//		bool bAddSecondBracket = false;
//
//		bool bShouldShowParenthesis = true;
//
//		if (m_d == 0)
//		{
//			bShouldShowParenthesis = false;
//		}
//
//		if (!(m_a == 1))
//		{
//			cout << m_a;
//		}
//
//		if (!(m_k == 1))
//		{
//			cout << "[" << m_k;
//			bAddSecondBracket = true;
//		}
//
//		if (bShouldShowParenthesis)
//		{
//			cout << "(x";
//		}
//		else
//		{
//			cout << "x";
//		}
//	
//
//		if (m_d == 0)
//		{
//
//		}
//		else
//		{
//			cout << " - " << m_d;
//		}
//
//		if (bShouldShowParenthesis)
//			cout << ")";
//		
//		if (bAddSecondBracket == true)
//			cout << "]";
//			
//		cout << "^" << m_n;
//
//		if (!(m_c == 0))
//		{
//			cout << " + " << m_c;
//		}
//			
//	}
//
//};


#endif
