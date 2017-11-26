#pragma once


#ifndef LOGARITHMICFUNCTION_H
#define LOGARITHMICFUNCTION_H

#include "TranscendentalFunction.h"
#include "FunctionEnums.h"
#include "MathConstants.h" // for eulere

#include <exception>
#include <cmath>
#include <iostream>

using std::exception;
using std::log10;
using std::log2;
using std::log;
using std::cout;
using std::endl;

// Since log and exponential are inverses...
// ln(e^x) = x and e^ln(x) = x

// Rules: http://www.onlinemathlearning.com/image-files/logarithm-rules.png
// Properties: http://www.onlinemathlearning.com/image-files/log-properties.png
class LogarithmicFunction : TranscendentalFunction
{
private:


	//double m_a;
	double m_b;

	Domain m_Domain;
	Range m_Range;


public:

	explicit LogarithmicFunction(const double& b) //const double& b)
		: m_b(b) //m_a(a)
	{


		//// Is this the correct rules? (is this true?)
		//if (b == 1)//a == 1 || )
		//{
		//	
		//}

		if (b <= 0 || b == 1)
			throw exception("b value of log func has invalid value");

		if (b > 0 && b != 1)
		{
			m_Domain = Domain::ExclusiveZeroToPosInfinity;
			m_Range = Range::NegInfinityToPosInfinity;

			// satisfies: where log_b(x) = y if and only if b^y = x
			const double InputForTest = 8;
			double y = operator()(InputForTest);

			if (std::pow(m_b, y) == InputForTest)
			{
				std::cout << "This satisfies the conditions for a logrithmic function\n";
			}
		}


	}

	double operator()(const double x) const
	{
		if (m_b == 10)
			return std::log10(x);

		if (m_b == 2)
			return std::log2(x);

		// ln(x) == log_e(x)
		if (m_b == Eulere)
			return std::log(x);
	}

	inline void PrintPropertiesOfLogaritims()
	{
		std::cout << "If a, b, c > 0, and b != 1 and r is any real number, then\n";
		std::cout << "log_b(a*c) = log_b(a) + log_b(c) (Product property)\n";
		std::cout << "log_b(a / c) = log_b(a) - log_b(c) (Quotient property\n";
		std::cout << "log_b(a^r) = rlog_b(a) (Power property)\n";
	}

	// evaluate an equation with a non-standard base
	// log_a(x) = ln(x)/ln(x)
	// with b == e
	inline double UseChangeBaseFormula(const double& x, const double& a, const double& b = Eulere)
	{
		if (a < 0 || b < 0)
			throw std::exception("invalid change of base formula input < 0");

		if (a == 1 || b == 1)
			throw std::exception("invalid change of base formula input == 1");

		if (x < 0)
			throw std::exception("invalid change of base formula inputx < 0");

		double OutRes{ 0.0 };

		if (b == Eulere)
		{
			OutRes = std::log(x) / std::log(a);
		}

		return OutRes;
	}

	// earthquake 1 currently unused
	inline double EvaulateTwoEarthquakesRichterScale(
		const double& R1, const double& R2)
	{
		double LHS = R1 - R2;

		double OutRes = std::pow(10, LHS);
		//OutRes = OutRes * EarthQuake2;

		// OutRes currently means that 
		// earthquake1 was ___ times less/more intense than earthquake2

		return OutRes;

	}

};




#endif
