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
private:

	ParabolaOpen m_ParabolaOpens;

	//vector<double> m_RealNumberZeros;
	int m_AmountOfRealNumberZeros;
	//AmountOfRealNumberSoltuions m_RealSolutionAmount;

	vector<double> m_AllZeros;

	double m_SymmetryLine;

	double m_MaxValueAtXIsEqualTo;
	double m_MinValueAtXIsEqualTo;

	void AutoSetHowManyRealZeroVariables();

	bool m_bJustABForm = false;
	bool m_bJustACForm = false;


	void PrintFunctionEndBehavior() const;
	void PrintParabolaOpensDirection() const;


	LinearFunction m_DerivativeFunction;// = LinearFunction(1, 0);

	

public:

	 QuadraticFunction() = default;

	 QuadraticFunction(const QuadraticFunction&) = default;

	//QuadraticFunction(QuadraticFunction&&) = default;

	//QuadraticFunction& operator =(const QuadraticFunction &) = default;

	explicit QuadraticFunction(double a, double b, double c = 0) : m_a(a), m_b(b), m_c(c)
	{
		m_PolyFunctionType = PolynomialFunctionType::QUADRATIC;

		if (a == 0)
			throw std::exception("a cannot == 0 for quadratic func initalization");

		
		// Qudratic degree is 2
		m_Degree = 2;

		// The degree of quadratic function is 2 so its even
		m_bIsEvenFunction = true;

		if (a > 0)
		{
			m_EndBehavior = EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToPosInfinity;
			m_ParabolaOpens = ParabolaOpen::UP;
		}
		else
		{
			// a < 0 at this point 
			m_EndBehavior = EndBehavior::AsXGoesToPosOrNegInfinityFOfXGoesToNegInfinity;
			m_ParabolaOpens = ParabolaOpen::DOWN;

		}

		if (c == 0)
		{
			m_bJustABForm = true;
		}

		if (b == 0)
		{
			m_bJustACForm = true;
		}

		AutoSetHowManyRealZeroVariables();

		SetZerosQuadraticFormula(*this);

		AutoSetDerivativeFunction(*this);
		
	}

	// assumes user has put in Vertex Form example y = a*(x - h)^2 + k
	// TODO: This is an awesome function, polish it and finish it up
	// TODO: need to add additional functionality for the different operator options
	// TODO: need to try different values for the "a" input and see what happens and add fixes 
	explicit QuadraticFunction(const std::string& FuncForm)
	{
		std::istringstream iss(FuncForm);

		double a = 0;
		char FirstParathensis;
		char x;
		char PlusOrMinusOpOne;
		double h;
		char PowerOp;
		char TwoChar;
		char SecondParathensis;
		char PlusOrMinusOpTwo;
		double k;

		// TODO: Add check for an implicit + k missing

		int MaybeA = iss.peek();
		if (MaybeA == EOF) return;
		if (!isdigit(MaybeA))
		{
			// Assume that we are dealing with an implicit 1 for the a variable
			a = 1;
			iss >> FirstParathensis >> x >> PlusOrMinusOpOne >> h >> PowerOp >> TwoChar >> SecondParathensis >>
				PlusOrMinusOpTwo >> k;
		}
		else
		{
			// otherwise read the function like normal
			iss >> a >> FirstParathensis >> x >> PlusOrMinusOpOne >> h >> PowerOp >> TwoChar >> SecondParathensis >>
				PlusOrMinusOpTwo >> k;
		}

		// TODO: remove debug code
		std::cout << a << std::endl;
		std::cout << FirstParathensis << std::endl;
		std::cout << x << std::endl;
		std::cout << PlusOrMinusOpOne << std::endl;
		std::cout << h << std::endl;
		std::cout << PowerOp << std::endl;
		std::cout << TwoChar << std::endl;
		std::cout << SecondParathensis << std::endl;
		std::cout << PlusOrMinusOpTwo << std::endl;
		std::cout << k << std::endl;

		double TermInOne, TermInTwo, TermInThree;

		if (PlusOrMinusOpOne == '-')
		{
			if (PlusOrMinusOpTwo == '+')
			{
				// y = a*(x-h)^2 + k
				TermInOne = a;
				TermInTwo = (((h) * (-1)) * 2);
				TermInThree = (std::pow(h, 2) * (a)) + k;

				m_a = TermInOne;
				m_b = TermInTwo;
				m_c = TermInThree;

			}
		}

		// make a helper function to handle all the misc constuctor stuff that is used in the main constuctor above

		AutoSetHowManyRealZeroVariables();

		SetZerosQuadraticFormula(*this);

		AutoSetDerivativeFunction(*this);

	}


	// multiplication
	QuarticFunction operator*(QuadraticFunction const& rhs) const;
	QuadraticFunction operator-(QuadraticFunction const& rhs) const;
	
	// The operator overloaded directly below didnt work
	//QuadraticTrigometric operator*(const TrigometricFunction<MPSIN>& rhs) const;
	/*QuarticFunction& operator*=(CubicFunction const& rhs);*/

	double operator()(const double x) const
	{
		double TermOne, TermTwo, TermThree;
		TermOne = (m_a * (std::pow(x, 2)));
		TermTwo = x*m_b;
		TermThree = m_c;
		return TermOne + TermTwo + TermThree;
	}

	inline std::tuple<double, double, double> GetABC() const
	{
		return std::tuple<double, double, double>(m_a, m_b, m_c);
	}

	void SetAllZeroVec(std::vector<double> InVec) { m_AllZeros = InVec; }
	//void SetRealNumberZeroVec(std::vector<double> InVec) { m_RealNumberZeros = InVec; }
	void SetLineOfSymmetry(double InNum) { m_SymmetryLine = InNum; }
	void SetTheMaxMinValue(double InNum);


	void PrintAllZeros() const;
	//void PrintAllRealNumberZeros() const;
	void PrintNumberOfRealNumberSoltions() const;

	void PrintBasicFunctionInfo() const;
	void PrintFunction() const;


	inline bool IsABForm() const { return m_bJustABForm; }
	inline bool IsACForm() const { return m_bJustACForm; }


	inline std::tuple<double, double> GetAB() const
	{
		return std::tuple<double, double>(m_a, m_b);
	}

	inline std::tuple<double, double> GetAC() const
	{
		return std::tuple<double, double>(m_a, m_c);
	}

	//std::vector<double> GetRealNumberZerosVec() const { return m_RealNumberZeros; }
	//inline int GetAmountOfRealZeros() const { return m_AmountOfRealNumberZeros; }

	std::vector<double> GetAllZerosVec() const { return m_AllZeros; }
	

	LinearFunction GetDerivativeFunction() const;
	void SetDerivativeFunction(LinearFunction& InFunc); 

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

