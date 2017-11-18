#pragma once


#include <vector>
#include <cmath>
#include <utility>
#include <string>
#include <sstream>
#include <functional>
#include <iostream>
#include <tuple>
#include <limits>
#include <iomanip>
#include <algorithm>
#include "Circle.h"

#define NEGINFINITY   ((float)(_HUGE_ENUF * _HUGE_ENUF) *(-1))

const double Eulere = 2.718282;


using std::pair;

typedef pair<double, double> Point;

long mygcd(long a, long b);

inline std::pair<double, double> OutputDecimalAsFract(const double& input)
{
	double integral = std::floor(input);
	double frac = input - integral;

	const long precision = 10; // This is the accuracy. // was 1 million or so

	long gcd_ = mygcd(round(frac * precision), precision);

	long denominator = precision / gcd_;
	long numerator = round(frac * precision) / gcd_;

	/*
	if (integral != 0)
	{
	std::cout << integral << " + ";
	}*/
	if (input < 0)
		numerator = numerator*-1;

	//std::cout << numerator << "/" << denominator << std::endl;
	return std::pair<double, double>(numerator, denominator);
}



// Inclusive uses bracket notation
enum class Domain
{
	NegInfinityToPosInfinity,
	InclusiveZeroToPosInfinity,
	ExclusiveZeroToPosInfinity,
	Reciprocal,
	ExactGiven
};

enum class Range
{
	NegInfinityToPosInfinity,
	InclusiveZeroToPosInfinity,
	ExclusiveZeroToPosInfinity,
	Reciprocal,
	ExactGiven
};

enum class LineBehavior
{
	Increasing,
	Decreasing,
	Horizontal
};

enum class EndBehavior
{
	// The two behaviors for quadratic functions
	AsXGoesToPosOrNegInfinityFOfXGoesToPosInfinity,
	AsXGoesToPosOrNegInfinityFOfXGoesToNegInfinity,

	// TODO: maybe add more behaviors for other more functions later(cubic etc)

	// Exponential Func Behaviors
	// b > 1
	AsXGoesToPosInfinityFOfXGoesToPosInfinity,
	AsXGoesToNegInfinityFOfXGoesToZero,
	// 0 < b < 1
	AsXGoesToPosInfinityFOfXGoesToZero,
	AsXGoesToNegInfinityFOfXGoesToPosInfinity,
};


enum class ParabolaOpen
{
	UP,
	DOWN
};


// Basic Important Notes on Functions
/* 1. (algebraic functions can only have powers that are rational numbers.)

*/



void PrintPointSlopeForm(const double& Slope, const Point& InPoint);

double GetSlope(const Point& FirstPoint, const Point& SecondPoint);

void PrintSlopeInterceptForm(const Point& Point, const double& Slope);

template<typename T>
double evaluate_at(double x, const T& Function)
{
	return Function(x);
}

enum class PolynomialFunctionType
{
	LINEAR,
	RATIONAL,
	CUBIC,
	QUADRATIC,
	ROOT,
};

class PolynomialFunction
{
private:
	

protected:
	int m_Degree;
	bool m_bIsEvenFunction;
	EndBehavior m_EndBehavior;

	Domain m_Domain;
	Range m_Range;

	PolynomialFunctionType m_PolyFunctionType;

public:

	PolynomialFunctionType GetCurrentFunctionType() const { return m_PolyFunctionType; }

};


// The function of the coefficient aa in the general equation is to make the graph 
// "wider" or "skinnier", or to reflect it (if negative):
// 
//  f(x)=ax^3+bx^2+cx+d

// its natural domain and range including all real numbers x and y respectively. 
// TODO: Firstly, a cubic equation always has its center point on the vertical line x = -b/3*a
// TODO: when the quadratic term bx^2is not present in the cubic, the inflection point is on the y axis.
// TODO: We also know that a cubic function can have at most three x-intercepts.
// TODO: The constant d in the equation is the y -intercept of the graph. at  x = d
// TODO: Fifthly, we know that if the coefficient of x^3 is negative, then the cubic, notwithstanding any local extrema, 
// is falling with a negative slope from top left to bottom right. If the coefficient of x^3 is positive, 
// the curve is generally rising with positive slope from bottom left to top right.

class CubicFunction : public PolynomialFunction
{
private:
	double m_a = 0;
	double m_b = 0;
	double m_c = 0;
	double m_d = 0;

	std::pair<double, double> m_JustAAndDCubic = std::pair<double,double>(0,0);
	std::tuple<double, double, double> m_ACDForm;


	bool m_bIsInAAndDForm = false;
	bool m_bIsInACDForm = false;

	inline std::tuple<double, double, double> GetACD() const
	{
		return std::tuple<double, double, double>(m_a, m_c, m_d);
	}

public:

	CubicFunction() = default;

	explicit CubicFunction(const double& a, const double& b, const double& c, const double& d) 
		: m_a(a), m_b(b), m_c(c), m_d(d)
	{
		m_PolyFunctionType = PolynomialFunctionType::CUBIC;

		// TODO: might have to change this later?
		// Linear functions have all real number domain/ranges
		m_Domain = Domain::NegInfinityToPosInfinity;
		m_Range = Range::NegInfinityToPosInfinity;
		
		if (b == 0 && c == 0)
		{
			m_JustAAndDCubic = std::pair<double,double>(a, d);
			m_bIsInAAndDForm = true;
		}

		if (b == 0)
		{
			m_ACDForm = GetACD();
			m_bIsInACDForm = true;
		}

	}

	inline std::pair<double, double> GetAAndDCubicFuncForm() const { return m_JustAAndDCubic; }

	inline std::tuple<double, double, double, double> GetABCD() const
	{
		return std::tuple<double, double, double, double>(m_a, m_b, m_c, m_d);
	}



	bool GetIsFuncInAAndDForm() const { return m_bIsInAAndDForm; }
	bool GetIsFuncInACDForm() const { return m_bIsInACDForm; }


	std::tuple<double, double, double> GetACDForm() const { return m_ACDForm; }

};

// f(x)=ax+b
class LinearFunction : public PolynomialFunction
{
private:
	double m_a;
	double m_b;

	LineBehavior m_LineBehavior;


	float m_Slope;

	pair<double, double> m_YIntercept;
	pair<double, double> m_XIntercept;

	bool m_bIsBOnlyForm = false;

public:

	// if a > b etc.. TODO:

	// TODO: finish filling out the specific information gathering functions

	LinearFunction() = default;

	explicit LinearFunction(double a, double b) : m_a(a), m_b(b)
	{
		m_PolyFunctionType = PolynomialFunctionType::LINEAR;

		// TODO: might have to change this later?
		// Linear functions have all real number domain/ranges
		m_Domain = Domain::NegInfinityToPosInfinity;
		m_Range = Range::NegInfinityToPosInfinity;

		m_Slope = a;

		if (m_Slope == 0)
		{
			m_Degree = 0;
			m_bIsBOnlyForm = true;
		}
		else if (m_Slope != 0)
		{
			m_Degree = 1;
		}
		else
		{
			throw std::exception("Something went wrong with assigning the degree");
		}

		m_XIntercept.first = (-1 * b) / m_Slope;
		m_XIntercept.second = 0;

		m_YIntercept.first = 0;
		m_YIntercept.second = b;

		if (a > 0)
		{
			// Line rises as x increases
			m_LineBehavior = LineBehavior::Increasing;


		}
		else if (a < 0)
		{
			// line falls as x increases
			m_LineBehavior = LineBehavior::Decreasing;
		}
		else
		{
			// a == 0
			// horizontal line
			m_LineBehavior = LineBehavior::Horizontal;
		}

	}

	double operator()(double x) const { return ((m_a*x) + m_b); }

	double GetA() const { return m_a; }
	double GetB() const { return m_b; }

	// TODO: add functioanlity for verticle lines... 
	//ax + by = c,
	//	ax + by = c,
	//where a, ba, b are both not zero, to denote the standard form of a line.

	inline void PrintLinearFunctionInfo() const { std::cout << " f(x) = " << m_a << "x^2 + " << m_b << std::endl; }

	inline bool IsBOnlyForm() const { return m_bIsBOnlyForm; }


};

class QuadraticFunction;

std::vector<double> GetZerosQuadraticFormula(QuadraticFunction& QuadraticFunc);

std::vector<double> GetZerosQuadraticFormula(const double& a, const double& b, const double& c);

// quadratic function has the form f(x)=ax^2+bx+c,where a≠0
class QuadraticFunction : public PolynomialFunction
{
private:
	double m_a;
	double m_b;
	double m_c;

	ParabolaOpen m_ParabolaOpens;

	std::vector<double> m_RealNumberZeros;
	int m_AmountOfRealNumberZeros;
	//AmountOfRealNumberSoltuions m_RealSolutionAmount;

	std::vector<double> m_AllZeros;

	double m_SymmetryLine;

	double m_MaxValueAtXIsEqualTo;
	double m_MinValueAtXIsEqualTo;

	void AutomaticSetRealZeroVariables();

	bool m_bJustABForm = false;
	bool m_bJustACForm = false;

public:

	QuadraticFunction() = default;

	explicit QuadraticFunction(double a, double b, double c = 0) : m_a(a), m_b(b), m_c(c)
	{
		m_PolyFunctionType = PolynomialFunctionType::QUADRATIC;

		if (a == 0)
			throw std::exception("a cannot == 0 for quadratic func initalization");

		AutomaticSetRealZeroVariables();

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

		GetZerosQuadraticFormula(*this);

	}
	
	inline std::tuple<double, double, double> GetABC() const
	{
		return std::tuple<double, double, double>(m_a, m_b, m_c);
	}

	void SetAllZeroVec(std::vector<double> InVec) { m_AllZeros = InVec; }
	void SetRealNumberZeroVec(std::vector<double> InVec) { m_RealNumberZeros = InVec; }
	void SetLineOfSymmetry(double InNum) { m_SymmetryLine = InNum; }
	void SetTheMaxMinValue(double InNum);


	void PrintFunctionEndBehavior() const;
	void PrintParabolaOpensDirection() const;
	void PrintAllZeros() const;
	void PrintAllRealNumberZeros() const;
	void PrintNumberOfRealNumberSoltions() const;

	void PrintBasicFunctionInfo() const;

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
	
	std::vector<double> GetRealNumberZerosVec() const { return m_RealNumberZeros; }
	inline int GetAmountOfRealZeros() const { return m_AmountOfRealNumberZeros; }

};

template <typename T>
bool is_close_to_zero(T x)
{
	return std::abs(x) < std::numeric_limits<T>::epsilon();
}



inline bool isEven(int n);

// f(x)= a*sqrt(x - b) + c
class RootFunction : public PolynomialFunction
{
private:
	double m_n;

	double m_a;
	double m_b;
	double m_c;

	// TODO: Set domain and range in constructor


public:
	RootFunction() = default;

	explicit RootFunction(const double& n, const double& a,
		const double& b, const double& c) : m_n(n), m_a(a), m_b(b), m_c(c) 
	{

		// ROOT FUNCTIONS ARE NOT POLYNOMIAL FUNCTIONS AHHHHHHHHHHH
		m_PolyFunctionType = PolynomialFunctionType::ROOT;

		if (n == 0)
			throw std::exception("Pretty sure this is invalid algebra");

		if (n >= 1 && !isEven(n))
		{
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_bIsEvenFunction = false;

		}

		if (n >= 2 && isEven(n))
		{
			m_Domain = Domain::InclusiveZeroToPosInfinity;
		}
	}

	double operator()(const double x) const
	{
		double Power = 1.0 / m_n;
		double Base = x - m_b;
		double RootRes = std::pow(Base, Power);


		std::cout << "Power: " << Power << std::endl;
		std::cout << "Base: " << Base << std::endl;
		std::cout << "RootRes: " << RootRes << std::endl;
		std::cout << "A: " << m_a << std::endl;
		std::cout << "C: " << m_c << std::endl;
		return (m_a * (RootRes)) + m_c;
	}

	inline std::tuple<double, double, double, double> GetNABC() const
	{
		return std::tuple<double, double, double, double>(m_n, m_a, m_b, m_c);
	}


};

inline bool isEven(int n)
{
	if (n % 2 == 0)
		return true;
	else
		return false;
}

// The power function f(x)=x^n is an even function if n is even and n≠0,
// and it is an odd function if  n  is odd
class PowerFunction : PolynomialFunction
{

};

class TranscendentalFunction
{

};

// Formally, a function f is periodic if there exists a number p such that f(x + p) = f(x) 
// for all x. The smallest possible value of p is the period.
// The reciprocal of period is frequency.
enum class Period
{
	OnePI,
	TwoPI
};


// PI/180 == 180/PI == 1
// Convert 225 degrees to radians = (225 * pi) / 180 = 5pi/4
// Convert 5pi/3 to degrees = 5pi/3 * 180/pi == 900 / 3 == 300

// A trigonometric function relates the ratios of two sides of a right triangle. 
// They are sinx,cosx,tanx,cotx,secx,andcscx.

// going to define this function using the unit circle
/*
Let P=(x,y) be a point on the unit circle centered at the origin O.
Let θ be an angle with an initial side along the positive x-axis 
and a terminal side given by the line segment OP.
The trigonometric functions are then defined as

sin(theta) = y
cos(theta) = x
tan(theta) = y/x

*/


// Trig Identity Link below
// http://www.mathwords.com/t/trig_identities.htm

//  Acute Angle
// An angle that has measure less than 90°.
// SOHCAHTOA Defining only works on acute angles 

// The sine, cosine, secant, and cosecant functions have a period of 2π.
// Since the tangent and cotangent functions repeat on an interval of length π

class TrigonometricFunction : TranscendentalFunction
{
private:
	Period m_Period;
	
	double m_PeriodB;

	double m_HorizontalShifta;

	double m_VerticleStretchA;

	double m_VerticleShiftC;


	bool m_SinFunctionInUse = false;
	bool m_CosFunctionInUse = false;


	std::vector<std::pair<double, double>> m_SinFunctionPoints
	{
		std::pair<double,double>(0,0),
		std::pair<double,double>(M_PIOverTwo, 1),
		std::pair<double,double>(M_PI, 0),
		std::pair<double,double>(ThreeM_PIOverTwo, -1),
		std::pair<double,double>(TwoM_PI, 0)
	};

	std::vector<std::pair<double, double>> m_CosFunctionKeyPoints
	{
		std::pair<double,double>(0,1),
		std::pair<double,double>(M_PIOverTwo, 0),
		std::pair<double,double>(M_PI, -1),
		std::pair<double,double>(ThreeM_PIOverTwo, 0),
		std::pair<double,double>(TwoM_PI, 1)
	};

	double FixAngleBetweenZeroAndTwoRad(const double& InAngle);


public:

	explicit TrigonometricFunction(const double& A, const double& B, const double& a, const double& C)
	//	: m_VerticleStretchA(A), m_PeriodB(B), m_HorizontalShifta(a), m_VerticleShiftC(C)
	{
		// TODO: This is currently setup to automatically evaluate a sin function
		// Add option for cos/tan/create empty etc..


		// why when i used these two variables in the formula it didnt work???
		m_PeriodB = TwoM_PI / std::abs(B);
		m_VerticleStretchA = std::abs(A);
		

		m_HorizontalShifta = a;
		m_VerticleShiftC = C;

		m_SinFunctionInUse = true;

		//for (const auto& Point : m_SinFunctionPoints)
		//{
		//	std::cout << "x: " << Point.first << " f(x): " << Point.second << std::endl;
		//}

		for (int i = 0; i < m_SinFunctionPoints.size(); ++i)
		{
			// This didnt work the one below did?
			//m_SinFunctionPoints[i].second = ((std::sin((m_SinFunctionPoints[i].first - m_HorizontalShifta)
			//		* m_PeriodB) * m_VerticleStretchA) + m_VerticleShiftC);

			m_SinFunctionPoints[i].second = 
				((std::sin((m_SinFunctionPoints[i].first - a) * B) * A) + C);
		}

		PrintCurrentSinFunctionKeyPoints();
	}

	double SinOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);
	double CosOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);
	double TanOfUnitCircleAngle(const UnitCircle& InUnitCircle, const double& InAngle);

	inline void PrintCurrentSinFunctionKeyPoints() const
	{
		for (const auto& Point : m_SinFunctionPoints)
			std::cout << "x: " << Point.first << " f(x): " << Point.second << std::endl;
	}

	inline double EvaluateInverseSinComposition(const double& InSinNum)
	{
		UnitCircle TempUnitCircle;
		auto Map = TempUnitCircle.GetRadianSinCosMap();

		//std::cout << InSinNum << std::endl;

		// Inverse sin domain [1, 1]
		// Inverse sin range [-90,+90] in deg

		//const double NegPIOverTwo = M_PIOverTwo * (-1);	
		
		double MaxRange = ThreeM_PIOverTwo;
		
		double OutResult(0);
		double OutRange(0);

		for (auto Num : Map)
		{
			if (Num.second.second == InSinNum) 
			{
				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
				//	<< Num.second.second << std::endl;

				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				//	<< InSinNum << std::endl;

				OutResult = Num.first;
				break;
			}
		}
		bool FirstRange = (OutResult >= ThreeM_PIOverTwo && OutResult <= 2 * M_PI);
		bool SecondRange = (OutResult >= 0 && OutResult <= M_PIOverTwo);

		bool InValidRange = (FirstRange || SecondRange);

		if (InValidRange)
		{

			std::cout << OutResult << std::endl;
			std::cout << " or \n";
			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)

			return OutResult;
		}
		else
		{
			std::cout << OutResult << std::endl;

			return OutResult;
		}
	}

	inline double EvaluateInverseCosComposition(const double& InCosNum)
	{
		UnitCircle TempUnitCircle;
		auto Map = TempUnitCircle.GetRadianSinCosMap();

		//std::cout << InSinNum << std::endl;

		// Inverse cos domain [-1, 1]
		// Inverse cos range [0,PI] in rad

		double OutResult(0);
		double OutRange(0);

		for (auto Num : Map)
		{
			if (Num.second.first == InCosNum)
			{
				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) 
				//	<< Num.second.second << std::endl;

				//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				//	<< InSinNum << std::endl;

				OutResult = Num.first;
				break;
			}
		}
		bool FirstRange = (OutResult <= M_PI && OutResult >= 0);
		bool InValidRange = (FirstRange); // || SecondRange);

		if (InValidRange)
		{

			std::cout << OutResult << std::endl;
			std::cout << " or \n";
			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)

			return OutResult;
		}
		else
		{
			std::cout << OutResult << std::endl;

			return OutResult;
		}
	}

	inline double EvaluateInverseTanComposition(const double& InTanNum)
	{
		UnitCircle TempUnitCircle;
		auto Map = TempUnitCircle.GetRadianTanMap();


		// Inverse tan domain R
		// Inverse tan range [-pi/2,pi/2] in rad

		double OutResult(0);
		double OutRange(0);

		for (auto Num : Map)
		{
			if (Num.second == InTanNum)
			{
				OutResult = Num.first;
				break;
			}
		}
		bool Valid = (OutResult <= ThreeM_PIOverTwo && OutResult >= M_PIOverTwo);
		bool InValidRange = (Valid);

		if (InValidRange)
		{

			std::cout << OutResult << std::endl;
			std::cout << " or \n";
			std::cout << OutResult - (2 * M_PI) << std::endl; // (-pi/3)

			return OutResult;
		}
		else
		{
			std::cout << OutResult << std::endl;

			return OutResult;
		}
	}

};

inline double EvaluateCommonAngleSin(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, Point>  Map = TestCircle.GetRadianSinCosMap();

	double TestSinRes;

	for (const auto& Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestSinRes = Val.second.second;
		}
	}
	return TestSinRes;
}

inline double EvaluateCommonAngleCos(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, Point>  Map = TestCircle.GetRadianSinCosMap();

	double TestCosRes;

	for (const auto& Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestCosRes = Val.second.first;
		}
	}
	return TestCosRes;
}

inline double EvaluateCommonAngleTan(const double& InAngle)
{
	UnitCircle TestCircle;
	std::map<double, double>  Map = TestCircle.GetRadianTanMap();

	double TestTanRes(0);

	for (auto Val : Map)
	{
		if (Val.first == InAngle)
		{
			TestTanRes = Val.second;
			
		}
	}
	return TestTanRes;
}


// Exponent Rules
// http://dailprice.weebly.com/uploads/1/3/2/7/13279202/2545731_orig.png

//  An exponential function is a function of the form f(x)= a*b^x, where a > 0
// where the base 0 < b < 1 or b > 1
class ExponentialFunction : TranscendentalFunction
{
private:
	double m_a;
	double m_b;

	// If not growth its decay
	bool m_bIsExponentialGrowthFunc = true;

	Domain m_Domain;
	Range m_Range;

	EndBehavior m_AsXGoesPosEndBehavior;
	EndBehavior m_AsXGoesNegEndBehavior;

public:
	// Exponential functions have constant bases and variable exponents
	explicit ExponentialFunction(const double& a, const double& b)
		: m_a(a), m_b(b)
	{
		// a must be positive
		const bool AIsNotGreaterThanZero = (!(a > 0));

		if (AIsNotGreaterThanZero)
			throw std::exception("a needs to be greater than 0");

		// exponential growth or exponential decay?
		if (b > 1)
		{
			// growth
			m_bIsExponentialGrowthFunc = true;
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_Range = Range::ExclusiveZeroToPosInfinity;
			// Increasing on (neg inf to pos inf)
			// b^x -> pos inf AS x -> pos inf 
			// b^x -> 0 AS x -> neg inf
			m_AsXGoesPosEndBehavior = EndBehavior::AsXGoesToPosInfinityFOfXGoesToPosInfinity;
			m_AsXGoesNegEndBehavior = EndBehavior::AsXGoesToNegInfinityFOfXGoesToZero;
		}
		else if (b > 0 && b < 1)
		{
			// decay
			m_bIsExponentialGrowthFunc = false;
			m_Domain = Domain::NegInfinityToPosInfinity;
			m_Range = Range::ExclusiveZeroToPosInfinity;
			// decreasing on (neg inf to pos inf)
			// b^x -> 0 AS x -> pos inf 
			// b^x -> pos inf AS x -> neg inf
			m_AsXGoesPosEndBehavior = EndBehavior::AsXGoesToPosInfinityFOfXGoesToZero;
			m_AsXGoesNegEndBehavior = EndBehavior::AsXGoesToNegInfinityFOfXGoesToPosInfinity;
		}
		else if (b == 1)
		{
			throw std::exception("b cannot be == to 1");
		}
		else
		{
			// b < 0
			throw std::exception("b needs to be greater than 0");
		}
		
	}

};

// TODO: maybe move these 2 formula functions to a different file later.
//http://www.mathwords.com/c/continuously_compounded_interest.htm
// Continous Compound Intrest
// P = principal - Starting Amount 
// r = rate of intrest per year // Example: .032 which is 3.2% as a decimal
// t = number of years
// return A = final amount 
double FindCompoundIntrest(
	const double& P,
	const double& r,
	const double& t);

// http://www.softschools.com/formulas/math/compound_interest_formula/138/
// Compound Intrest Formula (Exponential Func)
// P = principal
// r = intrest rate as a decimal // Example: .032 which is 3.2% as a decimal
// n = number of times compounded per year
// t = number of years
// return A = the future value a particular investment will have.
double FindCompoundIntrest(
	const double& P,
	const double& r, 
	const double& n,
	const double& t);



// Since log and exponential are inverses the line below is true
// ln(e^x) = x and e^ln(x) = x

// NOTE:
/*
Inverse Hyperbolic Functions
From the graphs of the hyperbolic functions,
we see that all of them are one-to-one except coshx and sechx. 
If we restrict the domains of these two functions to the interval [0,∞),
then all the hyperbolic functions are one-to-one, 
and we can define the inverse hyperbolic functions.
*/

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
			throw std::exception("b value of log func has invalid value");

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

inline double COSH(const double& x)
{
	return (std::exp(x) + std::exp(-x)) / 2;
}

inline double SINH(const double& x)
{
	return (std::exp(x) - std::exp(-x)) / 2;
}

inline double TANH(const double& x)
{
	return SINH(x) / COSH(x);
}

inline double CSCH(const double& x)
{
	return 1.0 / SINH(x);
}
inline double SECH(const double& x)
{
	return 1.0 / COSH(x);
}

inline double COTH(const double& x)
{
	return COSH(x) / SINH(x);
}

// HYPERBLOIC Identitys link below
// https://image.slidesharecdn.com/lesson3derivativeofhyperbolicfunctions-150718183502-lva1-app6892/95/lesson-3-derivative-of-hyperbolic-functions-5-638.jpg?cb=1437244604
// sinh^-1 -- arcsinhx -- ln(x+sqrt(x^2+1))
inline double FindInverseHyperbolicSin(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (x + (std::sqrt(std::pow(x, 2) + 1)));
	double LastStep = std::log(InsideParathesis);

	return LastStep;

}

inline double FindInverseHyperbolicCos(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (x + (std::sqrt(std::pow(x, 2) - 1)));
	double LastStep = std::log(InsideParathesis);

	return LastStep;

}

inline double FindInverseHyperbolicTan(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = ((1 + x) / (1 - x));
	double NextStep = std::log(InsideParathesis);
	double LastStep = NextStep * 0.5;

	return LastStep;
}

inline double FindInverseHyperbolicCot(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = ((x + 1) / (x - 1));
	double NextStep = std::log(InsideParathesis);
	double LastStep = NextStep * 0.5;

	return LastStep;
}

inline double FindInverseHyperbolicSec(const double& x)
{
	// TODO: should i block the invalid input or what?
	// undefined for x <= 0, x > 1

	// setup the formula and return result
	double InsideParathesis =  (1 + (std::sqrt(1 - std::pow(x, 2)))) / x;
	double LastStep = std::log(InsideParathesis);

	return LastStep;
}

inline double FindInverseHyperbolicCsc(const double& x)
{
	// setup the formula and return result
	double InsideParathesis = (1 / x) + ((std::sqrt(1 + std::pow(x, 2))) / std::abs(x));
	double LastStep = std::log(InsideParathesis);

	return LastStep;
}

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

public:

	//RationalFunction() = delete;

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


	inline PolynomialFunctionType GetNumeratorFunctionType() const
	{
		return m_NumeratorFuncType;
	}

	inline PolynomialFunctionType GetDenominatorFunctionType() const
	{
		return m_DenominatorFuncType;
	}


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


// TODO: Make this function print to std::cout
// TODO: Connect this functionality to work with the other function classes
void PrintFunctionTransformationInfo();

// Notes for making the range domain finders...

/*

Set the denominator equal to zero, if it’s a fraction.
denominator of this function is (x - 1).
Set it equal to zero and solve for x: x – 1 = 0, x = 1.
Write the domain: The domain of this function cannot include 1, but includes all real numbers except 1; 
therefore, the domain is (-∞, 1) U (1, ∞).

Set the terms inside the radical to be greater than or equal to zero,
if there’s a root function.
For example: Identify the domain of the function f(x) = √(x + 3).
The terms within the radical are (x + 3).
Set them greater than or equal to zero: (x + 3) ≥ 0.
Solve for x: x ≥ -3.
The domain of this function includes all real numbers greater than or equal to -3;
therefore, the domain is [-3, ∞).

// Range of quadratic func

Confirm that you have a quadratic function
Find the x-value of the vertex of the function.
For example, find the range of 3x^2 + 6x -2.
Calculate x-coordinate of vertex: x = -b/2a = -6/(2*3) = -1
Calculate y-coordinate: y = 3x2 + 6x – 2 = 3(-1)2 + 6(-1) -2 = -5.
The vertex of this function is (-1, -5).
Determine the direction of the parabola by plugging in at least one more x-value.
Choose any other x-value and plug it into the function to calculate the corresponding
y-value. If the y-value is above the vertex, the parabola continues to +∞. 
If the y-value is below the vertex, the parabola continues to -∞.
Use the x-value -2: y = 3x2 + 6x – 2 = y = 3(-2)2 + 6(-2) – 2 = 12 -12 -2 = -2.
This yields the coordinate (-2, -2).
This coordinate tells you that the parabola continues above the vertex (-1, -5); therefore, the range encompasses all y-values above -5.
The range of this function is [-5, ∞)

Minor Note involving >= on sqrt functions
solve sqrt(x+3) + 1  Implies -- sqrt(x+3) = -1 
since sqrt(x+3) >= 0 for all x 
this equation has no solutions, and therefore f has no zeros.

// In other words, f(−x)=f(x). If a function f has this property,
we say f is an even function, which has symmetry about the y-axis

if a function f is symmetric about the origin,
then whenever the point (x,y) is on the graph,
the point (−x,−y) is also on the graph.
In other words, f(−x)=−f(x)
If f has this property, we say ff is an odd function,
which has symmetry about the origin.

Domain and range of an absolute value function
the absolute value function is defined for all real numbers,
the domain of this function is (−∞,∞).

for the function f(x) = 2|x - 3| + 4
|x - 3| >= 0 for all x 
2|x - 3| + 4 >= 4 
The range is {y|y>=4}



*/

class BasicFunction
{
private:
	std::function<double(const double&)> m_Function;
	std::string m_InputRange; // used for piecewise func
	int m_RangeVariable;

public:
	BasicFunction(std::function<double(const double&)> InFunc, const std::string& InInputRange,
		int RangeVariable)
		: m_Function(InFunc), m_InputRange(InInputRange), m_RangeVariable(RangeVariable)
	{

	}


	int GetRange() const { return m_RangeVariable; }
	std::string GetInputRange() const { return m_InputRange; }
	void SetInputRange(const std::string& InRange) { m_InputRange = InRange; }
	void SetRangeVariable(const int& Range) { m_RangeVariable = Range; }
	std::function<double(const double&)> GetFunction() const { return m_Function; }

	double RunFunctionNoRangeCheck(const double& x)
	{
		return this->GetFunction()(x);
	}
};

class PiecewiseFunction
{
private:
	BasicFunction m_FirstFunc;
	BasicFunction m_SecondFunc;

	int m_RangeVariable{ 0 };
	//std::string m_RangeString;

public:
	PiecewiseFunction(const BasicFunction& InFirstFunc, const BasicFunction& InSecondFunc)
		: m_FirstFunc(InFirstFunc), m_SecondFunc(InSecondFunc)
	{
		m_RangeVariable = InFirstFunc.GetRange();
	}

	double RunFunction(const double& x)
	{
		std::string FirstFuncInputRange = m_FirstFunc.GetInputRange();
		std::string SecondFuncInputRange = m_SecondFunc.GetInputRange();

		if (FirstFuncInputRange == ">")
		{
			if (x > m_RangeVariable)
			{
				return m_FirstFunc.GetFunction()(x);
			}
			else
			{
				return m_SecondFunc.GetFunction()(x);
			}

		}
		else if (FirstFuncInputRange == ">=")
		{
			if (x >= m_RangeVariable)
			{
				return m_FirstFunc.GetFunction()(x);
			}
			else
			{
				return m_SecondFunc.GetFunction()(x);
			}
		}
		else if (FirstFuncInputRange == "<")
		{
			if (x < m_RangeVariable)
			{
				return m_FirstFunc.GetFunction()(x);
			}
			else
			{
				return m_SecondFunc.GetFunction()(x);
			}
		}
		else if (FirstFuncInputRange == "<=")
		{
			if (x <= m_RangeVariable)
			{
				return m_FirstFunc.GetFunction()(x);
			}
			else
			{
				return m_SecondFunc.GetFunction()(x);
			}
		}
	}
};

class CalculusFunction
{
private:
	bool m_IsPiecewise{ false };


public:
	CalculusFunction();


	~CalculusFunction();
};

// TODO: fix input ranges to be accepted elsewhere
double SimpleFunction(const double& x);

double SimpleFunction2(const double& x);

double SimpleFunction3(const double& x);


// Start of calculus chapter 
// TODO: Sort/organize this stuff later

// Slope of a secant line formula
// through Points (a,f(a)) and (x,f(x))
inline double GetSlopeOfSecantLineTwoPoints(const Point& FirstPoint, const Point& SecondPoint)
{
	double a = FirstPoint.first;
	double FOfa = FirstPoint.second;

	double x = SecondPoint.first;
	double FOfx = SecondPoint.second;

	// use formula return result
	double Numerator = FOfx - FOfa;
	double Denominator = x - a;

	double SlopeOfSecantLine = Numerator / Denominator;

	return SlopeOfSecantLine;
}


// TODO: This functions variables need to be better understood
//  average velocity of an object over a time period to be the change in its position divided by the length of the time period.
inline double GetAverageVelocity(const double& a, const double& t,
	std::function<double(const double&)> s)
{
	// x == t
	/*
	Let s(t)s(t) be the position of an object moving along a coordinate axis at time t. 
	The average velocity of the object over a time interval [a,t] where a<t
	(or [t,a] if t<a)t<a) 
	is
	
	s(t) - s(a)
	/
	t - a

	*/
	if (a < t)
	{
		std::cout << "a < t using iterval [a, t]\n";

		// use formula return result
		double Numerator = s(t) - s(a);
		double Denominator = t - a;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}

	if (a > t)
	{

		std::cout << "t < a using iterval [t, a]\n";
		// use formula return result
		double Numerator = s(a) - s(t);
		double Denominator = a - t;

		double AverageVelocity = Numerator / Denominator;

		return AverageVelocity;
	}
}

// Interval [0,3] under estimate means
// evaluate the function for 0,1,2 and add the results
// underestimate = false means evaluate for 1 2 3
inline double GetAreaUnderCurve(const int& IntervalStart, const int& IntervalEnd,
	std::function<double(const double&)> func,
	const bool& bIsUnderEstimate = true)
{
	double OutResult{ 0.0 };

	if (bIsUnderEstimate)
	{
		for (int i = IntervalStart; i < IntervalEnd; i++)
		{
			OutResult = OutResult + func(i);
		}
	}

	// TODO: is this bottom part correct?
	if (bIsUnderEstimate == false)
	{
		for (int i = IntervalStart + 1; i <= IntervalEnd; ++i)
		{
			OutResult = OutResult + func(i);
		}
	}

	return OutResult;
}

// The limit of a constant is a constant 
// The limit of x as x approaches a is a

// TODO: maybe make a limit class file later?
class Limit
{
private:

	double m_a;
	std::function<double(const double&)> m_Function;

	double m_L;

	// triggers to true if a root function is detected in the numerator or denominator
	// This disables the solution finding by factor/canceling
	bool m_bTryConjSolution = false;

	inline RationalFunction SolveByConjugateMultiplication(const RootFunction& Numerator, const LinearFunction& Denominator)
	{
		auto NumeratorVars = Numerator.GetNABC();

		auto N = std::get<0>(NumeratorVars);

		auto A = std::get<1>(NumeratorVars);

		auto B = std::get<2>(NumeratorVars);

		auto C = std::get<3>(NumeratorVars);
		
		LinearFunction NumeratorRes;
		RootFunction DenominatorRes;

		bool bBInputIsPositive = (B > 0);
		

		LinearFunction TopRes;
		RootFunction BottomRes;

		// find conjugate
		if (C < 0)
		{
			// add C in conjugate
			//double Conjugate = (A*(std::pow(m_a - B, (1.0 / N)))) + C;
			if (bBInputIsPositive)
			{
				 TopRes = LinearFunction(1, B*(-1) + C*(-C));// == x + TopRes
				 BottomRes = RootFunction(N, A, B*(1), (C*(-1))); // Bottom Res = 
			}
			else
			{
				TopRes = LinearFunction(1, B*(-1) + C*(-C));// == x + TopRes
				BottomRes = RootFunction(N, A, B*(1), (C*(-1))); // Bottom Res = 
			}

			std::cout << TopRes.GetA() << std::endl;
			std::cout << Denominator.GetA() << std::endl;
			std::cout << TopRes.GetB() << std::endl;
			std::cout << Denominator.GetB() << std::endl;

			if (TopRes.GetA() == Denominator.GetA() && TopRes.GetB() == Denominator.GetB())
			{
				

				std::cout << "TEST INSIDE CANCEL FUNCS" << std::endl;
				// Pretend these canceled Out
				NumeratorRes = LinearFunction(0, 1);


				DenominatorRes = BottomRes;
			}
		}
		else
		{
			// subtract C

		}

		return RationalFunction(NumeratorRes, DenominatorRes);

	}

	// helper function to help evaluate a limit if its a linear function
	inline double EvaluateLinearFuncLimit(LinearFunction& InLinearFunc)
	{
		double Tempa = InLinearFunc.GetA();
		double Tempb = InLinearFunc.GetB();

		double FirstLimit = ApplyLimitConstantMultipleLaw(Tempa);
		double SecondLimit = ApplyBasicLimitRuleForConstants(Tempb);

		return FirstLimit + SecondLimit;
	}

	// TODO: This function really needs cleaned up.
	inline double EvaluateRationalFuncLimit(const RationalFunction& InRationalFunc)
	{
		double NumeratorRes = 0;
		double DenominatorRes = 0;

		// TODO: For fraction functions that contain square roots. Conjugate multiplication. SET IT UP
		if (InRationalFunc.GetNumeratorFunctionType() == PolynomialFunctionType::ROOT &&
			InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::LINEAR)
		{
			RootFunction Numerator = InRationalFunc.GetNumeratorRoot();

			NumeratorRes = Numerator(m_a);

			//std::cout << NumeratorRes << std::endl;

			LinearFunction Denominator = InRationalFunc.GetDenominatorLinear();

			DenominatorRes = Denominator(m_a);

			// root function detected, disable factor/cancel solution option
			m_bTryConjSolution = true;


			if ((NumeratorRes == 0 && DenominatorRes == 0) && m_bTryConjSolution == true)
			{

				std::cout << "Testttttttttttttt" << std::endl;
				RationalFunction NewFactoredFunction = SolveByConjugateMultiplication(Numerator, Denominator);
			

				//NumeratorRes = NewFactoredFunction.GetNumeratorLinear()(m_a);
				DenominatorRes = NewFactoredFunction.GetDenominatorRoot()(m_a);
				



				if (NewFactoredFunction.GetNumeratorLinear().IsBOnlyForm())
				{
					NumeratorRes = NewFactoredFunction.GetNumeratorLinear().GetB();
				}
				//std::cout << "Test" <<  NumeratorRes << std::endl;

				//LinearFunction NewFactoredNumerator;
				//RootFunction NewFactoredDenominator;
			}
		}

		if (InRationalFunc.GetNumeratorFunctionType() == PolynomialFunctionType::QUADRATIC &&
			InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
		{
			QuadraticFunction Numerator = InRationalFunc.GetNumeratorQuadratic();

			if (Numerator.IsABForm())
			{
				std::tuple<double, double> AB = Numerator.GetAB();

				double A = std::get<0>(AB);
				double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
				double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);
				double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

				double B = std::get<1>(AB);
				double FinalB = ApplyBasicLimitRuleForX() * B;
				NumeratorRes = FinalA + FinalB;
			}
			else
			{

				// a b c form
				// you now know the form break up the limits
				std::tuple<double, double, double> ABC = Numerator.GetABC();

				double A = std::get<0>(ABC);

				double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
				double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

				double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

				double B = std::get<1>(ABC);
				double FinalB = ApplyBasicLimitRuleForX() * B;

				double C = std::get<2>(ABC);
				double FinalC = ApplyBasicLimitRuleForConstants(C);

				NumeratorRes = FinalA + FinalB + FinalC;
			}


			//if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
			//{
			QuadraticFunction Denominator = InRationalFunc.GetDenominatorQuadratic();

			if (Denominator.IsACForm()) // TODO: Set this up next
			{
				std::tuple<double, double> AC = Denominator.GetAC();

				double A = std::get<0>(AC);
				double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
				double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);
				double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

				//std::cout << FinalA << std::endl;

				double C = std::get<1>(AC);
				double FinalC = ApplyBasicLimitRuleForConstants(C);
				//std::cout << FinalC << std::endl;

				DenominatorRes = FinalA + FinalC;
				//std::cout << DenominatorRes << std::endl;
			}
			else
			{
				// a b c form
				// you now know the form break up the limits
				std::tuple<double, double, double> ABC = Denominator.GetABC();

				double A = std::get<0>(ABC);

				double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
				double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

				double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

				double B = std::get<1>(ABC);
				double FinalB = ApplyBasicLimitRuleForX() * B;

				double C = std::get<2>(ABC);
				double FinalC = ApplyBasicLimitRuleForConstants(C);

				DenominatorRes = FinalA + FinalB + FinalC;
			}


			//}


			if ((NumeratorRes == 0 && DenominatorRes == 0) && m_bTryConjSolution == false)
			{
				LinearFunction NewFactoredNumerator;
				LinearFunction NewFactoredDenominator;

				std::vector<double> NumeratorZeros = Numerator.GetRealNumberZerosVec();
				std::vector<double> DenominatorZeros = Denominator.GetRealNumberZerosVec();

				std::vector<double> NewNumerator;

				for (int i = 0; i < NumeratorZeros.size(); ++i)
				{
					auto result1 = std::find(std::begin(DenominatorZeros), std::end(DenominatorZeros), NumeratorZeros[i]);

					if (result1 != std::end(DenominatorZeros))
					{
						std::cout << "Numerator contains: " << NumeratorZeros[i] << '\n';
					}
					else
					{
						std::cout << "Numerator does not contain: " << NumeratorZeros[i] << '\n';
						NewNumerator.push_back(NumeratorZeros[i]);
					}
				}

				for (int i = 0; i < NewNumerator.size(); ++i)
				{
					std::cout << NewNumerator[i] << " ";
					std::cout << std::endl;
				}

				//double One = std::abs(std::floor(NewNumerator[0]));
				//double Two = NewNumerator[0];

				//std::cout << "abs " << One << std::endl;
				//std::cout << Two << std::endl;

				if (/*std::abs*/(std::floor(NewNumerator[0])) == NewNumerator[0])
				{
					std::cout << "NewNumerator[0] is whole\n";


					if (NewNumerator[0] == 0.0)
					{
						NewFactoredNumerator = LinearFunction(1, 0);
					}
					else
					{
						if (NewNumerator[0] < 0.0)
						{
							NewFactoredNumerator = LinearFunction(1, NewNumerator[0] * (-1));
						}
						else
						{
							NewFactoredNumerator = LinearFunction(1, NewNumerator[0]);
						}
					}
				}
				else
				{
					std::cout << "NewNumerator[0] is not whole\n";

					std::pair<double, double> NumeratorFract = OutputDecimalAsFract(NewNumerator[0]);

					std::cout << NumeratorFract.first << "/" << NumeratorFract.second << std::endl;

					if (NumeratorFract.first < 0)
					{
						NewFactoredNumerator = LinearFunction(NumeratorFract.second, NumeratorFract.first);
						//std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
					}
				}



				std::vector<double> NewDenominator;

				for (int i = 0; i < NumeratorZeros.size(); ++i)
				{
					auto result1 = std::find(std::begin(NumeratorZeros), std::end(NumeratorZeros), DenominatorZeros[i]);

					if (result1 != std::end(NumeratorZeros))
					{
						std::cout << "Numerator contains: " << DenominatorZeros[i] << '\n';
					}
					else
					{
						std::cout << "Numerator does not contain: " << DenominatorZeros[i] << '\n';
						NewDenominator.push_back(DenominatorZeros[i]);
					}
				}

				for (int i = 0; i < NewDenominator.size(); ++i)
				{
					std::cout << NewDenominator[i] << " ";
					std::cout << std::endl;
				}

				if (std::abs(std::floor(NewDenominator[0])) == NewDenominator[0])
				{
					std::cout << "NewDenominator[0] is whole\n";

					if (NewDenominator[0] == 0.0)
					{
						NewFactoredDenominator = LinearFunction(1, 0);
					}
					else
					{
						if (NewDenominator[0] < 0.0)
						{
							NewFactoredDenominator = LinearFunction(1, NewDenominator[0]);
						}
						else
						{
							NewFactoredDenominator = LinearFunction(1, NewDenominator[0] * (-1));
						}
					}

				}
				else
				{
					std::cout << "NewDenominator[0] is not whole\n";

					std::pair<double, double> DenominatorFract = OutputDecimalAsFract(NewDenominator[0]);
					if (DenominatorFract.first < 0)
					{
						std::cout << DenominatorFract.first << "/" << DenominatorFract.second << std::endl;
						NewFactoredDenominator = LinearFunction(DenominatorFract.second, DenominatorFract.first*(-1));

					}
				}

				NewFactoredNumerator.PrintLinearFunctionInfo();
				NewFactoredDenominator.PrintLinearFunctionInfo();

				double FinalFactoredNumerator = EvaluateLinearFuncLimit(NewFactoredNumerator);
				double FinalFactoredDenominator = EvaluateLinearFuncLimit(NewFactoredDenominator);

				std::cout << "FF Num: " << FinalFactoredNumerator << std::endl;
				std::cout << "FF Den: " << FinalFactoredDenominator << std::endl;

				return FinalFactoredNumerator / FinalFactoredDenominator;

			}
		}


		if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::CUBIC)
		{
			CubicFunction Denominator = InRationalFunc.GetDenominatorCubic();

			if (Denominator.GetIsFuncInAAndDForm())
			{

				std::pair<double, double> AAndD = Denominator.GetAAndDCubicFuncForm();

				double A = AAndD.first;
				
				double TempAPowerLawCubic = ApplyLimitRuleForPowers(3);
				double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);
				
				double FinalA = TempAConstMultipleLaw*TempAPowerLawCubic;
				
				double D = AAndD.second;
				double TempDConstRule = ApplyBasicLimitRuleForConstants(D);

				double FinalD = TempDConstRule;

				DenominatorRes = FinalA + FinalD;
			}
		}
		
		//if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::LINEAR)
		//{
		//	LinearFunction Denominator = InRationalFunc.GetDenominatorLinear();

		//	DenominatorRes = EvaluateLinearFuncLimit(Denominator);
		//}

		//if (InRationalFunc.GetDenominatorFunctionType() == PolynomialFunctionType::QUADRATIC)
		//{
		//	QuadraticFunction Denominator = InRationalFunc.GetDenominatorQuadratic();

		//	// a b c form
		//	// you now know the form break up the limits
		//	std::tuple<double, double, double> ABC = Denominator.GetABC();

		//	double A = std::get<0>(ABC);

		//	double TempAPowerLawQuadratic = ApplyLimitRuleForPowers(2);
		//	double TempAConstMultipleLaw = ApplyBasicLimitRuleForConstants(A);

		//	double FinalA = TempAConstMultipleLaw*TempAPowerLawQuadratic;

		//	double B = std::get<1>(ABC);
		//	double FinalB = ApplyBasicLimitRuleForX() * B;

		//	double C = std::get<2>(ABC);
		//	double FinalC = ApplyBasicLimitRuleForConstants(C);

		//	NumeratorRes = FinalA + FinalB + FinalC;
		//}

		// TODO: remove debug code
		std::cout << "Numerator:\t " << NumeratorRes << std::endl;
		std::cout << "Denominator:\t " << DenominatorRes << std::endl;


		//if (NumeratorRes == 0 && DenominatorRes == 0)
		//{

		//
		//	std::pair<double, double> Fract = OutputDecimalAsFract(NumeratorRes / DenominatorRes);

		//	std::cout << Fract.first << " / " << Fract.second << std::endl;

		//}


		return NumeratorRes / DenominatorRes;
	}

	inline double EvaluateCubicFunctionLimit(const CubicFunction& InCubicFunc)
	{
		if (InCubicFunc.GetIsFuncInACDForm())
		{
			std::tuple<double, double, double> ACDVariables = InCubicFunc.GetACDForm();

			double A = std::get<0>(ACDVariables);
			double C = std::get<1>(ACDVariables);
			double D = std::get<2>(ACDVariables);

			const double x = ApplyBasicLimitRuleForX();
			
			double PowerLawA = ApplyLimitRuleForPowers(3);
			double FinalA = PowerLawA * A;

			double FinalC = C*x;

			double FinalD = ApplyBasicLimitRuleForConstants(D);

			return FinalA + FinalC + FinalD;

		}

		return 0.0;

	}

	// TODO: add evaluations for other limit functions // logaritm/rationals etc.

	// This functions assumes you put into it the multiple of x of a linear function for the moment
	inline double ApplyLimitConstantMultipleLaw(const double& InA)
	{
		// m_a assigned during the constructor
		// This is moved outside of the limit here
		double c = InA;
		
		// Now inside of the limit all we have is an x
		// taking the limit of just an x is a!
		double LimitOfX = ApplyBasicLimitRuleForX();

		return c*LimitOfX;
	}

	inline double ApplyLimitRuleForPowers(const double& InN)
	{
		double XLimit = ApplyBasicLimitRuleForX();
		return std::pow(XLimit, InN);
	}

	inline double ApplyBasicLimitRuleForConstants(const double& InConstant)
	{
		double LimitOfConstant = InConstant;
		return LimitOfConstant;
	}

	inline double ApplyBasicLimitRuleForX()
	{
		return m_a;
	}



public:

	explicit Limit(std::function<double(const double&)> InFunc, const double& a)
		: m_Function(InFunc), m_a(a)
	{
		
		// automatically run the limit on construction
		m_L = operator()(a);
	}

	explicit Limit(LinearFunction& InLinearFunc, const double& a)
		: m_a(a)
	{
		m_L = EvaluateLinearFuncLimit(InLinearFunc);
		
		// TODO: remove debug code
		DisplayLimitResult();
	}

	explicit Limit(const CubicFunction& InCubicFunc, const double& a)
		: m_a(a)
	{
		m_L = EvaluateCubicFunctionLimit(InCubicFunc);

		// TODO: remove debug code
		DisplayLimitResult();
	}

	explicit Limit(const RationalFunction& InRationalFunc, const double& a)
		: m_a(a)
	{

		m_L = EvaluateRationalFuncLimit(InRationalFunc);

		// TODO: remove debug code
		DisplayLimitResult();
	}
	

	// TODO: Clean this function up by adding helper functions
	double operator()(double x) const 
	{

		std::vector<std::pair<double, double>> PosDirVec;
		// pos direction
		std::pair<double, double> PosPairFirst;
		PosPairFirst.first = m_a + 0.1;
		PosPairFirst.second = m_Function(m_a + 0.1);
		
		std::pair<double, double> PosPairSecond;
		PosPairSecond.first = m_a + 0.01;
		PosPairSecond.second = m_Function(m_a + 0.01);

		std::pair<double, double> PosPairThird;
		PosPairThird.first = m_a + 0.001;
		PosPairThird.second = m_Function(m_a + 0.001);
	
		std::pair<double, double> PosPairFourth;
		PosPairFourth.first = m_a + 0.0001;
		PosPairFourth.second = m_Function(m_a + 0.0001);
	
		PosDirVec.push_back(PosPairFirst);
		PosDirVec.push_back(PosPairSecond);
		PosDirVec.push_back(PosPairThird);
		PosDirVec.push_back(PosPairFourth);

		std::vector<std::pair<double, double>> NegDirVec;
		// neg direction
		std::pair<double, double> NegPairFirst;
		NegPairFirst.first = m_a - 0.1;
		NegPairFirst.second = m_Function(m_a - 0.1);

		std::pair<double, double> NegPairSecond;
		NegPairSecond.first = m_a - 0.01;
		NegPairSecond.second = m_Function(m_a - 0.01);

		std::pair<double, double> NegPairThird;
		NegPairThird.first = m_a - 0.001;
		NegPairThird.second = m_Function(m_a - 0.001);

		std::pair<double, double> NegPairFourth;
		NegPairFourth.first = m_a - 0.0001;
		NegPairFourth.second = m_Function(m_a - 0.0001);

		NegDirVec.push_back(NegPairFirst);
		NegDirVec.push_back(NegPairSecond);
		NegDirVec.push_back(NegPairThird);
		NegDirVec.push_back(NegPairFourth);

		std::cout << "Evaluating Limit: Please Wait...\n";

		for (auto & num : PosDirVec)
		{

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << num.first << " " << num.second << std::endl;
		}

		std::cout << std::endl;

		for (auto & num : NegDirVec)
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				 << std::setw(7) << num.first << " " << num.second << std::endl;
		}

		double TopRes =  PosDirVec[3].second;
		double BottomRes = NegDirVec[3].second;
		/*std::cout << TopRes << std::endl;
		std::cout << BottomRes << std::endl;*/
		
		
		std::string TopResStr = std::to_string(TopRes);
		// TODO: remove debug code later
		std::cout << TopResStr << std::endl;
	
		std::string BottomResStr = std::to_string(BottomRes);
		// TODO: remove debug code later
		std::cout << BottomResStr << std::endl;

		double LocalPosRes = std::stod(TopResStr) * 100;

		double LocalNegRes = std::stod(BottomResStr) * 100;

		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;
		


		std::cout << std::endl;

		// Are these limits infinite?
		bool bIsPosDirPosInfinity = false;
		bool bIsPosDirNegInfinity = false;
		bool bIsNegDirNegInfinity = false;
		bool bIsNegDirPosInfinity = false;

		if (PosDirVec[1].second > PosDirVec[0].second)
		{
			if (PosDirVec[2].second > PosDirVec[1].second)
			{
				bIsPosDirPosInfinity = true;
			}
		}

		if (NegDirVec[1].second < NegDirVec[0].second)
		{
			if (NegDirVec[2].second < NegDirVec[1].second)
			{
				bIsNegDirNegInfinity = true;
			}
		}

		if (PosDirVec[1].second < PosDirVec[0].second)
		{
			if (PosDirVec[2].second < PosDirVec[1].second)
			{
				bIsPosDirNegInfinity = true;
			}
		}
		if (NegDirVec[1].second > NegDirVec[0].second)
		{
			if (NegDirVec[2].second > NegDirVec[1].second)
			{
				bIsNegDirPosInfinity = true;
			}
		}

		if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
		{
			std::cout << "We have a working limit result! 1\n";
			

			if (bIsPosDirPosInfinity)
			{
				if (bIsNegDirPosInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << INFINITY << std::endl;
					return INFINITY;
				}
			}

			if (bIsPosDirNegInfinity)
			{
				if (bIsNegDirNegInfinity)
				{
					std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
					std::cout << NEGINFINITY << std::endl;
					return NEGINFINITY;
				}
			}

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";

			//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			//	// return either
			//return std::ceil(LocalPosRes) / 100;
		}

		// if you have two results that are returning negatives
		// you need to do some sort of swap with floor / ceil
		if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
		{
			std::cout << "We have a working limit result! 2\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;

			std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";

			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;


			// return either
			return std::floor(LocalPosRes) / 100;
		}


		std::cout << "Limit: DNE (does not exist): \n\n";

		std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
		if (bIsPosDirPosInfinity)
		{
			std::cout << INFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
		}

		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
		if (bIsNegDirNegInfinity)
		{
			std::cout << NEGINFINITY << std::endl;
		}
		else
		{
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
		}

		std::cout << std::endl;
		
		return 0;
	}

	//void TakeOneSidedLimits()


	//inline void TakeLimit(double x)
	//{
	//	m_L = operator()(x);
	//}

	inline void DisplayLimitResult()
	{
		std::cout << "Limit of f(x) as x --> " << m_a << " is " << m_L << std::endl;
	}
	

	
};



inline void EvaluateLimit(const Limit& InLimit)
{

}

