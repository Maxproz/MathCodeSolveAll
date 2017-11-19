#pragma once

#ifndef TRIGFUNC
#define TRIGFUNC

#include <utility>
#include <vector>
#include <iomanip>

#include "Conics.h" // For TrigFunctions enum

typedef double Radians;
typedef double Degrees;

constexpr double M_PI = 3.14159;
constexpr auto TwoM_PI = M_PI * 2;

enum class DoubleAngleReturnType
{
	SIN,
	COS,
	TAN
};

enum class Quadrants
{
	One,
	Two,
	Three,
	Four
};

enum class HalfAngleFormula
{
	Sin,
	Cos,
	Tan
};

enum class AngleType
{
	Radians,
	Degrees
};

enum class GivenSides
{
	AB,
	AC,
	BC
};

enum class GivenAngle
{
	Alpha,
	Beta,
	Gamma
};

enum class GivenSide
{
	A,
	B,
	C
};

double GetReferenceAngle(const Quadrants& Quad, const AngleType& Type,
	const double InAngle);

Degrees RadianToDegrees(const Radians InRad);

Radians DegreesToRadians(const Degrees InDeg);

// the division problem s/r -- how many r's are there in s?
Radians RadianMeasureOfArcInUnitCircle(const double ArcLength, const double Radius);

// circumference of circle
// C = 2(pi)(r) = pi(d)
double FindCircumferenceOfCircle(const double& Radius);
// area of circle
// A = (pi)(radius^2)
double FindAreaOfCircle(const double& Radius);

// Area of triangle ( these return square units )
// A = (1/2)(b)(h) or (1/2)(b)(c)(sinA) or (1/2)(a)(c)(sinB) or (1/2)(a)(b)(sinC)
double FindAreaOfTriangle(const double& Base, const double& Height);
double FindAreaOfTriangleAlpha(const double& B, const double& C, const double& AlphaAngle);
double FindAreaOfTriangleBeta(const double& A, const double& C, const double& BetaAngle);
double FindAreaOfTriangleGamma(const double& A, const double& B, const double& GammaAngle);

// Pythagorean Theorm
// a squared + bsquared = c squared
// user supplies 2 sides with a value and one side that they set to 0
// return true if success
bool FindMissingSidePythagoreanTheorm(double& base, double& height, double& hypotenuse);

// distance formula
// input is 4 points (2 cord pairs) 
// output is a double as the distance between them
double FindDistanceBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2);
// midpoint formula
// input is 4 points (2 cord pairs) 
// output is a double pair as the midpoint cord's between them
std::pair<double, double>
FindMidpointBetweenACordPair(const double& x1, const double& y1, const double& x2, const double& y2);

double FindSlopeWithTwoPoints(const double& x1, const double& y1, const double& x2, const double& y2);

void OutputSignsOfTrigFunctionQuadrants(const Quadrants& InQuad);

void OutputTrigEquivelantEquations();

// https://cnx.org/contents/E6wQevFf@6.38:aIPS8_HQ@8/Inverse-Trigonometric-Function
// Evaluating the Composition of a Sine with an Inverse Cosine
// We can use the Pythagorean identity to do this.
double CompositionofaSineWithAnInverseCosine(double InCos);

// TODO: This function needs to be set up to detect which quadrant the result is in
// Evaluating the Composition of a Sine with an Inverse Sine
// We can use the Pythagorean identity to do this.
double CompositionofaCosineWithAnInverseSine(double InSine);


// Reciprocal Identities

double GetReciprocalIdentity(const double Angle);

// SUM AND DIFFERENCE FORMULAS FOR COSINE

// Cosine of the Difference of Two Angles
// Assumes the user entered the angles in radians
double GetCosineOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle);

// Cosine of the Sum of Two Angles
// Assumes the user entered the angles in radians
double GetCosineOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle);


double GetProductOfSumCosACosB(const double& AlphaAngle, const double& BetaAngle);

double GetProductOfSumSinACosB(const double& AlphaAngle, const double& BetaAngle);

double GetProductOfSumSinASinB(const double& AlphaAngle, const double& BetaAngle);

double GetProductOfSumCosASinB(const double& AlphaAngle, const double& BetaAngle);

double SumToProductSinAddSin(const double& AlphaAngle, const double& BetaAngle);

double SumToProductCosSubtractCos(const double& AlphaAngle, const double& BetaAngle);

double SumToProductCosSubtractCosDegrees(const double& AlphaAngle, const double& BetaAngle, const bool& ReturnRadians = true);


// SUM AND DIFFERENCE FORMULAS FOR SINE

// SINE of the Difference of Two Angles
// Assumes the user entered the angles in radians
double GetSINEOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle);

// SINE of the Sum of Two Angles
// Assumes the user entered the angles in radians
double GetSINEOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle);

// https://cnx.org/contents/E6wQevFf@6.38:LKVQKWmT@5/Sum-and-Difference-Identities
// The sum and difference formulas for tangent are:
// tan(alpha + beta) = tan(alpha) + tan(beta) / 1 - tan(alpha)*tan(beta)
double GetTanOfTheSumOfTwoAngles(double AlphaAngle, double BetaAngle);

// tan(alpha - beta) = tan(alpha) - tan(beta) / 1 + tan(alpha)*tan(beta)
double GetTanOfTheDifferenceOfTwoAngles(double AlphaAngle, double BetaAngle);

// Finding Multiple Sums and Differences of Angles
/*
Given sin α = 3/5, 0<α<π2, cos β = −5/13, π<β<3π2,  find
sin(α + β)sin(α + β)
cos(α + β)cos(α + β)
tan(α + β)tan(α + β)
tan(α−β)
*/
void FindMultipleSumsAndDifferencesOfAngles(std::pair<double, double> InSineOfAlphaAngle,
	std::pair<double, double> InCosineOfBetaAngle);

// 68.622153857191046346385474553595
// 64.031242374328486864882176746218

// sinθ = 1/cscθ
// cosθ = 1/secθ
// tanθ = 1/cotθ
// csc θ = 1/sin θ
// sec θ = 1/cos θ
// cot  θ = 1/tan θ

double InSinOutCosCofunction(double InSinAngle);

double InTanOutCotCofunction(double InTanAngle);

double InSecOutCscCofunction(double InSecAngle);

double InCosOutSinCofunction(double InCosAngle);

double InCotOutTanCofunction(double InCotAngle);

double InCscOutSecCofunction(double InCscAngle);


// use The double-angle formulas given tan of an angle
double InTanOfAnAngleOutDoubleAngleSin(std::pair<double, double> InTanOfAngle,
	const DoubleAngleReturnType& TypeToUse);


void FindAltitudeOfObject(const double& AngleA, const double& AngleB,
	const double& DistanceBetweenC);

double LawOfSinesFindAngleAlpha(const double& SideA, const double& BetaAngle, const double& SideB);

double LawOfSinesFindAngleBeta(const double& SideA, const double& AlphaAngle, const double& SideB);

double LawOfSinesFindAngleGamma(const double& SideA, const double& AlphaAngle, const double& SideC);

void ConvertPolarToRectangle(const double& InRadius, const double InAngle, double& OutX, double& OutY);

void PrintCordinates(const double& x, const double& y);

void PrintCordinatesSTDPair(const std::pair<double, double>& Cords);

// Convert the rectangular coordinates (3, 3) to polar coordinates.
std::pair<double, double> ConvertToPolarFromRectangular(const double& x, const double &y);

template <typename T>
inline void PrintAllPolarPoints(const T& CONT)
{
	for (const auto& pair : CONT)
	{
		std::cout << std::fixed << "(" << pair.first << ", " << pair.second << ")" << std::endl;
	}
}



struct PairSort
{
	template <typename T>
	bool operator()(T a, T b) const
	{
		return a.first > b.first;
	}
};


template <typename T>
inline void PrintMaxPairForPolarContainer(const T& CONT)
{
	PairSort sort;

	auto cont = CONT;
	std::sort(cont.begin(), cont.end(), sort);

	std::cout << std::fixed << "(" << cont[0].first << ", " << cont[0].second << ")" << std::endl;
}


template <typename T>
inline T GetPairsOfTheZerosForPolarContainer(T& CONT)
{
	T OutCONT;

	for (T::iterator it = CONT.begin(); it != CONT.end(); ++it)
	{
		if ((*it).first <= 0.1 && (*it).first >= -0.1)
		{
			OutCONT.push_back(*it);
		}
	}

	return OutCONT;

}

// Recall that, to find the zeros of polynomial functions,
// we set the equation equal to zero and then solve for x.
// We use the same process for polar equations. Set r=0, and solve for θ.

std::vector<std::pair<double, double>> ReturnMaxMinPointsSin(const double& PoleDist);

std::vector<std::pair<double, double>> ReturnMaxMinPointsCos(const double& PoleDist,
	const double& PlusMinusVar = 0);

class PolarFuncSymmetryResult
{
public:

	void PrintTestResultsNewLines() const
	{
		if (GetFirstTestResult())
		{
			std::cout << "Symmetry around pi/2 Test Passes" << std::endl;
		}
		if (GetSecondTestResult())
		{
			std::cout << "Polar Axis Symmetry Test Passes" << std::endl;
		}
		if (GetThirdTestResult())
		{
			std::cout << "Respect to pole Symmetry Test Passes" << std::endl;
		}
	}

	void PassFirstTest()
	{
		m_PI_2Symmetry = true;
	}
	void PassSecondTest()
	{
		m_PolarAxisSymmetry = true;
	}
	void PassThirdTest()
	{
		m_RespectToPoleSymmetry = true;
	}

	bool GetFirstTestResult() const { return m_PI_2Symmetry; }
	bool GetSecondTestResult() const { return m_PolarAxisSymmetry; }
	bool GetThirdTestResult() const { return m_RespectToPoleSymmetry; }


private:
	bool m_PI_2Symmetry = false;
	bool m_PolarAxisSymmetry = false;
	bool m_RespectToPoleSymmetry = false;
};

PolarFuncSymmetryResult PolarSymmentryTest(//const double& r,
	const double& Theta,
	const TrigFunctions& Formula,
	const double& FuncMultiplier,
	const double& PlusMinusVar = 0);

// Used below in the following class
struct SecondAnswerTriangleDataLawOfSinesSSA
{
	typedef double Side;
	typedef double Angle;

	Side m_A;
	Side m_B;
	Side m_C;

	Angle m_Alpha;
	Angle m_Beta;
	Angle m_Gamma;

};

class LawOfSinesInfoSSA
{
public:
	typedef double Side;
	typedef double Angle;

	LawOfSinesInfoSSA(const Side& FirstSide, const Side& SecondSide,
		const GivenSides& InSides, const Angle& AngleAmount, const GivenAngle& InAngle)
	{

		switch (InSides)
		{
			case GivenSides::AB:
			{
				m_A = FirstSide;
				m_B = SecondSide;
				m_C = 0.0;

				SecondaryData.m_A = FirstSide;
				SecondaryData.m_B = SecondSide;
				SecondaryData.m_C = 0.0;


				if (m_B > m_A)
				{
					//m_IsObtuse = true;

				}
				else
				{
					m_IsObtuse = true;
					//throw std::exception("Error: You must enter an obtuse triangle to use this functionality");
				}

				break;
			}
			case GivenSides::AC:
			{
				m_A = FirstSide;
				m_C = SecondSide;
				m_B = 0.0;

				break;
			}
			case GivenSides::BC:
			{
				m_A = 0.0;
				m_B = FirstSide;
				m_C = SecondSide;

				SecondaryData.m_A = 0.0;
				SecondaryData.m_B = FirstSide;
				SecondaryData.m_C = SecondSide;

				if (m_C > m_B)
				{
					m_IsObtuse = true;
				}
				else
				{
					throw std::exception("Error: You must enter an obtuse triangle to use this functionality");
				}

				break;
			}
		}

		switch (InAngle)
		{
			case GivenAngle::Alpha:
			{
				m_Alpha = AngleAmount;
				SecondaryData.m_Alpha = AngleAmount;

				break;
			}
			case GivenAngle::Beta:
			{
				m_Beta = AngleAmount;
				break;
			}
			case GivenAngle::Gamma:
			{
				m_Gamma = AngleAmount;
				SecondaryData.m_Gamma = AngleAmount;
				break;
			}
		}



		CalculateMissingAngles(InAngle);
		CalculateMissingSide(InSides);

	}

	inline void CalculateMissingAngles(const GivenAngle& InAngle)
	{

		switch (InAngle)
		{
			case GivenAngle::Alpha:
			{
				double LHS = (m_B * std::sin(DegreesToRadians(m_Alpha))) / m_A;


				m_Beta = std::asin(LHS);
				m_Beta = RadianToDegrees(m_Beta);

				SecondaryData.m_Beta = m_Beta;


				if (m_IsObtuse == true)
				{
					m_Beta = 180.0 - m_Beta;
					m_Gamma = 180.0 - (m_Alpha)-(m_Beta);
				}

				SecondaryData.m_Gamma = 180.0 - SecondaryData.m_Alpha - SecondaryData.m_Beta;

				break;
			}
			case GivenAngle::Beta:
			{

				//m_Beta = AngleAmount;


				break;

			}
			case GivenAngle::Gamma:
			{

				double LHS = (m_B * std::sin(DegreesToRadians(m_Gamma))) / m_C;


				m_Beta = std::asin(LHS);
				m_Beta = RadianToDegrees(m_Beta);

				SecondaryData.m_Beta = m_Beta;


				if (m_IsObtuse == true)
				{
					m_Beta = 180.0 - m_Beta;
					m_Alpha = 180.0 - (m_Gamma)-(m_Beta);
				}

				SecondaryData.m_Alpha = 180.0 - SecondaryData.m_Gamma - SecondaryData.m_Beta;



				break;
			}
		}

	}

	inline void CalculateMissingSide(const GivenSides& InSides)
	{
		switch (InSides)
		{
			case GivenSides::AB:
			{

				m_C = (m_A * std::sin(DegreesToRadians(m_Gamma))) /
					std::sin(DegreesToRadians(m_Alpha));

				SecondaryData.m_C = (SecondaryData.m_A * std::sin(DegreesToRadians(SecondaryData.m_Gamma)) /
					std::sin(DegreesToRadians(SecondaryData.m_Alpha)));

				/*
				if (m_IsObtuse)
				{
				m_C =
				}*/


				break;
			}
			case GivenSides::AC:
			{
				//m_A = FirstSide;
				//	m_C = SecondSide;
				//m_B = 0.0;

				break;
			}
			case GivenSides::BC:
			{

				m_A = (m_C * std::sin(DegreesToRadians(m_Alpha))) /
					std::sin(DegreesToRadians(m_Gamma));

				SecondaryData.m_A = (SecondaryData.m_C * std::sin(DegreesToRadians(SecondaryData.m_Alpha)) /
					std::sin(DegreesToRadians(SecondaryData.m_Gamma)));

				//m_A = 0.0;
				//m_B = FirstSide;
				//m_C = SecondSide;

				break;
			}
		}
	}

	inline double GetAlphaAngle() const { return m_Alpha; }
	inline double GetGammaAngle() const { return m_Gamma; }
	inline double GetBetaAngle() const { return m_Beta; }

	inline double GetSideA() const { return m_A; }
	inline double GetSideB() const { return m_B; }
	inline double GetSideC() const { return m_C; }

	inline void PrintSideA() const
	{
		std::cout << std::setprecision(2) << "Side A: " << GetSideA() << std::endl;
		std::cout << std::fixed;
	}

	inline void PrintSideB() const
	{
		std::cout << std::setprecision(2) << "Side B: " << GetSideB() << std::endl;
		std::cout << std::fixed;
	}

	inline void PrintSideC() const
	{
		std::cout << std::setprecision(2) << "Side C: " << GetSideC() << std::endl;
		std::cout << std::fixed;
	}

	inline void PrintAllSidesNewLines() const
	{
		PrintSideA();
		PrintSideB();
		PrintSideC();
	}

	inline void PrintAngleAlpha() const
	{
		std::cout << "Alpha: " << GetAlphaAngle() << " Degrees" << std::endl;
	}

	inline void PrintAngleBeta() const
	{
		std::cout << "Beta: " << GetBetaAngle() << " Degrees" << std::endl;
	}

	inline void PrintAngleGamma() const
	{
		std::cout << "Gamma: " << GetGammaAngle() << " Degrees" << std::endl;
	}

	inline void PrintAllAnglesNewLines() const
	{
		PrintAngleAlpha();
		PrintAngleBeta();
		PrintAngleGamma();
	}

	inline void PrintAllSidesAndAngles() const
	{
		std::cout << "If there is any negative angles in these triangles.\nThen that is a failed triangle (ignore the results for it) " << std::endl << std::endl;

		PrintAllSidesNewLines();
		PrintAllAnglesNewLines();
		PrintOtherAnswerTriangle();
	}

	void PrintOtherAnswerTriangle() const
	{
		if (m_IsObtuse)
		{
			// next
			std::cout << std::endl;
			std::cout << "Printing next possible triangle. " << std::endl << std::endl;

			std::cout << std::setprecision(2) << "Side A: " << SecondaryData.m_A << std::endl;
			std::cout << std::fixed;

			std::cout << std::setprecision(2) << "Side B: " << SecondaryData.m_B << std::endl;
			std::cout << std::fixed;

			std::cout << std::setprecision(2) << "Side C: " << SecondaryData.m_C << std::endl;
			std::cout << std::fixed;

			std::cout << "Alpha: " << SecondaryData.m_Alpha << " Degrees" << std::endl;

			std::cout << "Beta: " << SecondaryData.m_Beta << " Degrees" << std::endl;

			std::cout << "Gamma: " << SecondaryData.m_Gamma << " Degrees" << std::endl;
		}
	}

private:
	Side m_A;
	Side m_B;
	Side m_C;

	Angle m_Alpha;
	Angle m_Beta;
	Angle m_Gamma;

	bool m_IsObtuse = false;
	//bool m_NotPossibleResult = false;

	SecondAnswerTriangleDataLawOfSinesSSA SecondaryData;
};

struct LawOfSinesInfoAAS
{
	LawOfSinesInfoAAS(const double Alpha,
		const double Beta,
		const double Gamma,
		const GivenSide& SideType,
		const double SideAmount)
	{
		m_Alpha = Alpha;
		m_Beta = Beta;
		m_Gamma = Gamma;


		switch (SideType)
		{
			case GivenSide::A:
			{
				SetSideA(SideAmount);
				break;
			}
			case GivenSide::B:
			{
				SetSideB(SideAmount);
				break;
			}
			case GivenSide::C:
			{
				SetSideC(SideAmount);
				break;
			}
		}

		if (m_Alpha == 0)
		{
			SetAlphaAngle();
		}
		else if (m_Beta == 0)
		{
			SetBetaAngle();
		}
		else
		{
			SetGammaAngle();
		}

		CalculateMissingSides(SideType);

	}

	void SetSideA(const double Length) { m_A = Length; }
	void SetSideB(const double Length) { m_B = Length; }
	void SetSideC(const double Length) { m_C = Length; }

	double GetAlphaAngle() const { return m_Alpha; }
	double GetGammaAngle() const { return m_Gamma; }
	double GetBetaAngle() const { return m_Beta; }

	double GetSideA() const { return m_A; }
	double GetSideB() const { return m_B; }
	double GetSideC() const { return m_C; }

	void PrintSideA() const
	{
		std::cout << std::setprecision(2) << "Side A: " << GetSideA() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideB() const
	{
		std::cout << std::setprecision(2) << "Side B: " << GetSideB() << std::endl;
		std::cout << std::fixed;
	}

	void PrintSideC() const
	{
		std::cout << std::setprecision(2) << "Side C: " << GetSideC() << std::endl;
		std::cout << std::fixed;
	}

	void PrintAllSidesNewLines() const
	{
		PrintSideA();
		PrintSideB();
		PrintSideC();
	}

	void PrintAngleAlpha() const
	{
		std::cout << "Alpha: " << GetAlphaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleBeta() const
	{
		std::cout << "Beta: " << GetBetaAngle() << " Degrees" << std::endl;
	}

	void PrintAngleGamma() const
	{
		std::cout << "Gamma: " << GetGammaAngle() << " Degrees" << std::endl;
	}

	void PrintAllAnglesNewLines() const
	{
		PrintAngleAlpha();
		PrintAngleBeta();
		PrintAngleGamma();
	}

	void PrintAllSidesAndAngles() const
	{
		PrintAllSidesNewLines();
		PrintAllAnglesNewLines();
	}


	void CalculateMissingSides(const GivenSide& Side)
	{
		switch (Side)
		{
			case GivenSide::A:
			{
				// To find an unknown side, we need to know the corresponding angle and a known ratio. 
				// We know that Angle Alpha = 50 degrees and its side a = 10
				// we can now use the law of sines proprotion
				// Sin(AlphaAngle)/SideA = Sin(GammaAngle)/C
				// then solve for C

				SetSideC((std::sin(DegreesToRadians(GetGammaAngle()))*GetSideA()) /
					std::sin(DegreesToRadians(GetAlphaAngle())));

				SetSideB((std::sin(DegreesToRadians(GetBetaAngle()))*GetSideA()) /
					std::sin(DegreesToRadians(GetAlphaAngle())));

				break;
			}
			case GivenSide::B:
			{
				SetSideC((std::sin(DegreesToRadians(GetGammaAngle()))*GetSideB()) /
					std::sin(DegreesToRadians(GetBetaAngle())));

				SetSideA((std::sin(DegreesToRadians(GetAlphaAngle()))*GetSideB()) /
					std::sin(DegreesToRadians(GetBetaAngle())));

				break;
			}
			case GivenSide::C:
			{
				SetSideA((std::sin(DegreesToRadians(GetAlphaAngle()))*GetSideC()) /
					std::sin(DegreesToRadians(GetGammaAngle())));

				SetSideB((std::sin(DegreesToRadians(GetBetaAngle()))*GetSideC()) /
					std::sin(DegreesToRadians(GetGammaAngle())));

				break;
			}
		}
	}

private:

	double m_Alpha{ 0.0 };
	double m_Beta{ 0.0 };
	double m_Gamma{ 0.0 };

	double m_A{ 0.0 };
	double m_B{ 0.0 };
	double m_C{ 0.0 };


	void SetAlphaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Alpha = 180.0 - GetBetaAngle() - GetGammaAngle();
	}

	void SetBetaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Beta = 180.0 - GetAlphaAngle() - GetGammaAngle();
	}

	void SetGammaAngle()
	{
		// The three angles must add up to 180 degrees. From this, we can determine that
		m_Gamma = 180.0 - GetAlphaAngle() - GetBetaAngle();
	}

};

double SinHalfAngleFormula2(const double& InSin, const std::pair<double, double>& InCos);

double CosHalfAngleFormula(const double& CosAngle, const std::pair<double, double>& InCos);

double TanHalfAngleFormula(const double& TanAngle, const std::pair<double, double>& InCos);



#endif TRIGFUNC