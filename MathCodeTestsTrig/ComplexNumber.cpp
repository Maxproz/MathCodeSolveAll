#include "stdafx.h"
#include "ComplexNumber.h"

constexpr double MATH_PI = 3.14159265359;

ComplexNumber::ComplexNumber(double real, bool AddImag, double imag)
{
	if (imag == 0)
	{
		if (AddImag == true)
		{
			m_Imaginary = 1;
		}
		else
		{
			m_Imaginary = -1;
		}
	}
	else
	{
		if (AddImag == true)
		{
			m_Imaginary = imag * 1;
		}
		else
		{
			m_Imaginary = imag * (-1);
		}
	}

	m_Real = real;

	UpdateAbsValue();
}
ComplexNumber::ComplexNumber(double real, double imag)
{
	m_Imaginary = imag;
	m_Real = real;

	UpdateAbsValue();
}



ComplexNumber::~ComplexNumber()
{

}

void ComplexNumber::UpdateAbsValue()
{
	m_AbsValue = std::sqrt(std::pow(m_Real, 2) + std::pow(m_Imaginary, 2));
}

ComplexNumber ComplexNumber::EvaluatePowerExpressionDeMovire(const double & n)
{
	// Let us find r.
	double r = GetAbsValue();

	// Then we find θ. Using the formula tan θ=y/x  gives ~
	// tanθ = 1/1
	// tanθ = 1
	// θ = tan^-1(1) == pi/4
	double ThetaAng = std::atan(m_Imaginary / m_Real);

	// Use De Moivre’s Theorem to evaluate the expression.
	//(a + bi)^n = r^n[cos(nθ) + isin(nθ)]
	double OutReal = std::pow(r, n) * (std::cos(n * ThetaAng));
	double OutImag = std::pow(r, n) * (std::sin(n * ThetaAng));

	return ComplexNumber(OutReal, OutImag);

}

ComplexNumber ComplexNumber::FindNthRoot(const double & n, const double& theta,
	const double& r,
	const double& k) const
{
	double Updatedr = std::pow(r, (1.0/n));
	//double ThetaAng = std::atan(std::cos(m_Imaginary) / std::sin(m_Real));
	//std::polar(m_Real, m_Imaginary)
	double TempOutReal = 0.0;
	double TempOutImag = 0.0;
		
	double InsideCos = (theta / n);
	double InsideSin = (theta / n);

	if (k == 0)
	{
		InsideCos = InsideCos + ((2 * 0 * MATH_PI) / n);
		InsideSin = InsideSin + ((2 * 0 * MATH_PI) / n);
		
		//OutAngle = r * (std::cos(InsideCos))  std::sin(InsideSin);

	}
	else if (k == 1)
	{
		InsideCos = (InsideCos + ((2.0 * 1.0 * MATH_PI) / n));
		InsideSin = (InsideSin + ((2.0 * 1.0 * MATH_PI) / n));

	}
	else if (k == 2)
	{
		InsideCos = (InsideCos + ((2.0 * 2.0 * MATH_PI) / n));
		InsideSin = (InsideSin + ((2.0 * 2.0 * MATH_PI) / n));


	}
	else if (k == 3)
	{
		InsideCos = (InsideCos + ((2.0 * 3.0 * MATH_PI) / n));
		InsideSin = (InsideSin + ((2.0 * 3.0 * MATH_PI) / n));


	}
	else if (k == 4)
	{
		InsideCos = (InsideCos + ((2.0 * 4.0 * MATH_PI) / n));
		InsideSin = (InsideSin + ((2.0 * 4.0 * MATH_PI) / n));


	}
	else
	{
		// not implemented
	}

	TempOutReal = Updatedr * std::cos(InsideCos);
	TempOutImag = Updatedr * std::sin(InsideSin);

	ComplexNumber OutComplex(TempOutReal, TempOutImag);

	return OutComplex;
}
