#pragma once

#include <cmath>

class ComplexNumber
{
public:

	explicit ComplexNumber(double real, bool AddImag, double imag);
	explicit ComplexNumber(double real, double imag);
	~ComplexNumber();

	// Given  z=x+yi, a complex number, the absolute value of z is defined as
	// |z| = sqrt(x^2+y^2)
	// It is the distance from the origin to the point (x,y).
	// Notice that the absolute value of a real number gives the distance of the number from 0, 
	// while the absolute value of a complex number gives the distance
	// of the number from the origin, (0, 0).
	void UpdateAbsValue();
	const double GetAbsValue() const { return m_AbsValue; }


	// Finding Powers of Complex Numbers in Polar Form
	// If z=r(cos θ+isin θ)  z=r(cos θ+isin θ) is a complex number, then
	// z^n = r^n((cis)ntheta)
	// Evaluate the expression (1 + i)^n using De Moivre’s Theorem.
	ComplexNumber EvaluatePowerExpressionDeMovire(const double& n);

	// Finding Roots of Complex Numbers in Polar Form
	// To find the nth root of a complex number in polar form, we use the nth  nth Root Theorem or
	// De Moivre’s Theorem and raise the complex number to a power with a rational exponent.
	ComplexNumber FindNthRoot(const double& n, const double& theta,
		const double& r, const double& k = 0) const;


	double GetReal() const { return m_Real; }
	double GetImag() const { return m_Imaginary; }

private:
	double m_Real{ 0.0 };
	double m_Imaginary{ 0.0 };

	double m_AbsValue{ 0.0 };
};

