#include "MPVector.h"

#include <cmath>
#include <iostream>
#include <string>

const double MathPIE = 3.14159;

MPVector::MPVector(const double& x, const double &y)
{
	m_x = x;
	m_y = y;
}

MPVector::MPVector(const double & x1, const double & y1, const double & x2, const double & y2)
{
	// The position vector begins at (0,0) and terminates at (m_x, m_y)
	m_x = x2 - x1;
	m_y = y2 - y1;
}

double MPVector::GetMagnitude() const
{
	double OutMagnitude{ 0.0 };
	OutMagnitude = std::sqrt(std::pow(m_x, 2) + std::pow(m_y, 2));
	return OutMagnitude;
}

double MPVector::GetDirection() const
{
	double OutDirection{ 0.0 };
	OutDirection = std::atan(m_y / m_x);
	OutDirection = (OutDirection * 180) / MathPIE;

	while (OutDirection < 0)
	{
		OutDirection = OutDirection + 360.0;
	}

	while (OutDirection > 360)
	{
		OutDirection = OutDirection - 360;
	}

	return OutDirection;
}

bool MPVector::operator!=(const MPVector& other)
{
	return ((m_x != other.m_x) || m_y != other.m_y);
}

bool MPVector::operator==(const MPVector& other)
{
	return !(operator!=(other));
}

MPVector MPVector::MultiplyByScaler(const double& Scalar)
{
	m_x = m_x*Scalar;
	m_y = m_y*Scalar;
	return *this;
}

MPVector MPVector::GetUnitVector() const
{
	MPVector OutUnitVector(m_x, m_y);
	OutUnitVector.m_x = m_x / GetMagnitude();
	OutUnitVector.m_y = m_y / GetMagnitude();

	return OutUnitVector;
}


std::string MPVector::GetVectorComponentInTermsOfMagAndDir()
{
	double Angle = std::atan(m_y / m_x);
	const double AngleInDeg = (Angle*180.0) / MathPIE;

	std::string OutString = "y = ";
	const std::string MagnitudeString = std::to_string(GetMagnitude());
	const std::string AngleString = std::to_string(AngleInDeg);

	std::string MagStringCos = MagnitudeString;
	std::string MagStringSin = MagnitudeString;

	std::string AngleStringCos = AngleString;
	std::string AngleStringSin = AngleString;

	std::string CosString = "cos";
	std::string SinString = "sin";
	
	OutString += MagStringCos += CosString += AngleStringCos += "i + ";
	OutString += MagStringSin += SinString += AngleStringSin += "j";

	return OutString;
}

double MPVector::GetDotProductOfTwoVectors(const MPVector & other) const
{
	double OutResult{ 0.0 };
	OutResult = (m_x * other.m_x) + (m_y * other.m_y);
	return OutResult;
}

double MPVector::GetAngleBetweenTwoVectors(const MPVector & other) const
{
	MPVector v = GetUnitVector();
	MPVector u = other.GetUnitVector();

	double InRad = (v.m_x * u.m_x) + (v.m_y * u.m_y);
	InRad = std::acos(InRad);

	// return in degrees
	double ReturnedDegrees = (InRad * 180.0) / MathPIE;

	return ReturnedDegrees;
}


std::string GetVectorComponentInTermsOfMagAndDir(
	const double & Angle,
	const double & Magnitude)
{
	std::string OutString = "y =";

	const double AngleInRadians = (MathPIE * Angle) / 180.0;

	const double x = Magnitude * (std::cos(AngleInRadians));
	OutString += std::to_string(x);
	OutString += "i +";

	const double y = Magnitude * (std::sin(AngleInRadians));
	OutString += std::to_string(y);
	OutString += "j";

	return OutString;
}

void OutputVectorRectangularForm(const double & x1, const double & y1, const double & x2, const double & y2)
{
	// v = ai + bj
	
	double a = x2 - x1;
	double b = y2 - y1;

	double OutMagnitude{ 0.0 };
	OutMagnitude = std::sqrt(std::pow(a, 2) + std::pow(b, 2));

	std::cout << "v = " << a << 'i' << "+" << b << "j" << std::endl;
	std::cout << "Magnitude = " << OutMagnitude << std::endl;
}

void GetComponentsFromInitalAndTerminalPoints(const double & x1, const double & y1, const double & x2, const double & y2, double & OutHorizComponentX, double & OutHorizComponentY, double & OutVertiComponentX, double & OutVertiComponentY)
{
	double StandardPosX = x2 - x1;
	double StandardPosY = y2 - y1;

	OutHorizComponentX = StandardPosX;
	OutHorizComponentY = 0;

	OutVertiComponentX = 0;
	OutVertiComponentY = StandardPosY;
}

std::ostream& operator<<(std::ostream& Stream, const MPVector &InVec)
{
	// TODO: insert return statement here
	Stream << "(" << InVec.m_x << "," << InVec.m_y << ")";
	return Stream;
}

MPVector operator+(const MPVector &v1, const MPVector &v2)
{
	return MPVector(v1.m_x + v2.m_x, v1.m_y + v2.m_y);
}

MPVector operator-(const MPVector &v1, const MPVector &v2)
{
	return MPVector(v1.m_x - v2.m_x, v1.m_y - v2.m_y);
}