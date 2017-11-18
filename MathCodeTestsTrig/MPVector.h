#pragma once

//#include <utility>

#include <ostream>

class MPVector
{
public:

	explicit MPVector(const double& x, const double &y);
	explicit MPVector(const double& x1, const double& y1, const double& x2, const double& y2);
	MPVector(const MPVector&) = default; // copy constructor
	MPVector(MPVector&&) = default; // move constructor
	MPVector& operator=(const MPVector&) = default; // copy assign
	MPVector& operator=(MPVector&&) = default; // move assign
	virtual ~MPVector() = default;

	double GetMagnitude() const;
	double GetDirection() const;

	friend std::ostream& operator<<(std::ostream& Stream, const MPVector& InVec);
	
	bool operator!=(const MPVector& other);
	bool operator==(const MPVector& other);

	// The sum of two vectors u and v, or vector addition,
	// produces a third vector u+ v, the resultant vector.
	friend MPVector operator+(const MPVector &v1, const MPVector &v2);

	// Vector subtraction is similar to vector addition.
	// To find u − v, view it as u + (−v). 
	friend MPVector operator-(const MPVector &v1, const MPVector &v2);

	// Multiplying a vector by a scalar, a constant,
	// changes only the magnitude of the vector or the length of the line
	// Only the magnitude changes,
	// unless k is negative, if that is true then the vector reverses direction.
	MPVector MultiplyByScaler(const double& Scalar);

	//  The horizontal unit vector is written as i=⟨1,0⟩
	// and is directed along the positive horizontal axis. 
	// he vertical unit vector is written j=⟨0,1⟩
	// and is directed along the positive vertical axis.
	// If v  is a nonzero vector, then v/|v| is a unit vector in the direction of v.
	MPVector GetUnitVector() const;

	// Some basic notes
	// v = ai + bj 
	// mag = sqrt(a^2 + b^2)
	// v+u=(a+c)i+(b+d)j
	// v−u=(a−c)i+(b−d)j

	std::string GetVectorComponentInTermsOfMagAndDir();

	double GetDotProductOfTwoVectors(const MPVector& other) const;

	double GetAngleBetweenTwoVectors(const MPVector& other) const;


private:
	double m_x;
	double m_y;

};

std::string GetVectorComponentInTermsOfMagAndDir(
	const double& Angle,
	const double & Magnitude);

void OutputVectorRectangularForm(
	const double& x1,
	const double& y1,
	const double& x2,
	const double& y2);

void GetComponentsFromInitalAndTerminalPoints(
	const double& x1,
	const double& y1,
	const double& x2,
	const double& y2,
	double& OutHorizComponentX,
	double& OutHorizComponentY,
	double& OutVertiComponentX,
	double& OutVertiComponentY);

std::ostream& operator<<(std::ostream& Stream, const MPVector &InVec);

MPVector operator+(const MPVector &v1, const MPVector &v2);

MPVector operator-(const MPVector &v1, const MPVector &v2);