#pragma once

#ifndef CIRCLE_H					// avoid repeated expansion
#define CIRCLE_H

#include <cmath>
#include <utility>
#include <map>
#include <tuple>


#include "TrigFunctions.h"


using std::tuple;



const double Zero = 0;
const double M_PIOverSix = M_PI / 6.0;
const double M_PIOverFour = M_PI / 4.0;
const double M_PIOverThree = M_PI / 3.0;
const double M_PIOverTwo = M_PI / 2.0;
const double TwoM_PIOverThree = (M_PI * 2) / 3.0;
const double ThreeM_PIOverFour = (M_PI * 3) / 4.0;
const double FiveM_PIOverSix = (M_PI * 5) / 6.0;
const double SevenM_PIOverSix = (M_PI * 7) / 6.0;
const double FiveM_PIOverFour = (M_PI * 5) / 4.0;
const double FourM_PIOverThree = (M_PI * 4) / 3.0;
const double ThreeM_PIOverTwo = (M_PI * 3) / 2.0;
const double FiveM_PIOverThree = (M_PI * 5) / 3.0;
const double SevenM_PIOverFour = (M_PI * 7) / 4.0;
const double ElevenM_PIOverSix = (M_PI * 11) / 6.0;

const double SqrtTwoOverTwo = std::sqrt(2.0) / 2;
const double NegSqrtTwoOverTwo = SqrtTwoOverTwo * (-1);
const double SqrtThreeOverTwo = std::sqrt(3.0) / 2;
const double NegSqrtThreeOverTwo = SqrtThreeOverTwo * (-1);
const double OneOverTwo = 1.0 / 2.0;
const double NegOneOverTwo = OneOverTwo * (-1);
const double OneOverSqrtThree = 1.0 / std::sqrt(3.0);
const double NegOneOverSqrtThree = (OneOverSqrtThree * (-1));


typedef std::pair<double, double> Point;



class UnitCircle
{
private:
	Point m_Center = Point(0, 0); // Origin
	unsigned int m_Radius = 1;

	std::map<double, Point> m_RadianCosSinMap;
	std::map<double, double> m_RadianTanMap;

public:

	UnitCircle()
	{
		m_RadianCosSinMap[Zero] = Point(1, 0);
		m_RadianCosSinMap[M_PIOverSix] = Point(SqrtThreeOverTwo, OneOverTwo);
		m_RadianCosSinMap[M_PIOverFour] = Point(SqrtTwoOverTwo, SqrtTwoOverTwo);
		m_RadianCosSinMap[M_PIOverThree] = Point(OneOverTwo, SqrtThreeOverTwo);
		m_RadianCosSinMap[M_PIOverTwo] = Point(0, 1);
		m_RadianCosSinMap[TwoM_PIOverThree] = Point(NegOneOverTwo, SqrtThreeOverTwo);
		m_RadianCosSinMap[ThreeM_PIOverFour] = Point(NegSqrtTwoOverTwo, SqrtTwoOverTwo);
		m_RadianCosSinMap[FiveM_PIOverSix] = Point(NegSqrtThreeOverTwo, OneOverTwo);
		m_RadianCosSinMap[M_PI] = Point(-1, 0);
		m_RadianCosSinMap[SevenM_PIOverSix] = Point(NegSqrtThreeOverTwo, NegOneOverTwo);
		m_RadianCosSinMap[FiveM_PIOverFour] = Point(NegSqrtTwoOverTwo, NegSqrtTwoOverTwo);
		m_RadianCosSinMap[FourM_PIOverThree] = Point(NegOneOverTwo, SqrtThreeOverTwo);
		m_RadianCosSinMap[ThreeM_PIOverTwo] = Point(0, -1);
		m_RadianCosSinMap[FiveM_PIOverThree] = Point(OneOverTwo, NegSqrtThreeOverTwo);
		m_RadianCosSinMap[SevenM_PIOverFour] = Point(SqrtTwoOverTwo, NegSqrtTwoOverTwo);
		m_RadianCosSinMap[ElevenM_PIOverSix] = Point(SqrtThreeOverTwo, NegOneOverTwo);

		///m_RadianTanMap[Zero] = double(m_RadianCosSinMap[Zero].second / m_RadianCosSinMap[Zero].first);
		m_RadianTanMap[M_PIOverSix] = double(m_RadianCosSinMap[M_PIOverSix].second / m_RadianCosSinMap[M_PIOverSix].first);
		m_RadianTanMap[M_PIOverFour] = double(m_RadianCosSinMap[M_PIOverFour].second / m_RadianCosSinMap[M_PIOverFour].first);
		m_RadianTanMap[M_PIOverThree] = double(m_RadianCosSinMap[M_PIOverThree].second / m_RadianCosSinMap[M_PIOverThree].first);
		//m_RadianTanMap[M_PIOverTwo] = double(m_RadianCosSinMap[M_PIOverTwo].second / m_RadianCosSinMap[M_PIOverTwo].first);
		m_RadianTanMap[TwoM_PIOverThree] = double(m_RadianCosSinMap[TwoM_PIOverThree].second / m_RadianCosSinMap[TwoM_PIOverThree].first);
		m_RadianTanMap[ThreeM_PIOverFour] = double(m_RadianCosSinMap[ThreeM_PIOverFour].second / m_RadianCosSinMap[ThreeM_PIOverFour].first);
		m_RadianTanMap[FiveM_PIOverSix] = double(m_RadianCosSinMap[FiveM_PIOverSix].second / m_RadianCosSinMap[FiveM_PIOverSix].first);
		//m_RadianTanMap[M_PI] = double(m_RadianCosSinMap[M_PI].second / m_RadianCosSinMap[M_PI].first);
		m_RadianTanMap[SevenM_PIOverSix] = double(m_RadianCosSinMap[SevenM_PIOverSix].second / m_RadianCosSinMap[SevenM_PIOverSix].first);
		m_RadianTanMap[FiveM_PIOverFour] = double(m_RadianCosSinMap[FiveM_PIOverFour].second / m_RadianCosSinMap[FiveM_PIOverFour].first);
		m_RadianTanMap[FourM_PIOverThree] = double(m_RadianCosSinMap[FourM_PIOverThree].second / m_RadianCosSinMap[FourM_PIOverThree].first);
		//m_RadianTanMap[ThreeM_PIOverTwo] = double(m_RadianCosSinMap[ThreeM_PIOverTwo].second / m_RadianCosSinMap[ThreeM_PIOverTwo].first);
		m_RadianTanMap[FiveM_PIOverThree] = double(m_RadianCosSinMap[FiveM_PIOverThree].second / m_RadianCosSinMap[FiveM_PIOverThree].first);
		m_RadianTanMap[SevenM_PIOverFour] = double(m_RadianCosSinMap[SevenM_PIOverFour].second / m_RadianCosSinMap[SevenM_PIOverFour].first);
		m_RadianTanMap[ElevenM_PIOverSix] = double(m_RadianCosSinMap[ElevenM_PIOverSix].second / m_RadianCosSinMap[ElevenM_PIOverSix].first);
	}

	std::map<double, Point> GetRadianSinCosMap() const { return m_RadianCosSinMap; }
	std::map<double, double> GetRadianTanMap() const { return m_RadianTanMap; }
};

class CircleFunction
{
private:
	//Point m_Center;


	double m_x;
	double m_h;
	
	double m_y;
	double m_k;

	double m_r;



public:
	CircleFunction() = default;

	explicit CircleFunction(const double& x, const double& h, const double& y, const double& k, const double &r)
		: m_x{ x }, m_h{ h }, m_y{ y }, m_k{ k }, m_r{ r }
	{


	}

	std::tuple<double, double, double, double, double> GetXHYKR() const;


	/* INPUT: Value whose arc sine is computed, in the interval[-1, +1].
	If the argument is out of this interval, a domain error occurs. */
	/* RETURN: Principal arc sine of x, in the interval [-pi/2,+pi/2] radians.
	One radian is equivalent to 180/PI degrees.*/
	double operator()(const double& x, const double& y)
	{

		double LocalX = std::pow(x - m_h, 2);
		double LocalY = std::pow(y - m_k, 2);

		double LHS = LocalX + LocalY;

		double RHS = std::pow(m_r, 2);

		LHS == RHS;

		return 0;
	}

};

// Let P=(x,y) be a point on the unit circle centered at the origin O.
// Let θ be an angle with an initial side along the positive x-axis 
// and a terminal side given by the line segment OP
// If x=0,secθ and tanθ are undefined. If y=0, then cotθ and cscθ are undefined.


#endif