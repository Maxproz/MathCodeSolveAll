#pragma once

#include <utility>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>

typedef std::pair<double, double> CordPair;
typedef std::map<std::string, std::pair<double, double>> HyperbolaInfo;

enum class HyperbolaMajorAxis
{
	XAxis,
	YAxis
};

class Hyperbola
{
private:
	typedef std::pair<double, double> CordPair;

	CordPair m_FirstVertex;
	CordPair m_SecondVertex;

	CordPair m_FirstFoci;
	CordPair m_SecondFoci;

	CordPair m_FirstCoVertex;
	CordPair m_SecondCoVertex;

	HyperbolaMajorAxis m_MajorAxis;

	CordPair m_Center;

	HyperbolaInfo m_CurrentInfo;

	// the length of the transverse axis is 2a
	double m_TransverseAxisLength;

	// the length of the conjugate axis is 2b
	double m_ConjugateAxisLength;

	// updates the center and returns it incase i need it
	CordPair GetCenter();

	double m_aSquared;
	double m_cSquared;
	double m_bSquared;

public:

	Hyperbola(CordPair FirstVert, CordPair SecondVert, CordPair FirstFoci, CordPair SecondFoci)
		: m_FirstVertex(FirstVert), m_SecondVertex(SecondVert), m_FirstFoci(FirstFoci), m_SecondFoci(SecondFoci)
	{

		CordPair Center = GetCenter();
		double h = Center.first;
		double k = Center.second;

		double x1, x2, x3, x4;
		x1 = FirstVert.first;
		x2 = SecondVert.first;
		x3 = FirstFoci.first;
		x4 = SecondFoci.first;

		double c{ 0.0 };

		if (x1 == x2 && x1 == x3 && x1 == x4)
		{
			m_MajorAxis = HyperbolaMajorAxis::YAxis;

			// Next, we find a2. The length of the major axis 2a,
			// is bounded by the vertices. 
			// We solve for a by finding the distance between
			// the y-coordinates of the vertices.
			double a = ((SecondVert.second - FirstVert.second) / 2.0);
			m_aSquared = std::pow(a, 2);

			if (k < 0)
			{
				k = k * -1;
				c = m_SecondFoci.second + k;
			}
			else
			{
				k = k * -1;
				c = m_FirstFoci.second + k;
				c = c * -1;
			}

			m_cSquared = std::pow(c, 2);
			m_bSquared = m_aSquared - m_cSquared;
		}
		else
		{
			m_MajorAxis = HyperbolaMajorAxis::XAxis;

			double a = ((SecondVert.first - FirstVert.first) / 2.0);
			m_aSquared = std::pow(a, 2);

			if (h < 0)
			{
				h = h * -1;
				c = m_SecondFoci.first + h;
			}
			else
			{
				h = h * -1;
				c = m_FirstFoci.first + h;
				c = c * -1;
			}

			m_cSquared = std::pow(c, 2);
			m_bSquared = m_aSquared - m_cSquared;
		}
	}

	explicit Hyperbola(HyperbolaInfo& Info)
	{
		m_FirstVertex = Info["FirstVertex"];
		m_SecondVertex = Info["SecondVertex"];

		m_FirstCoVertex = Info["FirstCoVertex"];
		m_SecondCoVertex = Info["SecondCoVertex"];

		m_FirstFoci = Info["FirstFoci"];
	    m_SecondFoci = Info["SecondFoci"];

		if (Info["MajorAxis"] == CordPair(1, 0))
		{
			m_MajorAxis = HyperbolaMajorAxis::XAxis;
		}
		else
		{
			m_MajorAxis = HyperbolaMajorAxis::YAxis;
		}

		m_Center = Info["Center"];

		// the length of the transverse axis is 2a
		m_TransverseAxisLength = Info["Transverse/ConjugateAxisLengths"].first;

		// the length of the conjugate axis is 2b
		m_ConjugateAxisLength = Info["Transverse/ConjugateAxisLengths"].second;

		m_CurrentInfo = Info;
	}

	~Hyperbola() = default;

	void OutputFunctionNotation();

	HyperbolaInfo GetCurrentInfo() const { return m_CurrentInfo; }
	//void PrintAsymtoteEquation();

};

template<typename Map>
void print_map(Map& m);
//{
//	//std::cout << '{';
//	for (auto& p : m)
//	{
//		std::cout << p.first << ": (" << p.second.first << ", " << p.second.second << ")\n";
//	}
//	//std::cout << "}\n";
//}


HyperbolaInfo GetHyperbolaInfoFromGivenEquationAtOrigin(const double& aSquared, const double& bSquared, const HyperbolaMajorAxis& InAxis);

//HyperbolaInfo GetEllipseInfoFromGivenEquationAtOrigin(const double& XCoffeicent, const double& YCoffeicent, const double& EqualTo);
//
HyperbolaInfo GetHyperbolaInfoFromGivenEquation(
	const double& H,
	const double& K,
	const double& aSquared,
	const double& bSquared,
	const HyperbolaMajorAxis& InAxis);








