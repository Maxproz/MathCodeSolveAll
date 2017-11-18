#include "Ellipse.h"


#include <string>
#include <iostream>

CordPair Ellipse::GetCenter()
{
	// the center is halfway between the verticies
	CordPair Center;
	Center.first = (m_FirstVertex.first + m_SecondVertex.first) / 2.0;
	Center.second = (m_FirstVertex.second + m_SecondVertex.second) / 2.0;
	
	m_Center = Center;
	return Center;
}


void Ellipse::OutputFunctionNotation()
{
	if (m_MajorAxis == MajorAxis::YAxis)
	{

		std::string outstring;
		std::string HString = std::to_string(m_Center.first);
		std::string KString = std::to_string(m_Center.second);
		std::string bSquaredString = std::to_string(m_bSquared);
		std::string aSquaredString = std::to_string(m_aSquared);

		outstring.append("(x - ");
		outstring.append(HString);
		outstring.append(")^2 / ");
		outstring.append(bSquaredString);
		outstring.append(" + (y - ");
		outstring.append(KString);
		outstring.append(")^2 / ");
		outstring.append(aSquaredString);
		outstring.append(" =  1");

		std::cout << outstring << std::endl;
	}
	else
	{
		std::string outstring;
		std::string HString = std::to_string(m_Center.first);
		std::string KString = std::to_string(m_Center.second);
		std::string bSquaredString = std::to_string(m_bSquared);
		std::string aSquaredString = std::to_string(m_aSquared);

		outstring.append("(x - ");
		outstring.append(HString);
		outstring.append(")^2 / ");
		outstring.append(aSquaredString);
		outstring.append(" + (y - ");
		outstring.append(KString);
		outstring.append(")^2 / ");
		outstring.append(bSquaredString);
		outstring.append(" =  1");

		std::cout << outstring << std::endl;


	}
}

// Centerd at origin
EllipseInfo GetEllipseInfoFromGivenEquationAtOrigin(const double& aSquared, const double& bSquared, const MajorAxis& InAxis)
{
	EllipseInfo NewInfo;
	
	CordPair Center(0,0);

	CordPair FirstVertex;
	CordPair SecondVertex;

	CordPair FirstCoVertex;
	CordPair SecondCoVertex;

	CordPair FirstFoci;
	CordPair SecondFoci;


	double a = std::sqrt(aSquared);
	double b = std::sqrt(bSquared);
	double c{ 0.0 };

	if (aSquared > bSquared)
	{
		if (InAxis == MajorAxis::XAxis) // Bigger number on left
		{

			FirstVertex.first = -1 * a;
			FirstVertex.second = 0;

			SecondVertex.first = a;
			SecondVertex.second = 0;

			FirstCoVertex.first = 0;
			FirstCoVertex.second = -1 * b;

			SecondCoVertex.first = 0;
			SecondCoVertex.second = b;

			c = std::sqrt(aSquared - bSquared);

			FirstFoci.first = c * -1;
			FirstFoci.second = 0;

			SecondFoci.first = c;
			SecondFoci.second = 0;


			NewInfo["Center"] = Center;
			NewInfo["FirstVertex"] = FirstVertex;
			NewInfo["SecondVertex"] = SecondVertex;
			NewInfo["FirstCoVertex"] = FirstCoVertex;
			NewInfo["SecondCoVertex"] = SecondCoVertex;
			NewInfo["FirstFoci"] = FirstFoci;
			NewInfo["SecondFoci"] = SecondFoci;

		}
		else
		{
			// MajorAxis ==  YAxis // bigger number on right (swap)

			// Vertex are at (0, +-a)
			FirstVertex.first = 0;
			FirstVertex.second = -1 * a;

			SecondVertex.first = 0;
			SecondVertex.second = a;

			// CoVertex are at (+-b,0)
			FirstCoVertex.first = -1 * b;
			FirstCoVertex.second = 0;

			SecondCoVertex.first = b;
			SecondCoVertex.second = 0;

			// Foci are at (0, +-c)
			// find c
			// c^2 = a^2 - b^2
			c = std::sqrt(aSquared - bSquared);

			FirstFoci.first = 0;
			FirstFoci.second = c * -1;

			SecondFoci.first = 0;
			SecondFoci.second = c;


			NewInfo["Center"] = Center;
			NewInfo["FirstVertex"] = FirstVertex;
			NewInfo["SecondVertex"] = SecondVertex;
			NewInfo["FirstCoVertex"] = FirstCoVertex;
			NewInfo["SecondCoVertex"] = SecondCoVertex;
			NewInfo["FirstFoci"] = FirstFoci;
			NewInfo["SecondFoci"] = SecondFoci;
		}
	}
	else
	{
		throw std::exception("a is not greater than b");
	}


	return NewInfo;
}

EllipseInfo GetEllipseInfoFromGivenEquationAtOrigin(const double & XCoffeicent, const double & YCoffeicent, const double & EqualTo)
{

	EllipseInfo OutRes;

	double Temp1 = EqualTo / XCoffeicent;
	double Temp2 = EqualTo / YCoffeicent;

	if (Temp1 > Temp2)
	{
		OutRes = GetEllipseInfoFromGivenEquationAtOrigin(Temp1, Temp2, MajorAxis::XAxis);
	}
	else
	{
		OutRes = GetEllipseInfoFromGivenEquationAtOrigin(Temp2, Temp1, MajorAxis::YAxis);
	}


	return OutRes;
}

EllipseInfo GetEllipseInfoFromGivenEquation(const double & H, const double & K, const double & aSquared, const double & bSquared, const MajorAxis & InAxis)
{
	EllipseInfo NewInfo;

	CordPair Center(H, K);

	CordPair FirstVertex;
	CordPair SecondVertex;

	CordPair FirstCoVertex;
	CordPair SecondCoVertex;

	CordPair FirstFoci;
	CordPair SecondFoci;


	double a = std::sqrt(aSquared);
	double b = std::sqrt(bSquared);
	double c{ 0.0 };

	if (aSquared > bSquared)
	{
		if (InAxis == MajorAxis::XAxis) // Bigger number on left
		{
			FirstVertex.first = H + a;
			FirstVertex.second = K;

			SecondVertex.first = H - a;
			SecondVertex.second = K;

			FirstCoVertex.first = H;
			FirstCoVertex.second = K + b;

			SecondCoVertex.first = H;
			SecondCoVertex.second = K - b;

			c = std::sqrt(aSquared - bSquared);

			FirstFoci.first = H + c;
			FirstFoci.second = K;

			SecondFoci.first = H - c;
			SecondFoci.second = K;

			NewInfo["Center"] = Center;
			NewInfo["FirstVertex"] = FirstVertex;
			NewInfo["SecondVertex"] = SecondVertex;
			NewInfo["FirstCoVertex"] = FirstCoVertex;
			NewInfo["SecondCoVertex"] = SecondCoVertex;
			NewInfo["FirstFoci"] = FirstFoci;
			NewInfo["SecondFoci"] = SecondFoci;
		}
		else
		{
			// MajorAxis == Y::Axis
			FirstVertex.first = H;
			FirstVertex.second = K + a;

			SecondVertex.first = H;
			SecondVertex.second = K - a;

			FirstCoVertex.first = H + b;
			FirstCoVertex.second = K;

			SecondCoVertex.first = H - b;
			SecondCoVertex.second = K;

			c = std::sqrt(aSquared - bSquared);

			FirstFoci.first = H;
			FirstFoci.second = K + c;

			SecondFoci.first = H;
			SecondFoci.second = K - c;

			NewInfo["Center"] = Center;
			NewInfo["FirstVertex"] = FirstVertex;
			NewInfo["SecondVertex"] = SecondVertex;
			NewInfo["FirstCoVertex"] = FirstCoVertex;
			NewInfo["SecondCoVertex"] = SecondCoVertex;
			NewInfo["FirstFoci"] = FirstFoci;
			NewInfo["SecondFoci"] = SecondFoci;
		}
	}
	else
	{
		throw std::exception("a is not greater than b");
	}

	return NewInfo;
}
