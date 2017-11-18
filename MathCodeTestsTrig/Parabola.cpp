#include "Parabola.h"

#include <exception>
#include <iostream>
#include <iomanip>

void Parabola::PrintParabolaInfo() const
{
	std::cout << "Printing basic parabola function info with vertex at origin." << std::endl;
	std::cout << "Axis Of Symmetry: " << GetAxisOfSymmetryString() << std::endl;
	std::cout << "The Parabola Opens: " << GetParabolaOpensString() << std::endl;
    std::cout << "The Directrix Function Equation is: " << m_DirectrixEquationString << std::endl;
	std::cout << "The Focus is at: " << PrintCordPair(m_Focus) << std::endl;
	std::cout << "The First Latus Rectum point is at: " << PrintCordPair(m_LatusRectumFirst) << std::endl;
	std::cout << "The Second Latus Rectum point is at: " << PrintCordPair(m_LatusRectumSecond) << std::endl;
	std::cout << "End of current data" << std::endl;
}

Parabola::Parabola(const double & InCoffecient, const AxisOfSymmetry & FunctionForm)
{
	m_SymmetryAxis = FunctionForm;

	double p{ 0.0 };

	// origin
	m_Vertex = CordPair(0, 0);

	if (FunctionForm == AxisOfSymmetry::XAxis)
	{
		p = InCoffecient / 4.0;

		if (p > 0)
		{
			m_ParabolaOpens = ParabolaOpens::RIGHT;
		}
		else
		{
			m_ParabolaOpens = ParabolaOpens::LEFT;
		}

		m_Focus = CordPair(p, 0);

		m_DirectrixEquationString = std::string("x = ");
		double negp = p*-1;
		std::string NegPString = std::to_string(negp);
		while (NegPString.size() >= 5)
		{
			NegPString.pop_back();
		}

		m_DirectrixEquationString.append(NegPString);

		m_LatusRectumFirst = CordPair(p, 2 * p);
		m_LatusRectumSecond = CordPair(p, ((2 * p) * -1));

	}
	else if (FunctionForm == AxisOfSymmetry::YAxis)
	{
		p = InCoffecient / 4.0;

		if (p > 0)
		{
			m_ParabolaOpens = ParabolaOpens::UP;
		}
		else
		{
			m_ParabolaOpens = ParabolaOpens::DOWN;
		}

		m_Focus = CordPair(0, p);

		m_DirectrixEquationString = std::string("y = ");
		double negp = p*-1;
		std::string NegPString = std::to_string(negp);
		while (NegPString.size() >= 5)
		{
			NegPString.pop_back();
		}

		m_DirectrixEquationString.append(NegPString);

		m_LatusRectumFirst = CordPair(2 * p, p);
		m_LatusRectumSecond = CordPair(((2 * p) * -1), p);

	}
	else
	{
		throw std::exception("Invalid FunctionForm param");
	}

}

Parabola::Parabola(const double & h, const double & k, const double & InCoffecient, const AxisOfSymmetry & FunctionForm)
{
	m_Vertex = CordPair(h, k);

	m_SymmetryAxis = FunctionForm;

	double p{ 0.0 };

	p = InCoffecient / 4.0;

	if (FunctionForm == AxisOfSymmetry::XAxis)
	{
		m_SymmetryAxisString = "y = ";
		m_SymmetryAxisString.append(std::to_string(k));

		m_DirectrixEquationString = std::string("x = ");
		double negp = h - p;
		std::string NegPString = std::to_string(negp);
		while (NegPString.size() >= 5)
		{
			NegPString.pop_back();
		}

		m_DirectrixEquationString.append(NegPString);

		if (p > 0)
		{
			m_ParabolaOpens = ParabolaOpens::RIGHT;
		}
		else
		{
			m_ParabolaOpens = ParabolaOpens::LEFT;
		}


		m_Focus = CordPair(h + p, k);

		m_LatusRectumFirst = CordPair(h + p, k + (2*p));
		m_LatusRectumSecond = CordPair(h + p, k - (2*p));
	}
	else if (FunctionForm == AxisOfSymmetry::YAxis)
	{
		m_SymmetryAxisString = std::string("x = ");
		m_SymmetryAxisString.append(std::to_string(h));

		m_DirectrixEquationString = std::string("y = ");
		double negp = k - p;
		std::string NegPString = std::to_string(negp);
		while (NegPString.size() >= 5)
		{
			NegPString.pop_back();
		}

		m_DirectrixEquationString.append(NegPString);

		if (p > 0)
		{
			m_ParabolaOpens = ParabolaOpens::UP;
		}
		else
		{
			m_ParabolaOpens = ParabolaOpens::DOWN;
		}

		m_Focus = CordPair(h, k + p);

		m_LatusRectumFirst = CordPair(h + (2 * p), k + p);
		m_LatusRectumSecond = CordPair(h - (2 * p), k + p);
	}
	else
	{
		throw std::exception("Invalid FunctionForm param");
	}
}

Parabola::~Parabola()
{

}

std::string Parabola::GetParabolaOpensString() const
{
	switch (m_ParabolaOpens)
	{
		case ParabolaOpens::UP:
		{
			return "Up";
		}
		case ParabolaOpens::DOWN:
		{
			return "Down";
		}
		case ParabolaOpens::LEFT:
		{
			return "Left";
		}
		case ParabolaOpens::RIGHT:
		{
			return "Right";
		}
	}
}

std::string Parabola::GetAxisOfSymmetryString() const
{
	return m_SymmetryAxisString;
}

std::string Parabola::PrintCordPair(const CordPair PairToPrint) const
{
	std::string OutString;
	OutString.append("(");
	std::string FirstNumString = std::to_string(PairToPrint.first);
	while (FirstNumString.size() >= 5)
	{
		FirstNumString.pop_back();
	}
	OutString.append(FirstNumString);
	OutString.append(", ");
	std::string SecondNumString = std::to_string(PairToPrint.second);
	while (SecondNumString.size() >= 5)
	{
		SecondNumString.pop_back();
	}
	OutString.append(SecondNumString);
	OutString.append(")");

	return OutString;
}

void PrintParabolaEquationInStandardForm(const CordPair & Focus, const AxisOfSymmetry & InAxis)
{
	if (InAxis == AxisOfSymmetry::XAxis)
	{
		double TempVar = Focus.first * 4;

		std::cout << "y^2 = " << TempVar << "x" << std::endl;
	}
	else
	{
		double TempVar = Focus.second * 4;

		std::cout << "x^2 = " << TempVar << "y" << std::endl;
	}
}
