#pragma once

#include <utility>
#include <string>


typedef std::pair<double, double> CordPair;


enum class AxisOfSymmetry
{
	XAxis,
	YAxis
};

enum class ParabolaOpens
{
	UP,
	DOWN,
	LEFT,
	RIGHT
};

class Parabola
{
private:

	// if form == y^2 = 4px then == XAxis
	AxisOfSymmetry m_SymmetryAxis;
	std::string m_SymmetryAxisString;

	ParabolaOpens m_ParabolaOpens;

	CordPair m_Focus;

	CordPair m_Vertex;

	std::string m_DirectrixEquationString;

	CordPair m_LatusRectumFirst;
	CordPair m_LatusRectumSecond;

	std::string GetParabolaOpensString() const;
	std::string GetAxisOfSymmetryString() const;
	std::string PrintCordPair(const CordPair PairToPrint) const; 


public:
	explicit Parabola(const double& InCoffecient, const AxisOfSymmetry& FunctionForm);
	explicit Parabola(
		const double & h,
		const double& k,
		const double& InCoffecient,
		const AxisOfSymmetry& FunctionForm);

	~Parabola();

	void PrintParabolaInfo() const;

	CordPair GetFocus() const { return m_Focus; }
	AxisOfSymmetry GetSymmetryAxis() const { return m_SymmetryAxis; }


};


void PrintParabolaEquationInStandardForm(const CordPair& Focus, const AxisOfSymmetry& InAxis);