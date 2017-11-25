#pragma once

#ifndef ABSOLUTEVALUEFUNCTION_H
#define ABSOLUTEVALUEFUNCTION_H

#include <cmath>
#include <utility>

using std::pair;

// form y = a|x-h| + k
class AbsoluteValueFunction
{
private:

	double m_a;
	double m_h;
	double m_k;


	// x = h
	double m_AxisOfSymmetry; 

	pair<double, double> m_Vertex;

public:

	explicit AbsoluteValueFunction(const double& a, const double& h, const double& k)
		: m_a(a), m_h(h), m_k(k) 
	{

		m_Vertex.first = h;
		m_Vertex.second = k;

		if (a > 0)
		{
			// The domain of the graph is set of all real numbers and the range is y≥k  when a>0
			// The Function graph opens up
		}
		else if (a < 0)
		{
			// The domain of the graph is set of all real numbers and the range is y≤k when a<0a
			// The function Graph opens down
		}

		m_AxisOfSymmetry = m_h;

		// Note: The graph of a|x| is wider than the graph of |x| if |a| < 1, and narrower if |a| > 1

	}
	
	
	~AbsoluteValueFunction() = default;


	double operator()(const double& x) const
	{
		double FirstPart = std::abs(x - m_h);
		FirstPart = FirstPart  * m_a;
		FirstPart = FirstPart + m_k;

		return FirstPart;
	}


};

#endif