#include "Circle.h"

std::tuple<double, double, double, double, double> CircleFunction::GetXHYKR() const
{
	return std::tuple<double, double, double, double, double>(m_x, m_h, m_y, m_k, m_r);
}
