#include "PolynomialFunction.h"

#include <iostream>

void PolynomialFunction::PrintEndBehaviours() const
{
	std::cout << "as x goes to " << GetEndBehaviourNegDir().first << " f(x) goes to " << GetEndBehaviourNegDir().second << std::endl;
	std::cout << "as x goes to " << GetEndBehaviourPosDir().first << " f(x) goes to " << GetEndBehaviourPosDir().second << std::endl;
}
