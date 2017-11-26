#include "FunctionEnums.h"

#include <iostream>

using std::cout;
//using std::endl;

void PrintDiscontinunityType(const DiscontinunityType& InType)
{
	switch (InType)
	{
		case DiscontinunityType::REMOVEABLE:
		{
			std::cout << "REMOVEABLE\n";
			break;
		}
		case DiscontinunityType::JUMP:
		{
			std::cout << "JUMP\n";
			break;
		}
		case DiscontinunityType::INFINITE:
		{
			std::cout << "INFINITE\n";
			break;
		}
		default:
		{
			throw std::logic_error("discontinuity assignment error in the print func");
		}
	}
}