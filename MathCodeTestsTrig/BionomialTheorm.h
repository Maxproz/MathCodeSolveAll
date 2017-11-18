#pragma once

#include <string>

// https://cnx.org/contents/E6wQevFf@6.38:Bezm7bOH@9/Binomial-Theorem
// A polynomial with two terms is called a binomial.

uint64_t BTFactorial(uint64_t n);

uint64_t GetPossibleCombinationsBT(const int & nTotalAmount, const int & rSelectAmount);

// (3x - y)^4
// This happens because (−y) raised to odd powers is negative, 
// but (−y) raised to even powers is positive. 
// This will occur whenever the binomial contains a subtraction sign.
// + - + - ... etc
std::string DisplaySimpleBionomialTherom(const int& n, const bool& bIsPlus = true);

std::string FindASingleTerm(const int& n, const int& r);


