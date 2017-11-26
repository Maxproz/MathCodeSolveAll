#pragma once


#ifndef FUNCTIONENUMS_H
#define FUNCTIONENUMS_H




// Inclusive uses bracket notation
enum class Domain
{
	NegInfinityToPosInfinity,
	InclusiveZeroToPosInfinity,
	ExclusiveZeroToPosInfinity,
	Reciprocal,
	ExactGiven
};

enum class Range
{
	NegInfinityToPosInfinity,
	InclusiveZeroToPosInfinity,
	ExclusiveZeroToPosInfinity,
	Reciprocal,
	ExactGiven
};

enum class LineBehavior
{
	Increasing,
	Decreasing,
	Horizontal
};

enum class EndBehavior
{
	// The two behaviors for quadratic functions
	AsXGoesToPosOrNegInfinityFOfXGoesToPosInfinity,
	AsXGoesToPosOrNegInfinityFOfXGoesToNegInfinity,

	// TODO: maybe add more behaviors for other more functions later(cubic etc)

	// Exponential Func Behaviors
	// b > 1
	AsXGoesToPosInfinityFOfXGoesToPosInfinity,
	AsXGoesToNegInfinityFOfXGoesToZero,
	// 0 < b < 1
	AsXGoesToPosInfinityFOfXGoesToZero,
	AsXGoesToNegInfinityFOfXGoesToPosInfinity,
};


enum class ParabolaOpen
{
	UP,
	DOWN
};

enum class DiscontinunityType
{
	REMOVEABLE,
	JUMP,
	INFINITE,
};

void PrintDiscontinunityType(const DiscontinunityType& InType);




#endif