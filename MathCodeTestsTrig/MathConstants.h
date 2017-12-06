#pragma once


#ifndef MATHCONSTANTS_H
#define MATHCONSTANTS_H

#include "FunctionEnums.h"

#include <utility>

using std::pair;


#define NEGINFINITY   ((float)(_HUGE_ENUF * _HUGE_ENUF) *(-1))

const double Eulere = 2.718282;
const double my_gravityfeet(32);  // f/s^2

typedef pair<double, double> Point;
typedef std::tuple<float, float, IntervalType> Interval;
typedef std::pair<float, float> FuncEndBehavior;


#endif