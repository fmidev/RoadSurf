#pragma once

#include <string>

double CalcTdewOrRH(double t, double tdew, double RH);
float ConvertPrecPhase(const std::string& precConversion, float value);
