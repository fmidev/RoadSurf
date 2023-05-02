#include "MeteorologyTools.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
//-----------------------------------------------------------------------
/*!
 * \brief Convert dewpoint temperature to relative humidity or vice verse
 */
// ----------------------------------------------------------------------

double CalcTdewOrRH(double t, double tdew, double rh)
{
  //Set some parameters
  const double Alphaw = 17.269;  // over water
  const double Alphai = 21.875;  // over ice
  const double Betaw = 237.3;    // over water
  const double Betai = 265.5;    // over ice
  const double AFact = 0.61078;  // e in kPa
  double Alpha = 0.0;
  double Beta = 0.0;

  //Checi if temperature is over 0.0
  if (t >= 0.0)
  {
    Alpha = Alphaw;
    Beta = Betaw;
  }
  else
  {
    Alpha = Alphai;
    Beta = Betai;
  }
  //Saturation temperature
  double EsatT = AFact * exp(Alpha * t / (t + Beta));
  //if tdew is not missing, calculate relative humidity
  if (!std::isnan(tdew) and tdew > -1000)
  {
    double EsatTD = AFact * exp(Alpha * tdew / (tdew + Beta));
    return std::min((EsatTD / EsatT) * 100.0, 100.0);
  }
  //if relative humidity is not missing, calculate dew point temperature
  if (!std::isnan(rh) and rh > -1)
  {
    double Epr = 0.01 * rh * EsatT;
    double XX = log(Epr / AFact);
    return Beta * XX / (Alpha - XX);
  }

  return std::numeric_limits<double>::quiet_NaN();
}

