#pragma once

#include <smartmet/newbase/NFmiMetTime.h>
#include <limits>

/*
 *! \brief A C++ data holder for the input weather at some point & time
 */

struct Weather
{
  NFmiMetTime date = NFmiMetTime::gMissingTime;
  double troad = std::numeric_limits<double>::quiet_NaN();
  double t2m = std::numeric_limits<double>::quiet_NaN();
  double tdew2m = std::numeric_limits<double>::quiet_NaN();
  double rh2m = std::numeric_limits<double>::quiet_NaN();
  double vz = std::numeric_limits<double>::quiet_NaN();
  double lwdn = std::numeric_limits<double>::quiet_NaN();
  double swdn = std::numeric_limits<double>::quiet_NaN();
  double lwnet = std::numeric_limits<double>::quiet_NaN();
  double swdir = std::numeric_limits<double>::quiet_NaN();
  double simuprec = std::numeric_limits<double>::quiet_NaN();
  // TODO? use int instead
  double phase = std::numeric_limits<double>::quiet_NaN();
  double depth = std::numeric_limits<double>::quiet_NaN();

  Weather() = default;

  explicit Weather(const NFmiMetTime& d) : date(d) {}

  Weather(const NFmiMetTime& d,
          double tr,
          double t,
          double tdew,
          double rh,
          double v,
          double lw,
          double sw,
          double lwn,
          double swd,
          double simu,
          double ph,
          double dth)
      : date(d),
        troad(tr),
        t2m(t),
        tdew2m(tdew),
        rh2m(rh),
        vz(v),
        lwdn(lw),
        swdn(sw),
        lwnet(lwn),
        swdir(swd),
        simuprec(simu),
        phase(ph),
        depth(dth)
  {
  }
};

const double MISSING = std::numeric_limits<double>::quiet_NaN();

inline bool is_missing(double pValue)
{
  return (std::isnan(pValue) || pValue < -9000);
}
