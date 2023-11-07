#pragma once

#include <macgyver/DateTime.h>
#include <smartmet/newbase/NFmiMetTime.h>
#include <vector>

// ----------------------------------------------------------------------
/*!
 * \brief Storage for common simulation timesteps
 */
// ----------------------------------------------------------------------

struct SimulationTime
{
  explicit SimulationTime(const Fmi::DateTime& pTime) : pt(pTime), t(pTime)
  {
    const auto ymd = pt.date().year_month_day();
    const auto tday = pt.time_of_day();

    year = ymd.year;
    month = ymd.month;
    day = ymd.day;
    hour = tday.hours();
    minute = tday.minutes();
    second = tday.seconds();
  }

  Fmi::DateTime pt;
  NFmiMetTime t = NFmiMetTime::gMissingTime;
  int year = 0;
  int month = 0;
  int day = 0;
  int hour = 0;
  int minute = 0;
  int second = 0;
};

using SimulationTimes = std::vector<SimulationTime>;
