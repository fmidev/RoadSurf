#pragma once
#include <roadsurf/Constants.h>
#include <json/json.h>
#include <cmath>

struct InputSettings;

//!> input parameters given by modelRunner.cpp

struct LocalParameters
{
  LocalParameters() = delete;

  explicit LocalParameters(const InputSettings& pSettings);

  double tair_relax = -9999.0;     //!< tair for relaxation
  double VZ_relax = -9999.0;       //!< wind speed for relaxation
  double RH_relax = -9999.0;       //!< relative humidity for relaxation
  int couplingIndexI = -9999;      //!< index in input data where coupling starts
  double couplingTsurf = -9999.0;  //!< surface temperature value used in coupling
  double lat = -9999.0;            //!< latitude
  double lon = -9999.0;            //!< longitude
  double sky_view = 1.0;           //!< Sky view factor
  int InitLenI = 0;                //!< The length of initialization period,
                                   //!< input timestepsa
};
