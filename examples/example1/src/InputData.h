#pragma once
#include "DataPointers.h"
#include <ctime>
#include <iostream>
#include <vector>
extern "C"
{
  //!> input data vectors

  struct InputData
  {
    InputData() = delete;
    explicit InputData(int pSize);

    std::vector<double> tair;      //!< air temperature (C)
    std::vector<double> tdew;      //!< dew point temperature (C)
    std::vector<double> VZ;        //!< wind speed (m/s)
    std::vector<double> Rhz;       //!< relative humidity (%)
    std::vector<double> prec;      //!< precipitation (mm/h)
    std::vector<double> SW;        //!< incoming short wave
                                   //!< radiation (W/m2)
    std::vector<double> LW;        //!< incoming long wave radiation (W/m2)
    std::vector<double> SW_dir;    //!< direct short wave radiation (W/m2)
    std::vector<double> LW_net;    //!< net long wave radiation (W/m2)
    std::vector<double> TSurfObs;  //!< surface temperature  (C)
    std::vector<int> PrecPhase;    //!< precipitation phase (Hail = 6;
                                   //!< FreezingRain = 5;
                                   //!< FreezingDrizzle = 4; Snow = 3;
                                   //!< Sleet = 2;
                                   //!< Rain = 1; Drizzle = 0;)

    std::vector<double> local_horizons;  //!< local horizon angles
    std::vector<double> Depth;           //!< Depth to calculate surface temperatre

    std::vector<int> year;    //!< Simulated year
    std::vector<int> month;   //!< Simulated month
    std::vector<int> day;     //!< Simulated day
    std::vector<int> hour;    //!< Simulated hour
    std::vector<int> minute;  //!< Simulated minute
    std::vector<int> second;  //!< Simulated second
    // Return pointers to vectors

    DataPointers pointers() const;
  };
}

std::ostream& operator<<(std::ostream& out, const InputData& data);
