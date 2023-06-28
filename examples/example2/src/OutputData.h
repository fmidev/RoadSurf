#pragma once
#include "OutputPointers.h"
#include <iostream>
#include <vector>

//!> Arrays for model output data

struct OutputData
{
  explicit OutputData(int pSize);

  std::vector<double> TsurfOut;    //!< Road surface temperature (C)
  std::vector<double> SnowOut;     //!< snow storage (water
                                   //!< equivalent mm)
  std::vector<double> WaterOut;    //!< water storage (mm)
  std::vector<double> IceOut;      //!< ice storage (water equivalent mm)
  std::vector<double> DepositOut;  //!< deposit storage (water
                                   //!< equivalent mm)
  std::vector<double> Ice2Out;     //!< secondary ice storage
                                   //!< (water equivalent mm)
  OutputPointers pointers() const;
};

std::ostream& operator<<(std::ostream& out, const OutputData& data);
