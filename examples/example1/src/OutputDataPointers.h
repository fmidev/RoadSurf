#pragma once

//!> pointers to output data arrays given by modelRunner.cpp

struct OutputDataPointers
{
  int outputLen;                        //!< output data length
  const double *c_TsurfOut;             //!< pointer to output surface temperature
                                        //!< vector
  const double *c_SnowOut;              //!< pointer to output snow vector
  const double *c_WaterOut;             //!< pointer to output water vector
  const double *c_IceOut;               //!< pointer to output ice vector
  const double *c_DepositOut;           //!< pointer to output deposit vector
  const double *c_Ice2Out;              //!< pointer to output secondary ice
                                        //!< vector
};
