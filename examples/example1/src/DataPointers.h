#pragma once

extern "C"
{
  //!> pointers to input data arrays given

  struct DataPointers
  {
    int inputLen;        //!< input data length

    const double *c_tair;            //!< pointer to input air temperature vector
    const double *c_tdew;            //!< pointer to input dew point temperature vector
    const double *c_VZ;              //!< pointer to input wind speed vector
    const double *c_Rhz;             //!< pointer to input relative humidity
    const double *c_prec;            //!< pointer to input precipitation vector
    const double *c_SW;              //!< pointer to input short wave radiation vector
    const double *c_LW;              //!< pointer to input long wave vector
    const double *c_SW_dir;          //!< pointer to input direct sw radiation
    const double *c_LW_net;          //!< oiubter ti ínput net lw radiation
    const double *c_TSurfObs;        //!< pointer to input tsurf obs vector
    const int *c_PrecPhase;          //!< pointer to input prec phase vector
    const double *c_local_horizons;  //!< pointer to local horizon angle vector
    const double *c_Depth;           //!< pointer to input depth vector
    const int *c_year;               //!< pointer to input simulated year vector
    const int *c_month;              //!< pointer to input simulated month vector
    const int *c_day;                //!< pointer to input simulated day vector
    const int *c_hour;               //!< pointer to input simulated hour vector
    const int *c_minute;             //!< pointer to input simulated minute vector
    const int *c_second;             //!< pointer to input simulated second vector
  };
}
