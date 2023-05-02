#pragma once
//STructure for storing sky view factor
//and local horizon angles
struct SkyVariables
{
  double sky_view = 1.0;
  std::vector<double> local_horizons;
};
