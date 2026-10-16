#pragma once
#include "InputData.h"
#include "LocalParameters.h"
#include "JsonTools.h"
#include "SkyVariables.h"
#include <smartmet/macgyver/NearTree.h>
#include <smartmet/macgyver/NearTreeLatLon.h>
#include <smartmet/newbase/NFmiPoint.h>
#include <smartmet/newbase/NFmiStation.h>

struct InputSettings;

class SkyView
{
 public:
  explicit SkyView(const Json::Value& pJson);
  void GetSkyVariables(LocalParameters& parameters, InputData& pData, const NFmiPoint& pLatLon);

 private:
  void add_local_horizons(int id,
                          double lat,
                          double lon,
                          const std::vector<double>& local_horizons);

  using Station = Fmi::NearTreeLatLon<int>;
  using StationTree = Fmi::NearTree<Station, Fmi::NearTreeLatLonDistance<Station>>;
  using DataMap = std::map<int, SkyVariables>;
  mutable StationTree mStationTree;
  mutable DataMap mData;

};  // class SkyView
