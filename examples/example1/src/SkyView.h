#pragma once
#include "InputData.h"
#include "LocalParameters.h"
#include "JsonTools.h"
#include "SkyVariables.h"
#include <smartmet/macgyver/NearTree.h>
#include <smartmet/macgyver/NearTreeLatLon.h>

struct InputSettings;
//Class for reading sky view and local horizon angles
class SkyView
{
 public:
  //Reads sky view and local horizon angles
  explicit SkyView(const Json::Value& pJson);
  //Get sky view factor and local horizon angles for point, by its ID or else
  //from the nearest station
  void GetSkyVariables(LocalParameters& parameters, InputData& pData, int pointID);

 private:
  //Data map containing SkyVariables
  using DataMap = std::map<int, SkyVariables>;
  mutable DataMap mData;
  //Station locations for finding the nearest station
  using Station = Fmi::NearTreeLatLon<int>;
  using StationTree = Fmi::NearTree<Station, Fmi::NearTreeLatLonDistance<Station>>;
  mutable StationTree mStationTree;

};  // class SkyView
