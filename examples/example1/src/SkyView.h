#pragma once
#include "InputData.h"
#include "LocalParameters.h"
#include "JsonTools.h"
#include "SkyVariables.h"

struct InputModelSettings;
//Class for reading sky view and local horizon angles
class SkyView
{
 public:
  //Reads sky view and local horizon angles
  explicit SkyView(const Json::Value& pJson);
  //Get sky view factor and local horizon angles for point
  void GetSkyVariables(LocalParameters& parameters, InputData& pData, int pointID);

 private:
  //Data map containing SkyVariables
  using DataMap = std::map<int, SkyVariables>;
  mutable DataMap mData;

};  // class SkyView
