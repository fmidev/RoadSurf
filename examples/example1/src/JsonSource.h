#pragma once

#include "GenericSource.h"
#include "LocalParameters.h"
#include "InputSettings.h"
#include <json/json.h>
#include <memory>
#include <string>

//Class for reading JSon data, Is based on Generic Source class
class JsonSource : public GenericSource
{
 public:
  ~JsonSource() override = default;  //Destructor

  JsonSource() = delete;
  JsonSource(const JsonSource& pOther) = delete;
  JsonSource& operator=(const JsonSource& pOther) = delete;
  JsonSource(JsonSource&& pOther) = delete;
  JsonSource& operator=(JsonSource&& pOther) = delete;

  //Crates JsonSource and reads input data
  JsonSource(const time_t& pWallClock,
                 const time_t& pStartTime,
                 const time_t& pEndTime,
                 const Json::Value& pSource,
                 const Json::Value& pConfig,
                 const int SimLen,
                 const int DTSecs);
  //Get data for given point
  void GetWeather(InputData& pData,
                  int pointID) const final;
  //Get vector containgin LocalParameters 
  std::vector<LocalParameters> GetLocations() const final;
  //Get wector of pointIDs
  std::vector<int> GetPointIDs() const final;
  //Get index of last observation
  int GetLatestObsIndex(int statPos) const final;

 private:
  class Impl;
  std::unique_ptr<Impl> impl;
};
