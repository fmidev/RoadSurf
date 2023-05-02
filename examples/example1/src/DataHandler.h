#pragma once

#include "GenericSource.h"
#include "InputData.h"
#include "LocalParameters.h"
#include <json/json.h>
#include <time.h>

//!> Class for managing input data
//It is possible to read in different datatyps as long as
//they have implementation corresponding GenericSource class
class DataHandler
{
 public:
  ~DataHandler(); //Destructor
  DataHandler() = default;

  DataHandler(const DataHandler& pOther) = delete;
  DataHandler& operator=(const DataHandler& pOther) = delete;
  DataHandler(DataHandler&& pOther) = delete;
  DataHandler& operator=(DataHandler&& pOther) = delete;

  //Creates DataHanlder and reads in input data
  void init(const time_t& pWallClock,
            const time_t& pStartTime,
            const time_t& pEndTime,
            const Json::Value& pConfig,
            const int SimLen,
            const int DTSecs);

  //Get data for one point
  void GetWeather(InputData& pData, int pointID) const;
  //Get vecotr of station locations
  std::vector<LocalParameters> GetLocations() const;
  //Get Vector of pointIDS
  std::vector<int> GetPointIDs() const;
  //Get Index of last observation
  int GetLatestObsIndex(int pointID) const;

 private:
  std::vector<GenericSource*> mDataSources;
  time_t mForecastTime;
  time_t mStartTime;
  time_t mEndTime;

};  // class DataHandler
