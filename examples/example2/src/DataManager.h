#pragma once

#include "DataSource.h"
#include "InputData.h"
#include "Weather.h"
#include <optional>
#include <boost/ptr_container/ptr_vector.hpp>
#include <json/json.h>
#include <smartmet/newbase/NFmiMetTime.h>

class NFmiPoint;
//!> Class for managing input data
class DataManager
{
 public:
  ~DataManager() = default;
  DataManager() = default;

  DataManager(const DataManager& pOther) = delete;
  DataManager& operator=(const DataManager& pOther) = delete;
  DataManager(DataManager&& pOther) = delete;
  DataManager& operator=(DataManager&& pOther) = delete;

  void init(const NFmiMetTime& pWallClock,
            const NFmiMetTime& pStartTime,
            const NFmiMetTime& pEndTime,
            const Json::Value& pConfig);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon) const;
  double MaxShortWaveRadiation(const NFmiPoint& pLatLon) const;
  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const;

 private:
  boost::ptr_vector<DataSource> mDataSources;
  NFmiMetTime mForecastTime;
  NFmiMetTime mStartTime;
  NFmiMetTime mEndTime;

};  // class DataManager
