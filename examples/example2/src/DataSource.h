#pragma once

#include "InputData.h"
#include "SimulationTime.h"
#include "Weather.h"
#include <optional>
#include <json/json.h>

class NFmiPoint;
class NFmiMetTime;

class DataSource
{
 public:
  DataSource() = default;

  DataSource(const DataSource& pOther) = delete;
  DataSource& operator=(const DataSource& pOther) = delete;
  DataSource(DataSource&& pOther) = delete;
  DataSource& operator=(DataSource&& pOther) = delete;

  virtual ~DataSource() = default;

  virtual void GetWeather(InputData& pData,
                          const SimulationTimes& pTimes,
                          const NFmiPoint& pLatLon) const = 0;

  virtual std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                        const std::string& variable) const = 0;
  bool is_observation = false;

 private:
};  // class DataSource
