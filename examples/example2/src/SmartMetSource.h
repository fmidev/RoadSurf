#pragma once

#include "DataSource.h"
#include <json/json.h>
#include <memory>
#include <string>

class SmartMetSource : public DataSource
{
 public:
  ~SmartMetSource() override = default;

  SmartMetSource() = delete;
  SmartMetSource(const SmartMetSource& pOther) = delete;
  SmartMetSource& operator=(const SmartMetSource& pOther) = delete;
  SmartMetSource(SmartMetSource&& pOther) = delete;
  SmartMetSource& operator=(SmartMetSource&& pOther) = delete;

  SmartMetSource(const NFmiMetTime& pWallClock,
                 const NFmiMetTime& pStartTime,
                 const NFmiMetTime& pEndTime,
                 const Json::Value& pSource,
                 const Json::Value& pConfig);

  void GetWeather(InputData& pData,
                  const SimulationTimes& pTimes,
                  const NFmiPoint& pLatLon) const final;

  boost::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const final;

 private:
  class Impl;
  std::unique_ptr<Impl> impl;
};
