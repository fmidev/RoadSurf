#pragma once

#include "DataSource.h"
#include <json/json.h>
#include <memory>
#include <string>


//!>Class for reading data in ascii format
//!>The format is same as with previous version
//!>of the road weather model
class AsciiSource : public DataSource
{
 public:
  ~AsciiSource() override = default;

  AsciiSource() = delete;
  AsciiSource(const AsciiSource& pOther) = delete;
  AsciiSource& operator=(const AsciiSource& pOther) = delete;
  AsciiSource(AsciiSource&& pOther) = delete;
  AsciiSource& operator=(AsciiSource&& pOther) = delete;

  AsciiSource(const NFmiMetTime& pWallClock,
              const NFmiMetTime& pStartTime,
              const NFmiMetTime& pEndTime,
              const Json::Value& pSource);

  void GetWeather(InputData& pData,
                  const SimulationTimes& pTimes,
                  const NFmiPoint& pLatLon) const final;

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const final;

 private:
  class Impl;
  std::shared_ptr<Impl> impl;
};
