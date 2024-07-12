#pragma once

#include "DataSource.h"
#include <json/json.h>
#include <smartmet/newbase/NFmiParameterName.h>
#include <memory>
#include <string>

class QueryDataSource : public DataSource
{
 public:
  ~QueryDataSource() override = default;

  QueryDataSource() = delete;
  QueryDataSource(const QueryDataSource& pOther) = delete;
  QueryDataSource& operator=(const QueryDataSource& pOther) = delete;
  QueryDataSource(QueryDataSource&& pOther) = delete;
  QueryDataSource& operator=(QueryDataSource&& pOther) = delete;

  QueryDataSource(const NFmiMetTime& pWallClock,
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
  std::unique_ptr<Impl> impl;
};
