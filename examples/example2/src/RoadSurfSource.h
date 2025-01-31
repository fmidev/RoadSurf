#pragma once

#include "DataSource.h"
#include <json/json.h>
#include <memory>
#include <string>

class RoadSurfSource : public DataSource
{
 public:
  ~RoadSurfSource() override;

  RoadSurfSource() = delete;
  RoadSurfSource(const RoadSurfSource& pOther) = delete;
  RoadSurfSource& operator=(const RoadSurfSource& pOther) = delete;
  RoadSurfSource(RoadSurfSource&& pOther) = delete;
  RoadSurfSource& operator=(RoadSurfSource&& pOther) = delete;

  RoadSurfSource(const NFmiMetTime& pWallClock,
                 const NFmiMetTime& pStartTime,
                 const NFmiMetTime& pEndTime,
                 const Json::Value& pSource,
                 const Json::Value& pConfig);

  void GetWeather(InputData& pData,
                  const SimulationTimes& pTimes,
                  const NFmiPoint& pLatLon) const final;

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const final;

 private:
  class Impl;
  std::unique_ptr<Impl> impl;
};
