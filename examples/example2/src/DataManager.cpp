#include "DataManager.h"
#include "DataSourceFactory.h"
#include <smartmet/newbase/NFmiMetTime.h>
#include <smartmet/newbase/NFmiPoint.h>

// ----------------------------------------------------------------------
/*!
 * \brief Construct configured data sources
 *
 * Sample settings:
 *
 *   "time":
 *   {
 *       "now": "20180101T1200", // for overriding wall clock
 *       "analysis": 48,
 *       "forecast": 120
 *   },
 *   "input":
 *   [
 *     {
 *        "name": "laps",
 *        "path": "/smartmet/data/laps/skandinavia/pinta/querydata",
 *        "type": "directory",
 *        "params": ["airtemperature","humidity"]
 *     },
 *     ....
 *   ]
 *
 */
// ----------------------------------------------------------------------

void DataManager::init(const NFmiMetTime& pWallClock,
                       const NFmiMetTime& pStartTime,
                       const NFmiMetTime& pEndTime,
                       const Json::Value& pConfig)
{
  mForecastTime = pWallClock;
  mStartTime = pStartTime;
  mEndTime = pEndTime;

  // Read data source settings

  Json::Value nulljson;

  auto json = pConfig.get("input", nulljson);
  if (json.isNull())
    throw std::runtime_error("Config variable 'input' must be set");
  if (!json.isArray())
    throw std::runtime_error(
        "Config variable 'input' must be an array of JSON objects defining data sources");

  for (const auto& source : json)
  {
    mDataSources.push_back(
        DataSourceFactory::create(mForecastTime, mStartTime, mEndTime, source, pConfig));
    if (source.get("source", nulljson) == "observations")
      mDataSources.back().is_observation = true;
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Retrieve available observations/forecasts for the coordinate
 */
// ----------------------------------------------------------------------

void DataManager::GetWeather(InputData& pData,
                             const SimulationTimes& pTimes,
                             const NFmiPoint& pLatLon) const
{
  for (const auto& source : mDataSources)
  {
    source.GetWeather(pData, pTimes, pLatLon);
//   std::cout<<pData<<std::endl;
  }
  //  std::exit(EXIT_FAILURE);
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish latest road observation available based on troad
 */
// ----------------------------------------------------------------------

boost::optional<NFmiMetTime> DataManager::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                           const std::string& variable) const
{
  boost::optional<NFmiMetTime> max_time;

  for (const auto& source : mDataSources)
  {
    if (source.is_observation)
    {
      auto tmp = source.GetLatestObsTime(pLatLon, variable);
      if (tmp)
      {
        if (!max_time || (*tmp > *max_time))
          max_time = tmp;
      }
    }
  }

  return max_time;
}
